#include "haero/processes/mam_calcsize_hostcxx_process.hpp"

#include <algorithm>
#include <cmath>
#include <ekat/ekat.hpp>
#include <ekat/ekat_pack_math.hpp>
#include <ekat/logging/ekat_logger.hpp>
#include <haero/constants.hpp>
#include <limits>

namespace {

using haero::Real;

template <typename Iter>
static inline auto fmt_vec_1d(Iter it, Iter end) -> std::string {
  std::stringstream ss;
  ss << "[";
  while (it++ != end) ss << *it << ",";
  ss << "]";
  return ss.str();
}

template <typename Container>
static inline auto fmt_vec_1d(Container c) -> std::string {
  std::stringstream ss;
  ss << "[";
  for (const auto &e : c) ss << e << ",";
  ss << "]";
  return ss.str();
}

template <typename Container>
static inline auto fmt_vec_2d(Container c) -> std::string {
  std::stringstream ss;
  ss << "[";
  for (const auto &e : c) ss << fmt_vec_1d(e) << ",";
  ss << "]";
  return ss.str();
}

template <typename ViewType, typename T = typename ViewType::value_type>
static inline auto view_to_vector_2d(ViewType view) {
  static_assert(
      std::is_same<Kokkos::HostSpace, typename ViewType::memory_space>::value,
      "CalcsizeHostCXX only runs in host space");
  std::vector<std::vector<T>> vec;
  vec.resize(view.extent(0));
  for (int i = 0; i < view.extent(0); i++) {
    vec[i].resize(view.extent(1));
    for (int j = 0; j < view.extent(1); j++) vec[i][j] = view(i, j);
  }
  return vec;
}

template <typename ViewType, typename T = typename ViewType::value_type>
static inline auto view_to_vector_1d(ViewType view) {
  static_assert(
      std::is_same<Kokkos::HostSpace, typename ViewType::memory_space>::value,
      "CalcsizeHostCXX only runs in host space");
  std::vector<T> vec;
  vec.resize(view.extent_0());
  for (int i = 0; i < view.extent_0(); i++) vec[i] = view(i);
  return vec;
}

template <typename T = Real>
static inline auto zero_vec(std::size_t sz) {
  return std::vector<T>(sz, T{0});
}

template <typename T = Real>
static inline auto zero_vec_2d(std::size_t s0, std::size_t s1) {
  return std::vector<std::vector<T>>(s0, std::vector<T>(s1, T{0}));
}

}  // namespace

namespace haero {

MAMCalcsizeHostCXXProcess::MAMCalcsizeHostCXXProcess()
    : DeviceAerosolProcess<MAMCalcsizeHostCXXProcess>(
          CalcsizeProcess, "MAMCalcsizeHostCXXProcess") {
  logger->set_pattern("MAMCalcsizeHostCXXProcess[%^%l%$]: %v");
}
//------------------------------------------------------------------------
//                                Accessors
//------------------------------------------------------------------------

void MAMCalcsizeHostCXXProcess::init_(
    const ModalAerosolConfig &modal_aerosol_config) {
  logger->debug("entering init_");

  nmodes = modal_aerosol_config.num_modes();
  logger->debug("nmodes=modal_aerosol_config.num_modes={}",
                modal_aerosol_config.num_modes());

  num_populations = modal_aerosol_config.num_aerosol_populations;
  logger->debug(
      "num_populations=modal_aerosol_config.num_aerosol_populations={}",
      modal_aerosol_config.num_aerosol_populations);

  num_mode_species.resize(nmodes);
  for (int i = 0; i < nmodes; i++)
    num_mode_species[i] =
        modal_aerosol_config.aerosol_species_for_mode(i).size();
  logger->debug("num mode species={}", fmt_vec_1d(num_mode_species));

  max_nspec = std::accumulate(num_mode_species.begin(), num_mode_species.end(),
                              num_mode_species.front(),
                              [](int a, int b) { return std::max(a, b); });
  logger->debug("max nspec={}", max_nspec);

  spec_density.resize(max_nspec * nmodes);
  const auto all_species = modal_aerosol_config.aerosol_species;
  logger->debug("modal_aerosol_config.aerosol_species.size={}",
                modal_aerosol_config.aerosol_species.size());
  for (int i = 0; i < all_species.size(); i++)
    logger->debug("spec={} has density={}", i, all_species[i].density);
  std::transform(all_species.begin(), all_species.end(), spec_density.begin(),
                 [](const auto &spec) -> int { return spec.density; });
  logger->debug("spec_density={}", fmt_vec_1d(spec_density));

  population_offsets.resize(nmodes);
  for (int i = 0; i < nmodes; i++)
    population_offsets[i] = modal_aerosol_config.population_index(i, 0);
  logger->debug("population offsets={}", fmt_vec_1d(population_offsets));

  v2nmin_nmodes = zero_vec(nmodes);
  v2nmax_nmodes = zero_vec(nmodes);
  dgnmin_nmodes = zero_vec(nmodes);
  dgnmax_nmodes = zero_vec(nmodes);
  common_factor = zero_vec(nmodes);

  for (int i = 0; i < nmodes; i++) {
    logger->debug("setting volume/num ratios and diameters mode={}", i);
    const auto &mode = modal_aerosol_config.aerosol_modes[i];
    using T = decltype(v2nmin_nmodes)::value_type;
    v2nmin_nmodes[i] = mode.min_vol_to_num_ratio<T>();
    v2nmax_nmodes[i] = mode.max_vol_to_num_ratio<T>();
    dgnmin_nmodes[i] = mode.min_diameter;
    dgnmax_nmodes[i] = mode.max_diameter;
    common_factor[i] =
        std::exp(4.5 * std::log(std::pow(mode.mean_std_dev, 2.0))) *
        Constants::pi_sixth;
  }

  aitken_idx = modal_aerosol_config.aerosol_mode_index("aitken");
  accum_idx = modal_aerosol_config.aerosol_mode_index("accum");
  logger->debug("leaving init_");
}

KOKKOS_FUNCTION
void MAMCalcsizeHostCXXProcess::run_(const TeamType &team, Real t, Real dt,
                                     const Prognostics &prognostics,
                                     const Atmosphere &atmosphere,
                                     const Diagnostics &diagnostics,
                                     const Tendencies &tendencies) const {
  logger->debug("entering run_");

  nlevels = prognostics.num_levels();

  std::size_t num_vert_packs = nlevels / HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < nlevels) {
    num_vert_packs++;
  }
  logger->debug("num levels={}, num vertical packs={}", nlevels,
                num_vert_packs);

  // interstitial mass and number mixing ratios
  const auto q_i = view_to_vector_2d(prognostics.interstitial_aerosols);
  logger->debug("prognostics.interstitial_aerosols.size={}", q_i.size());
  const auto n_i = view_to_vector_2d(prognostics.interstitial_num_mix_ratios);
  logger->debug("prognostics.interstitial_num_mix_ratios.size={}", n_i.size());

  // cloud-borne mass and number mixing ratios
  const auto q_c = view_to_vector_2d(prognostics.cloud_aerosols);
  logger->debug("prognostics.cloud_aerosols.size={}", q_c.size());
  const auto n_c = view_to_vector_2d(prognostics.cloud_num_mix_ratios);
  logger->debug("prognostics.cloud_num_mix_ratios.size={}", n_c.size());

  // tendencies for interstitial number mixing ratios
  const auto dnidt = view_to_vector_2d(tendencies.interstitial_num_mix_ratios);
  logger->debug("tendencies.interstitial_num_mix_ratios.size={}", n_c.size());

  // tendencies for cloud-borne number mixing ratios
  const auto dncdt = view_to_vector_2d(tendencies.cloud_num_mix_ratios);
  logger->debug("tendencies.cloud_num_mix_ratios.size={}", dncdt.size());

  auto dryvol_a = std::vector<PackType>(num_vert_packs, PackType{0});
  auto dryvol_c = std::vector<PackType>(num_vert_packs, PackType{0});

  auto dgncur_a = zero_vec_2d<PackType>(nmodes, num_vert_packs);
  auto dgncur_c = zero_vec_2d<PackType>(nmodes, num_vert_packs);
  auto v2ncur_a = zero_vec_2d<PackType>(nmodes, num_vert_packs);
  auto v2ncur_c = zero_vec_2d<PackType>(nmodes, num_vert_packs);

  // specie density array for each mode [kg/m3]
  auto density = std::vector<Real>(max_nspec, 0.0);

  // Loop through each mode and find particle diameter
  for (int imode = 0; imode < nmodes; imode++) {
    logger->debug("Finding particle diameter for mode={}", imode);

    set_initial_sz_and_volumes_(imode, top_level, nlevels, dgncur_a, v2ncur_a,
                                dryvol_a, num_vert_packs);
    set_initial_sz_and_volumes_(imode, top_level, nlevels, dgncur_c, v2ncur_c,
                                dryvol_c, num_vert_packs);

    // species starting index in the population (q_i and q_c) arrays for a mode
    const auto start_spec_idx = population_offsets[imode];

    // end index of species for all modes expect the last mode
    const auto end_spec_idx = ((1 + imode) == nmodes)
                                  ? num_populations
                                  : population_offsets[imode + 1] - 1;

    const auto nspec = num_mode_species[imode];

    logger->debug("species start idx={},end idx={},nspec={}", start_spec_idx,
                  end_spec_idx, nspec);

    // capture densities for each specie in this mode
    std::fill_n(density.begin(), max_nspec,
                std::numeric_limits<decltype(density)::value_type>::max());
    std::copy_n(spec_density.begin() + population_offsets[imode], nspec,
                density.begin());
    // logger->debug("densities for species in mode={}",
    // fmt_vec_1d(density.begin(), density.begin() + nspec));

    compute_dry_volume(imode, top_level, nlevels, start_spec_idx, end_spec_idx,
                       density, q_i, q_c, dryvol_a, dryvol_c, num_vert_packs);
  }
  logger->debug("leaving run_");
}

void MAMCalcsizeHostCXXProcess::set_initial_sz_and_volumes_(
    const int imode, const int top_lev, const int nlevs,
    std::vector<std::vector<PackType>> &dgncur,
    std::vector<std::vector<PackType>> &v2ncur, std::vector<PackType> &dryvol,
    const std::size_t num_vert_packs) const {
  logger->debug("set_initial_sz_and_volumes_:imode={},top_lev={},nlevs={}",
                imode, top_lev, nlevs);
  logger->debug(
      "set_initial_sz_and_volumes_:dgncur.size={},v2ncur.size={},dryvol.size={"
      "}",
      dgncur.size(), v2ncur.size(), dryvol.size());
  logger->debug(
      "set_initial_sz_and_volumes_:dgncur[0].size={},v2ncur[0].size={}",
      dgncur[0].size(), v2ncur[0].size());
  for (int pack_idx = 0; pack_idx < num_vert_packs; pack_idx++) {
    dgncur[imode][pack_idx] = PackType::scalar{0.0};
    v2ncur[imode][pack_idx] = PackType::scalar{0.0};
    dryvol[pack_idx] = PackType::scalar{0.0};
  }
}

void MAMCalcsizeHostCXXProcess::compute_dry_volume(
    const int imode, const int top_lev, const int nlevs, const int s_spec_ind,
    const int e_spec_ind, const std::vector<Real> &density,
    const std::vector<std::vector<PackType>> &q_i,
    const std::vector<std::vector<PackType>> &q_c,
    std::vector<PackType> &dryvol_a, std::vector<PackType> &dryvol_c,
    const std::size_t num_vert_packs) const {
  using namespace ekat;
  EKAT_REQUIRE_MSG(top_lev == 0, "top level must be zero");
  logger->debug("compute_dry_volume:imode={}", imode);
  for (int ispec = s_spec_ind; ispec < e_spec_ind; ispec++) {
    const auto density_ind = ispec - s_spec_ind;
    const PackType::scalar inv_density = 1.0 / density[density_ind];
    logger->debug("compute_dry_volume:spec={},density_ind={}", imode, ispec,
                  density_ind);
    for (int pack_idx = 0; pack_idx < num_vert_packs; pack_idx++) {
      dryvol_a[pack_idx] += max(0.0, q_i[ispec][pack_idx]) * inv_density;
      dryvol_c[pack_idx] += max(0.0, q_i[ispec][pack_idx]) * inv_density;
    }
  }
}

}  // namespace haero
