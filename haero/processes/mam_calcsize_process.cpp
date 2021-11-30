#include "haero/processes/mam_calcsize_process.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>

namespace haero {

MAMCalcsizeProcess::MAMCalcsizeProcess()
    : DeviceAerosolProcess<MAMCalcsizeProcess>("MAMCalcsizeProcess") {}
//------------------------------------------------------------------------
//                                Accessors
//------------------------------------------------------------------------

void MAMCalcsizeProcess::init_(const ModalAerosolConfig &modal_aerosol_config) {
  nmodes = modal_aerosol_config.num_modes();

  num_populations = modal_aerosol_config.num_aerosol_populations;

  Kokkos::resize(num_mode_species, nmodes);
  {
    auto nm_spec = Kokkos::create_mirror_view(num_mode_species);
    for (int i = 0; i < nmodes; i++)
      nm_spec[i] = modal_aerosol_config.aerosol_species_for_mode(i).size();

    max_nspec = nm_spec(0);

    // TODO: this should be a parallel reduce
    for (int i = 0; i < nm_spec.extent(0); i++)
      max_nspec = std::max(max_nspec, nm_spec(i));

    Kokkos::deep_copy(num_mode_species, nm_spec);
  }

  Kokkos::resize(spec_density, max_nspec * nmodes);
  {
    const auto h_spec_density = Kokkos::create_mirror_view(spec_density);
    const auto all_species = modal_aerosol_config.aerosol_species;
    for (int i = 0; i < all_species.size(); i++)
      h_spec_density(i) = all_species[i].density;
  }

  Kokkos::resize(population_offsets, nmodes);
  {
    const auto pop_offsets = Kokkos::create_mirror_view(population_offsets);
    for (int i = 0; i < nmodes; i++)
      pop_offsets[i] = modal_aerosol_config.population_index(i, 0);
    Kokkos::deep_copy(population_offsets, pop_offsets);
  }

  Kokkos::resize(v2nmin_nmodes, nmodes);
  Kokkos::resize(v2nmax_nmodes, nmodes);
  Kokkos::resize(dgnmin_nmodes, nmodes);
  Kokkos::resize(dgnmax_nmodes, nmodes);
  Kokkos::resize(common_factor_nmodes, nmodes);

  {
    const auto h_v2nmin_nmodes = Kokkos::create_mirror_view(v2nmin_nmodes);
    const auto h_v2nmax_nmodes = Kokkos::create_mirror_view(v2nmax_nmodes);
    const auto h_dgnmin_nmodes = Kokkos::create_mirror_view(dgnmin_nmodes);
    const auto h_dgnmax_nmodes = Kokkos::create_mirror_view(dgnmax_nmodes);
    const auto h_common_factor_nmodes =
        Kokkos::create_mirror_view(common_factor_nmodes);

    for (int i = 0; i < nmodes; i++) {
      const auto &mode = modal_aerosol_config.aerosol_modes[i];
      using T = decltype(v2nmin_nmodes)::value_type;
      h_v2nmin_nmodes[i] = mode.min_vol_to_num_ratio<T>();
      h_v2nmax_nmodes[i] = mode.max_vol_to_num_ratio<T>();
      h_dgnmin_nmodes[i] = mode.min_diameter;
      h_dgnmax_nmodes[i] = mode.max_diameter;
      h_common_factor_nmodes[i] =
          std::exp(4.5 * std::log(std::pow(mode.mean_std_dev, 2.0))) *
          Constants::pi_sixth;
    }

    Kokkos::deep_copy(v2nmin_nmodes, h_v2nmin_nmodes);
    Kokkos::deep_copy(v2nmax_nmodes, h_v2nmax_nmodes);
    Kokkos::deep_copy(dgnmin_nmodes, h_dgnmin_nmodes);
    Kokkos::deep_copy(dgnmax_nmodes, h_dgnmax_nmodes);
    Kokkos::deep_copy(common_factor_nmodes, h_common_factor_nmodes);
  }

  aitken_idx = modal_aerosol_config.aerosol_mode_index("aitken");
  accum_idx = modal_aerosol_config.aerosol_mode_index("accum");
}

KOKKOS_FUNCTION
void MAMCalcsizeProcess::run_(const TeamType &team, Real t, Real dt,
                              const Prognostics &prognostics,
                              const Atmosphere &atmosphere,
                              const Diagnostics &diagnostics,
                              const Tendencies &tendencies) const {
  const int nlevels = prognostics.num_levels();

  std::size_t num_vert_packs = nlevels / HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < nlevels) {
    num_vert_packs++;
  }

  // interstitial mass and number mixing ratios
  const auto q_i = prognostics.interstitial_aerosols;
  const auto n_i = prognostics.interstitial_num_mix_ratios;

  // cloud-borne mass and number mixing ratios
  const auto q_c = prognostics.cloud_aerosols;
  const auto n_c = prognostics.cloud_num_mix_ratios;

  // tendencies for interstitial number mixing ratios
  const auto dnidt = tendencies.interstitial_num_mix_ratios;

  // tendencies for cloud-borne number mixing ratios
  const auto dncdt = tendencies.cloud_num_mix_ratios;

  ColumnView dryvol_a("dryvol_a", num_vert_packs);
  ColumnView dryvol_c("dryvol_c", num_vert_packs);

  SpeciesColumnView dgncur_a("dgncur_a", nmodes, num_vert_packs);
  SpeciesColumnView dgncur_c("dgncur_c", nmodes, num_vert_packs);
  SpeciesColumnView v2ncur_a("v2ncur_a", nmodes, num_vert_packs);
  SpeciesColumnView v2ncur_c("v2ncur_c", nmodes, num_vert_packs);

  DeviceType::view_1d<Real> density("density", max_nspec);
  for (int imode = 0; imode < nmodes; imode++) {
    set_initial_sz_and_volumes_(imode, top_level, nlevels, dgncur_a, v2ncur_a,
                                dryvol_a, num_vert_packs);
    set_initial_sz_and_volumes_(imode, top_level, nlevels, dgncur_c, v2ncur_c,
                                dryvol_c, num_vert_packs);

    // species starting index in the population (q_i and q_c) arrays for a mode
    const auto start_spec_idx = population_offsets(imode);

    // end index of species for all modes expect the last mode
    const auto end_spec_idx = ((1 + imode) == nmodes)
                                  ? num_populations
                                  : population_offsets(imode + 1) - 1;

    const auto nspec = num_mode_species(imode);

    Kokkos::parallel_for(
        max_nspec, KOKKOS_LAMBDA(int i) {
          density(i) =
              std::numeric_limits<decltype(density)::value_type>::max();
        });
    Kokkos::parallel_for(
        nspec, KOKKOS_LAMBDA(int i) {
          density(i) = spec_density(population_offsets(imode) + i);
        });

    compute_dry_volume(imode, top_level, nlevels, start_spec_idx, end_spec_idx,
                       density, q_i, q_c, dryvol_a, dryvol_c, num_vert_packs);
  }
}

void MAMCalcsizeProcess::compute_dry_volume(
    const int imode, const int top_lev, const int nlevs, const int s_spec_ind,
    const int e_spec_ind, const DeviceType::view_1d<Real> &density,
    const SpeciesColumnView q_i, const SpeciesColumnView q_c,
    ColumnView dryvol_a, ColumnView dryvol_c,
    const std::size_t num_vert_packs) const {
  using namespace ekat;
  EKAT_REQUIRE_MSG(top_lev == 0, "top level must be zero");
  Kokkos::parallel_for(e_spec_ind, KOKKOS_LAMBDA(int ispec) {
        const auto density_ind = ispec - s_spec_ind;
        const PackType::scalar inv_density = 1.0 / density[density_ind];
        for (int pack_idx = 0; pack_idx < num_vert_packs; pack_idx++) {
          dryvol_a(pack_idx) += max(0.0, q_i(ispec, pack_idx)) * inv_density;
          dryvol_c(pack_idx) += max(0.0, q_i(ispec, pack_idx)) * inv_density;
        }
      });
}

void MAMCalcsizeProcess::set_initial_sz_and_volumes_(
    const int imode, const int top_lev, const int nlevs,
    SpeciesColumnView dgncur, SpeciesColumnView v2ncur, ColumnView dryvol,
    const std::size_t num_vert_packs) const {
  Kokkos::parallel_for(
      num_vert_packs, KOKKOS_LAMBDA(int pack_idx) {
        dgncur(imode, pack_idx) = PackType::scalar{0.0};
        v2ncur(imode, pack_idx) = PackType::scalar{0.0};
        dryvol(pack_idx) = PackType::scalar{0.0};
      });
}

}  // namespace haero
