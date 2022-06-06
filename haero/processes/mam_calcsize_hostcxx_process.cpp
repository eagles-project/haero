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
          "MAMCalcsizeHostCXXProcess") {
  logger->set_pattern("MAMCalcsizeHostCXXProcess[%^%l%$]: %v");
}
//------------------------------------------------------------------------
//                                Accessors
//------------------------------------------------------------------------

void MAMCalcsizeHostCXXProcess::init_(
    const ModalAerosolConfig &modal_aerosol_config) {
  logger->debug("entering init_");

  nmodes = modal_aerosol_config.num_aerosol_modes();
  logger->debug("nmodes=modal_aerosol_config.num_modes={}",
                modal_aerosol_config.num_aerosol_modes());

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
  common_factor_nmodes = zero_vec(nmodes);

  for (int i = 0; i < nmodes; i++) {
    logger->debug("setting volume/num ratios and diameters mode={}", i);
    const auto &mode = modal_aerosol_config.aerosol_modes[i];
    using T = decltype(v2nmin_nmodes)::value_type;
    v2nmin_nmodes[i] = mode.min_vol_to_num_ratio<T>();
    v2nmax_nmodes[i] = mode.max_vol_to_num_ratio<T>();
    dgnmin_nmodes[i] = mode.min_diameter;
    dgnmax_nmodes[i] = mode.max_diameter;
    common_factor_nmodes[i] =
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

  const int nlevels = prognostics.num_levels();

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

    compute_dry_volume(imode, top_level, nlevels, start_spec_idx, end_spec_idx,
                       density, q_i, q_c, dryvol_a, dryvol_c, num_vert_packs);

    // compute upper and lower limits for volume to num (v2n) ratios and
    // diameters (dgn)
    auto v2nmin = v2nmin_nmodes[imode];
    auto v2nmax = v2nmax_nmodes[imode];

    Real v2nminrl, v2nmaxrl;

    get_relaxed_v2n_limits(do_aitacc_transfer, imode == aitken_idx,
                           imode == accum_idx, v2nmin, v2nmax, v2nminrl,
                           v2nmaxrl);

    for (int pack_idx = 0; pack_idx < num_vert_packs; pack_idx++) {
      // Interstitial aerosols
      // dryvolume pack for this pack of levels
      auto dryvol_a_lvl = dryvol_a[pack_idx];

      // initial value of num interstitial for this pack and mode
      auto init_num_a = n_i[imode][pack_idx];

      // `adjust_num_sizes` will use the initial value, but other calculations
      // require this to be nonzero.
      auto num_a = PackType(init_num_a < 0, PackType(0.0), init_num_a);

      // Same calculations as above, but for cloudborne aerosols
      auto dryvol_c_lvl = dryvol_c[pack_idx];
      auto init_num_c = n_c[imode][pack_idx];
      auto num_c = PackType(init_num_c < 0, PackType(0.0), init_num_c);

      if (do_adjust) {
        /*------------------------------------------------------------------
         *  Do number adjustment for interstitial and activated particles
         *------------------------------------------------------------------
         * Adjustments that are applied over time-scale deltat
         * (model time step in seconds):
         *
         *   1. make numbers non-negative or
         *   2. make numbers zero when volume is zero
         *
         *
         * Adjustments that are applied over time-scale of a day (in seconds)
         *   3. bring numbers to within specified bounds
         *
         * (Adjustment details are explained in the process)
         *------------------------------------------------------------------*/

        // number tendencies to be updated by adjust_num_sizes subroutine
        auto interstitial_tend = dnidt[imode][pack_idx];
        auto cloudborne_tend = dncdt[imode][pack_idx];

        adjust_num_sizes(dryvol_a_lvl, dryvol_c_lvl, init_num_a, init_num_c, dt,
                         v2nmin, v2nmax, v2nminrl, v2nmaxrl, num_a, num_c,
                         interstitial_tend, cloudborne_tend);

        logger->debug(
            "performing number "
            "adjustment:interstitial_tend={},cloudborne_tend={}",
            interstitial_tend, cloudborne_tend);
      }
    }
  }
  logger->debug("leaving run_");
}

void MAMCalcsizeHostCXXProcess::adjust_num_sizes(
    const PackType &drv_a, const PackType &drv_c, const PackType &init_num_a,
    const PackType &init_num_c, const Real &dt, const Real &v2nmin,
    const Real &v2nmax, const Real &v2nminrl, const Real &v2nmaxrl,
    PackType &num_a, PackType &num_c, PackType &dqdt, PackType &dqqcwdt) const {
  logger->debug("enter adjust_num_sizes");

  static constexpr Real close_to_one = 1.0 + 1.0e-15;
  static constexpr Real seconds_in_a_day = 86400.0;

  /*
   *
   * The logic behind the number adjustment is described in detail in the "else"
   * section of the following "if" condition.
   *
   * We accomplish number adjustments in 3 steps:
   *
   *   1. Ensure that number mixing ratios are either zero or positive to begin
   *      with. If both of them are zero (or less), we make them zero and update
   *      tendencies accordingly (logic in the first "if" block")
   *   2. In this step, we use "relaxed" bounds for bringing number mixing
   *      ratios in their bounds. This is accomplished in three sub-steps [(a),
   *      (b) and (c)] described in "Step 2" below.
   *   3. In this step, we use the actual bounds for bringing number mixing
   *      ratios in their bounds. This is also accomplished in three sub-steps
   *      [(a), (b) and (c)] described in "Step 3" below.
   *
   * If the number mixing ratio in a mode is out of mode's min/max range, we
   * re-balance interstitial and cloud borne aerosols such that the number
   * mixing ratio comes within the range. Time step for such an operation is
   * assumed to be one day (in seconds). That is, it is assumed that number
   * mixing ratios will be within range in a day. "adj_tscale" represents that
   * time scale
   *
   */

  // time scale for number adjustment
  const auto adj_tscale = std::max(seconds_in_a_day, dt);

  // inverse of the adjustment time scale
  const auto adj_tscale_inv = 1.0 / (adj_tscale * close_to_one);

  // fraction of adj_tscale covered in the current time step "dt"
  const auto frac_adj_in_dt = std::max(0.0, std::min(1.0, dt * adj_tscale_inv));

  // inverse of time step
  const auto dtinv = 1.0 / (dt * close_to_one);

  logger->debug("adj_tscale={},adj_tscale_inv={},frac_adj_in_dt={},dtinv={}",
                adj_tscale, adj_tscale_inv, frac_adj_in_dt, dtinv);

  /*
   * The masks below represent four if-else conditions in the original fortran
   * code. The masks represent whether a given branch should be traversed for a
   * given element of the pack, and this pack is passed to the function
   * invokations.
   */
  const auto drva_le_zero = drv_a <= 0.0;
  num_a.set(drva_le_zero, 0.0);

  const auto drvc_le_zero = drv_c <= 0.0;
  num_c.set(drvc_le_zero, 0.0);

  /* If both interstitial (drv_a) and cloud borne (drv_c) dry volumes are zero
   * (or less) adjust numbers(num_a and num_c respectively) for both of them to
   * be zero for this mode and level
   */
  const auto drv_a_c_le_zero = drva_le_zero && drvc_le_zero;
  dqdt.set(drv_a_c_le_zero,
           update_number_mixing_ratio_tendencies(num_a, init_num_a, dtinv));
  dqqcwdt.set(drv_a_c_le_zero,
              update_number_mixing_ratio_tendencies(num_c, init_num_c, dtinv));

  /* if cloud borne dry volume (drv_c) is zero(or less), the interstitial
   * number/volume == total/combined apply step 1 and 3, but skip the relaxed
   * adjustment (step 2, see below)
   */
  const auto only_drvc_le_zero = !drva_le_zero && drvc_le_zero;
  {
    const auto numbnd = min_max_bounded(drv_a, v2nmin, v2nmax, num_a);
    num_a.set(only_drvc_le_zero, num_a + (numbnd - num_a) * frac_adj_in_dt);
  }

  /* interstitial volume is zero, treat similar to above */
  const auto only_drva_le_zero = !drvc_le_zero && drva_le_zero;
  {
    const auto numbnd = min_max_bounded(drv_c, v2nmin, v2nmax, num_c);
    num_c.set(only_drva_le_zero, num_c + (numbnd - num_c) * frac_adj_in_dt);
  }

  /* Note that anything in this scope that touches a pack outside this scope, it
   * must also refer to `drv_a_c_gt_zero`. eg
   * `pk.set(drv_a_c_gt_zero && some_other_cond, val);`
   */
  const auto drv_a_c_gt_zero = !drvc_le_zero && !drva_le_zero;
  if (drv_a_c_gt_zero.any()) {
    /*
     * The number adjustment is done in 3 steps:
     *
     * Step 1: assumes that num_a and num_c are non-negative (nothing to be done
     * here)
     */
    const auto num_a_stp1 = num_a;
    const auto num_c_stp1 = num_c;

    /*
     * Step 2 [Apply relaxed bounds] has 3 parts (a), (b) and (c)
     *
     * Step 2: (a) Apply relaxed bounds to bound num_a and num_c within
     * "relaxed" bounds.
     */
    auto numbnd = min_max_bounded(drv_a, v2nminrl, v2nmaxrl, num_a_stp1);

    /*
     * 2(b) Ideally, num_* should be in range. If they are not, we assume
     * that they will reach their maximum (or minimum)for this mode
     * within a day (time scale). We then compute how much num_* will
     * change in a time step by multiplying the difference between num_*
     * and its maximum(or minimum) with "frac_adj_in_dt".
     */
    const auto delta_num_a_stp2 = (numbnd - num_a_stp1) * frac_adj_in_dt;

    // change in num_a in one time step
    auto num_a_stp2 = num_a_stp1 + delta_num_a_stp2;

    // bounded to relaxed min and max
    numbnd = min_max_bounded(drv_c, v2nminrl, v2nmaxrl, num_c_stp1);
    const auto delta_num_c_stp2 = (numbnd - num_c_stp1) * frac_adj_in_dt;

    // change in num_a in one time step
    auto num_c_stp2 = num_c_stp1 + delta_num_c_stp2;

    /*
     * 2(c) We now also need to balance num_* incase only one among the
     * interstitial or cloud- borne is changing. If interstitial stayed the same
     * (i.e. it is within range) but cloud-borne is predicted to reach its
     * maximum(or minimum), we modify interstitial number (num_a), so as to
     * accomodate change in the cloud-borne aerosols (and vice-versa). We try to
     * balance these by moving the num_* in the opposite direction as much as
     * possible to conserve num_a + num_c (such that num_a+num_c stays close to
     * its original value)
     */
    const auto delta_num_a_stp2_eq0 = delta_num_a_stp2 == 0.0;
    const auto delta_num_c_stp2_eq0 = delta_num_c_stp2 == 0.0;

    num_a_stp2.set(delta_num_a_stp2_eq0 && !delta_num_c_stp2_eq0,
                   min_max_bounded(drv_a, v2nminrl, v2nmaxrl,
                                   num_a_stp1 - delta_num_c_stp2));

    num_c_stp2.set(delta_num_c_stp2_eq0 && !delta_num_a_stp2_eq0,
                   min_max_bounded(drv_c, v2nminrl, v2nmaxrl,
                                   num_c_stp1 - delta_num_a_stp2));

    /*
     * Step3[apply stricter bounds] has 3 parts (a), (b) and (c)
     * Step 3:(a) compute combined total of num_a and num_c
     */
    const auto total_drv = drv_a + drv_c;
    const auto total_num = num_a_stp2 + num_c_stp2;

    /*
     * 3(b) We now compute amount of num_* to change if total_num
     *     is out of range. If total_num is within range, we don't do anything
     * (i.e. delta_numa3 and delta_num_c_stp3 remain zero)
     */
    auto delta_num_a_stp3 = PackType(0.0);
    auto delta_num_c_stp3 = PackType(0.0);

    /*
     * "total_drv*v2nmin" represents minimum number for this mode, and
     * "total_drv*v2nmxn" represents maximum number for this mode
     */
    const auto min_number_bound = total_drv * v2nmin;
    const auto max_number_bound = total_drv * v2nmax;

    const auto total_lt_lowerbound = total_num < min_number_bound;
    {
      // change in total_num in one time step
      const auto delta_num_t3 = (min_number_bound - total_num) * frac_adj_in_dt;

      /*
       * Now we need to decide how to distribute "delta_num" (change in number)
       * for num_a and num_c.
       *
       * if both num_a and num_c are less than the lower bound distribute
       * "delta_num" using weighted ratios
       */
      const auto do_dist_delta_num =
          (num_a_stp2 < drv_a * v2nmin) && (num_c_stp2 < drv_c * v2nmin);

      delta_num_a_stp3.set(total_lt_lowerbound && do_dist_delta_num,
                           delta_num_t3 * (num_a_stp2 / total_num));

      delta_num_c_stp3.set(total_lt_lowerbound && do_dist_delta_num,
                           delta_num_t3 * (num_c_stp2 / total_num));

      // if only num_c is less than lower bound, assign total change to num_c
      delta_num_c_stp3.set(total_lt_lowerbound && (num_c_stp2 < drv_c * v2nmin),
                           delta_num_t3);

      // if only num_a is less than lower bound, assign total change to num_a
      delta_num_a_stp3.set(total_lt_lowerbound && (num_a_stp2 < drv_a * v2nmin),
                           delta_num_t3);
    }

    const auto total_gt_upperbound = total_num > max_number_bound;
    {
      // change in total_num in one time step
      const auto delta_num_t3 = (max_number_bound - total_num) * frac_adj_in_dt;

      // decide how to distribute "delta_num"(change in number) for num_a and
      // num_c
      const auto do_dist_delta_num =
          (num_a_stp2 > drv_a * v2nmax) && (num_c_stp2 > drv_c * v2nmax);

      /*
       * if both num_a and num_c are more than the upper bound distribute
       * "delta_num" using weighted ratios
       */
      delta_num_a_stp3.set(total_gt_upperbound && do_dist_delta_num,
                           delta_num_t3 * (num_a_stp2 / total_num));
      delta_num_c_stp3.set(total_gt_upperbound && do_dist_delta_num,
                           delta_num_t3 * (num_c_stp2 / total_num));

      // if only num_c is more than the upper bound, assign total change to
      // num_c
      delta_num_c_stp3.set(total_gt_upperbound && (num_c_stp2 > drv_c * v2nmax),
                           delta_num_t3);

      // if only num_a is more than the upper bound, assign total change to
      // num_a
      delta_num_a_stp3.set(total_gt_upperbound && (num_a_stp2 > drv_a * v2nmax),
                           delta_num_t3);
    }

    // Update num_a/c
    num_a.set(drv_a_c_gt_zero, num_a_stp2 + delta_num_a_stp3);
    num_c.set(drv_a_c_gt_zero, num_c_stp2 + delta_num_c_stp3);
  }

  // Update tendencies
  dqdt = update_number_mixing_ratio_tendencies(num_a, init_num_a, dtinv);
  dqqcwdt = update_number_mixing_ratio_tendencies(num_c, init_num_c, dtinv);

  logger->debug("exit adjust_num_sizes");
}

void MAMCalcsizeHostCXXProcess::get_relaxed_v2n_limits(
    const bool do_aitacc_transfer, const bool is_aitken_mode,
    const bool is_accum_mode, Real &v2nmin, Real &v2nmax, Real &v2nminrl,
    Real &v2nmaxrl) const {
  logger->debug(
      "get_relaxed_v2n_limits:do_aitacc_transfer={},is_aitken_mode={},is_accum_"
      "mode={}",
      do_aitacc_transfer, is_aitken_mode, is_accum_mode);

  /*
   * Relaxation factor is currently assumed to be a factor of 3 in diameter
   * which makes it 3**3=27 for volume.  i.e. dgnumlo_relaxed = dgnumlo/3 and
   * dgnumhi_relaxed = dgnumhi*3; therefore we use 3**3=27 as a relaxation
   * factor for volume.
   *
   * \see get_relaxed_v2n_limits
   */
  static constexpr Real relax_factor = 27.0;

  // factor to artifically inflate or deflate v2nmin and v2nmax
  static constexpr Real szadj_block_fac = 1.0e6;

  if (do_aitacc_transfer) {
    if (is_aitken_mode) v2nmin /= szadj_block_fac;

    if (is_accum_mode) v2nmax *= szadj_block_fac;
  }

  v2nminrl = v2nmin / relax_factor;
  v2nmaxrl = v2nmax * relax_factor;

  logger->debug("exit get_relaxed_v2n_limits:v2nminrl={},v2nmaxrl={}", v2nminrl,
                v2nmaxrl);
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
