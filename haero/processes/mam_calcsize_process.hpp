#ifndef HAERO_MAM_CALCSIZE_PROCESS_HPP
#define HAERO_MAM_CALCSIZE_PROCESS_HPP

#include <ekat/ekat.hpp>
#include <iomanip>

#include "haero/aerosol_process.hpp"

namespace haero {

/// @class MAMCalcsizeProcess
class MAMCalcsizeProcess final
    : public DeviceAerosolProcess<MAMCalcsizeProcess> {
 public:
  using RealView = DeviceType::view_1d<Real>;
  using IntView = DeviceType::view_1d<int>;

  MAMCalcsizeProcess();

  MAMCalcsizeProcess(const std::string &name, const ModalAerosolConfig &config)
      : DeviceAerosolProcess<MAMCalcsizeProcess>(name) {}

  /// Destructor.
  KOKKOS_INLINE_FUNCTION
  ~MAMCalcsizeProcess() {}

  /// Default copy constructor. For use in moving host instance to device.
  KOKKOS_INLINE_FUNCTION
  MAMCalcsizeProcess(const MAMCalcsizeProcess &pp)
      : DeviceAerosolProcess<MAMCalcsizeProcess>(pp),
        nmodes{pp.nmodes},
        max_nspec{pp.max_nspec},
        num_populations{pp.num_populations},
        aitken_idx{pp.aitken_idx},
        accum_idx{pp.accum_idx},
        population_offsets{pp.population_offsets},
        num_mode_species{pp.num_mode_species},
        spec_density{pp.spec_density},
        v2nmin_nmodes{pp.v2nmin_nmodes},
        v2nmax_nmodes{pp.v2nmax_nmodes},
        dgnmin_nmodes{pp.dgnmin_nmodes},
        dgnmax_nmodes{pp.dgnmax_nmodes},
        common_factor_nmodes{pp.common_factor_nmodes},
        density{pp.density} {}

  /// MAMCalcsizeProcess objects are not assignable.
  AerosolProcess &operator=(const MAMCalcsizeProcess &) = delete;

 protected:
  //------------------------------------------------------------------------
  //                                Overrides
  //------------------------------------------------------------------------

  void init_(const ModalAerosolConfig &modal_aerosol_config) override;

  KOKKOS_INLINE_FUNCTION
  void run_(const TeamType &team, Real t, Real dt,
            const Prognostics &prognostics, const Atmosphere &atmosphere,
            const Diagnostics &diagnostics,
            const Tendencies &tendencies) const override {
    const int nlevels = diagnostics.num_levels();

    const auto num_vert_packs = PackInfo::num_packs(nlevels);

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

    for (int k_pack = 0; k_pack < num_vert_packs; k_pack++) {
      for (int imode = 0; imode < nmodes; imode++) {
        // species starting index in the population (q_i and q_c) arrays for a
        // mode
        const auto s_spec_idx = population_offsets(imode);

        // end index of species for all modes expect the last mode
        const auto e_spec_idx = ((1 + imode) == nmodes)
                                    ? num_populations
                                    : population_offsets(imode + 1) - 1;

        const auto nspec = num_mode_species(imode);

        for (int i = 0; i < max_nspec; i++)
          density(i) =
              std::numeric_limits<decltype(density)::value_type>::max();
        for (int i = 0; i < nspec; i++)
          density(i) = spec_density(population_offsets(imode) + i);
        auto v2nmin = v2nmin_nmodes(imode);
        auto v2nmax = v2nmax_nmodes(imode);
        const auto dgnmin = dgnmin_nmodes(imode);
        const auto dgnmax = dgnmax_nmodes(imode);
        const auto common_factor = common_factor_nmodes(imode);

        Real v2nminrl, v2nmaxrl;

        // compute upper and lower limits for volume to num (v2n) ratios and
        // diameters (dgn)
        get_relaxed_v2n_limits(do_aitacc_transfer, imode == aitken_idx,
                               imode == accum_idx, v2nmin, v2nmax, v2nminrl,
                               v2nmaxrl);

        // Interstitial aerosols
        // dryvolume pack for this pack of levels
        auto dryvol_a_k = compute_dry_volume_k(s_spec_idx, e_spec_idx, density,
                                               q_i(imode, k_pack));

        // initial value of num interstitial for this pack and mode
        auto init_num_a = n_i(imode, k_pack);

        // `adjust_num_sizes` will use the initial value, but other
        // calculations require this to be nonzero.
        auto num_a = PackType(init_num_a < 0, PackType(0.0), init_num_a);

        // Same calculations as above, but for cloudborne aerosols
        auto dryvol_c_k = compute_dry_volume_k(s_spec_idx, e_spec_idx, density,
                                               q_c(imode, k_pack));

        auto init_num_c = n_c(imode, k_pack);
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
           * Adjustments that are applied over time-scale of a day (in
           *seconds)
           *   3. bring numbers to within specified bounds
           *
           * (Adjustment details are explained in the process)
           *------------------------------------------------------------------*/

          // number tendencies to be updated by adjust_num_sizes subroutine
          auto interstitial_tend = dnidt(imode, k_pack);
          auto cloudborne_tend = dncdt(imode, k_pack);

          adjust_num_sizes(dryvol_a_k, dryvol_c_k, init_num_a, init_num_c, dt,
                           v2nmin, v2nmax, v2nminrl, v2nmaxrl, num_a, num_c,
                           interstitial_tend, cloudborne_tend);
        }
      }
    }
  }

 private:
  /**
   * \brief Set initial defaults for the dry diameter, volume to num and dry
   * volume.
   */
  KOKKOS_INLINE_FUNCTION
  void set_initial_sz_and_volumes_(PackType &dgncur, PackType &v2ncur) const {
    dgncur = PackType::scalar{0.0};
    v2ncur = PackType::scalar{0.0};
  }

  /**
   * \brief Compute initial dry volume based on bulk mass mixing ratio (mmr) and
   * specie density volume = mmr/density
   */
  KOKKOS_INLINE_FUNCTION
  void compute_dry_volume(const int imode, const int top_lev, const int nlevs,
                          const int s_spec_idx, const int e_spec_idx,
                          const RealView &density, const SpeciesColumnView q_i,
                          const SpeciesColumnView q_c, ColumnView dryvol_a,
                          ColumnView dryvol_c,
                          const std::size_t num_vert_packs) const {
    using namespace ekat;
    for (int ispec = 0; ispec < e_spec_idx; ispec++) {
      const auto density_idx = ispec - s_spec_idx;
      const PackType::scalar inv_density = 1.0 / density[density_idx];
      for (int k_pack = 0; k_pack < num_vert_packs; k_pack++) {
        dryvol_a(k_pack) += max(0.0, q_i(ispec, k_pack)) * inv_density;
        dryvol_c(k_pack) += max(0.0, q_i(ispec, k_pack)) * inv_density;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  PackType compute_dry_volume_k(const int s_spec_idx, const int e_spec_idx,
                                const RealView &density,
                                const PackType &q_k) const {
    using namespace ekat;
    PackType dryvol = 0;
    for (int ispec = s_spec_idx; ispec < e_spec_idx; ispec++) {
      const auto density_idx = ispec - s_spec_idx;
      const PackType::scalar inv_density = 1.0 / density[density_idx];
      dryvol += max(0.0, q_k) * inv_density;
    }
    return dryvol;
  }

  /*
   * \brief Get relaxed limits for volume_to_num (we use relaxed limits for
   * aerosol number "adjustment" calculations via "adjust_num_sizes" subroutine.
   * Note: The relaxed limits will be artifically inflated (or deflated) for the
   * aitken and accumulation modes if "do_aitacc_transfer" flag is true to
   * effectively shut-off aerosol number "adjustment" calculations for these
   * modes because we do the explicit transfer (via "aitken_accum_exchange"
   * subroutine) from one mode to another instead of adjustments for these
   * modes).
   *
   * \note v2nmin and v2nmax are only updated for aitken and accumulation modes.
   */
  KOKKOS_INLINE_FUNCTION
  void get_relaxed_v2n_limits(const bool do_aitacc_transfer,
                              const bool is_aitken_mode,
                              const bool is_accum_mode, Real &v2nmin,
                              Real &v2nmax, Real &v2nminrl,
                              Real &v2nmaxrl) const {
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
  }

  /*
   * \brief number adjustment routine. See the implementation for more detailed
   * comments.
   */
  KOKKOS_INLINE_FUNCTION
  void adjust_num_sizes(const PackType &drv_a, const PackType &drv_c,
                        const PackType &init_num_a, const PackType &init_num_c,
                        const Real &dt, const Real &v2nmin, const Real &v2nmax,
                        const Real &v2nminrl, const Real &v2nmaxrl,
                        PackType &num_a, PackType &num_c, PackType &dqdt,
                        PackType &dqqcwdt) const {
    static constexpr Real close_to_one = 1.0 + 1.0e-15;
    static constexpr Real seconds_in_a_day = 86400.0;

    /*
     *
     * The logic behind the number adjustment is described in detail in the
     * "else" section of the following "if" condition.
     *
     * We accomplish number adjustments in 3 steps:
     *
     *   1. Ensure that number mixing ratios are either zero or positive to
     * begin with. If both of them are zero (or less), we make them zero and
     * update tendencies accordingly (logic in the first "if" block")
     *   2. In this step, we use "relaxed" bounds for bringing number mixing
     *      ratios in their bounds. This is accomplished in three sub-steps
     * [(a), (b) and (c)] described in "Step 2" below.
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
    const auto frac_adj_in_dt =
        std::max(0.0, std::min(1.0, dt * adj_tscale_inv));

    // inverse of time step
    const auto dtinv = 1.0 / (dt * close_to_one);

    /*
     * The masks below represent four if-else conditions in the original fortran
     * code. The masks represent whether a given branch should be traversed for
     * a given element of the pack, and this pack is passed to the function
     * invokations.
     */
    const auto drva_le_zero = drv_a <= 0.0;
    num_a.set(drva_le_zero, 0.0);

    const auto drvc_le_zero = drv_c <= 0.0;
    num_c.set(drvc_le_zero, 0.0);

    /* If both interstitial (drv_a) and cloud borne (drv_c) dry volumes are zero
     * (or less) adjust numbers(num_a and num_c respectively) for both of them
     * to be zero for this mode and level
     */
    const auto drv_a_c_le_zero = drva_le_zero && drvc_le_zero;
    dqdt.set(drv_a_c_le_zero, compute_tendency(num_a, init_num_a, dtinv));
    dqqcwdt.set(drv_a_c_le_zero, compute_tendency(num_c, init_num_c, dtinv));

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

    /* Note that anything in this scope that touches a pack outside this scope,
     * it must also refer to `drv_a_c_gt_zero`. eg `pk.set(drv_a_c_gt_zero &&
     * some_other_cond, val);`
     */
    const auto drv_a_c_gt_zero = !drvc_le_zero && !drva_le_zero;
    if (drv_a_c_gt_zero.any()) {
      /*
       * The number adjustment is done in 3 steps:
       *
       * Step 1: assumes that num_a and num_c are non-negative (nothing to be
       * done here)
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
       * interstitial or cloud- borne is changing. If interstitial stayed the
       * same (i.e. it is within range) but cloud-borne is predicted to reach
       * its maximum(or minimum), we modify interstitial number (num_a), so as
       * to accomodate change in the cloud-borne aerosols (and vice-versa). We
       * try to balance these by moving the num_* in the opposite direction as
       * much as possible to conserve num_a + num_c (such that num_a+num_c stays
       * close to its original value)
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
        const auto delta_num_t3 =
            (min_number_bound - total_num) * frac_adj_in_dt;

        /*
         * Now we need to decide how to distribute "delta_num" (change in
         * number) for num_a and num_c.
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
        delta_num_c_stp3.set(
            total_lt_lowerbound && (num_c_stp2 < drv_c * v2nmin), delta_num_t3);

        // if only num_a is less than lower bound, assign total change to num_a
        delta_num_a_stp3.set(
            total_lt_lowerbound && (num_a_stp2 < drv_a * v2nmin), delta_num_t3);
      }

      const auto total_gt_upperbound = total_num > max_number_bound;
      {
        // change in total_num in one time step
        const auto delta_num_t3 =
            (max_number_bound - total_num) * frac_adj_in_dt;

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
        delta_num_c_stp3.set(
            total_gt_upperbound && (num_c_stp2 > drv_c * v2nmax), delta_num_t3);

        // if only num_a is more than the upper bound, assign total change to
        // num_a
        delta_num_a_stp3.set(
            total_gt_upperbound && (num_a_stp2 > drv_a * v2nmax), delta_num_t3);
      }

      // Update num_a/c
      num_a.set(drv_a_c_gt_zero, num_a_stp2 + delta_num_a_stp3);
      num_c.set(drv_a_c_gt_zero, num_c_stp2 + delta_num_c_stp3);
    }

    // Update tendencies
    dqdt = compute_tendency(num_a, init_num_a, dtinv);
    dqqcwdt = compute_tendency(num_c, init_num_c, dtinv);
  }

  KOKKOS_INLINE_FUNCTION
  static PackType compute_tendency(const PackType &num, const PackType &num0,
                                   const PackType &dt_inverse) {
    return (num - num0) * dt_inverse;
  }

  KOKKOS_INLINE_FUNCTION
  static PackType min_max_bounded(const PackType &drv, const PackType &v2nmin,
                                  const PackType &v2nmax, const PackType &num) {
    return ekat::max(drv * v2nmin, ekat::min(drv * v2nmax, num));
  }

 private:
  static constexpr int top_level = 0;
  static constexpr bool do_adjust = true;
  static constexpr bool do_aitacc_transfer = true;

  std::size_t nmodes;
  std::size_t max_nspec;
  std::size_t num_populations;
  std::size_t aitken_idx;
  std::size_t accum_idx;

  IntView population_offsets;
  IntView num_mode_species;

  // NOTE: this has been flattened just like the aero species/modes held in
  // the modal_aerosol_config.
  IntView spec_density;

  RealView v2nmin_nmodes;
  RealView v2nmax_nmodes;
  RealView dgnmin_nmodes;
  RealView dgnmax_nmodes;

  // There is a common factor calculated over and over in the core loop of this
  // process. This factor has been pulled out so the calculation only has to be
  // performed once.
  RealView common_factor_nmodes;

  // The following are only used in run
  mutable RealView density;
};

}  // namespace haero

#endif
