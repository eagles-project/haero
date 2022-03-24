#ifndef HAERO_MAM_CALCSIZE_PROCESS_HPP
#define HAERO_MAM_CALCSIZE_PROCESS_HPP

#include <cmath>
#include <ekat/ekat.hpp>
#include <ekat/ekat_pack.hpp>
#include <iomanip>

#include "haero/aerosol_process.hpp"

namespace haero {

/// @class MAMCalcsizeProcess
class MAMCalcsizeProcess final
    : public DeviceAerosolProcess<MAMCalcsizeProcess> {
 public:
  using RealView = DeviceType::view_1d<Real>;
  using IntView = DeviceType::view_1d<int>;
  using IntView2D = DeviceType::view_2d<int>;

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
        v2nnom_nmodes{pp.v2nnom_nmodes},
        dgnmin_nmodes{pp.dgnmin_nmodes},
        dgnmax_nmodes{pp.dgnmax_nmodes},
        dgnnom_nmodes{pp.dgnnom_nmodes},
        common_factor_nmodes{pp.common_factor_nmodes},
        density{pp.density},
        srcmode_csizxf{pp.srcmode_csizxf},
        destmode_csizxf{pp.destmode_csizxf},
        ait_mode_inter{pp.ait_mode_inter},
        acc_mode_inter{pp.acc_mode_inter},
        ait_mode_cldbrn{pp.ait_mode_cldbrn},
        acc_mode_cldbrn{pp.acc_mode_cldbrn},
        do_aitacc_transfer_allowed{pp.do_aitacc_transfer_allowed},
        ait_spec_in_acc_inter{pp.ait_spec_in_acc_inter},
        acc_spec_in_ait_inter{pp.acc_spec_in_ait_inter},
        ait_spec_in_acc_cldbrn{pp.ait_spec_in_acc_cldbrn},
        acc_spec_in_ait_cldbrn{pp.acc_spec_in_ait_cldbrn},
        dgncur_a{pp.dgncur_a},
        dgncur_c{pp.dgncur_c},
        v2ncur_a{pp.v2ncur_a},
        v2ncur_c{pp.v2ncur_c},
        dryvol_a{pp.dryvol_a},
        dryvol_c{pp.dryvol_c},
        drv_a_aitsv{pp.drv_a_aitsv},
        drv_a_accsv{pp.drv_a_accsv},
        drv_c_aitsv{pp.drv_c_aitsv},
        drv_c_accsv{pp.drv_c_accsv},
        num_a_aitsv{pp.num_a_aitsv},
        num_a_accsv{pp.num_a_accsv},
        num_c_aitsv{pp.num_c_aitsv},
        num_c_accsv{pp.num_c_accsv},
        drv_a_sv{pp.drv_a_sv},
        drv_c_sv{pp.drv_c_sv},
        num_a_sv{pp.num_a_sv},
        num_c_sv{pp.num_c_sv},
        nspec_common{pp.nspec_common} {}

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

    // See declaration of num_levels_upper_bound for its documentation
    EKAT_KERNEL_ASSERT(nlevels <= num_levels_upper_bound);

    const auto num_vert_packs = PackInfo::num_packs(nlevels);

    // interstitial mass and number mixing ratios
    const auto q_i = prognostics.interstitial_aerosols;
    const auto n_i = prognostics.interstitial_num_mix_ratios;

    // cloud-borne mass and number mixing ratios
    const auto q_c = prognostics.cloud_aerosols;
    const auto n_c = prognostics.cloud_num_mix_ratios;

    // tendencies for interstitial number mixing ratios
    const ModeColumnView dnidt = tendencies.interstitial_num_mix_ratios;

    // tendencies for cloud-borne number mixing ratios
    const ModeColumnView dncdt = tendencies.cloud_num_mix_ratios;

    // tendencies for interstitial mass mixing ratios
    const SpeciesColumnView didt = tendencies.interstitial_aerosols;

    // tendencies for cloudborne mass mixing ratios
    const SpeciesColumnView dcdt = tendencies.cloud_aerosols;

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
          density(i) = i < nspec ? spec_density(population_offsets(imode) + i)
                                 : max_real;
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
        auto num_a_k = PackType(init_num_a < 0, PackType(0.0), init_num_a);

        // Same calculations as above, but for cloudborne aerosols
        auto dryvol_c_k = compute_dry_volume_k(s_spec_idx, e_spec_idx, density,
                                               q_c(imode, k_pack));

        auto init_num_c = n_c(imode, k_pack);
        auto num_c_k = PackType(init_num_c < 0, PackType(0.0), init_num_c);

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
                           v2nmin, v2nmax, v2nminrl, v2nmaxrl, num_a_k, num_c_k,
                           interstitial_tend, cloudborne_tend);
        }

        // update diameters and volume to num ratios for interstitial aerosols
        update_diameter_and_vol2num(k_pack, imode, dryvol_a_k, num_a_k, v2nmin,
                                    v2nmax, dgnmin, dgnmax, common_factor,
                                    dgncur_a, v2ncur_a);

        // update diameters and volume to num ratios for cloudborne aerosols
        update_diameter_and_vol2num(k_pack, imode, dryvol_c_k, num_c_k, v2nmin,
                                    v2nmax, dgnmin, dgnmax, common_factor,
                                    dgncur_c, v2ncur_c);

        if (do_aitacc_transfer) {
          if (imode == aitken_idx) {
            drv_a_aitsv(k_pack) = dryvol_a_k;
            num_a_aitsv(k_pack) = num_a_k;
            drv_c_aitsv(k_pack) = dryvol_c_k;
            num_c_aitsv(k_pack) = num_c_k;
          } else if (imode == accum_idx) {
            drv_a_accsv(k_pack) = dryvol_a_k;
            num_a_accsv(k_pack) = num_a_k;
            drv_c_accsv(k_pack) = dryvol_c_k;
            num_c_accsv(k_pack) = num_c_k;
          }
        }
        drv_a_sv(k_pack, imode) = dryvol_a_k;
        num_a_sv(k_pack, imode) = num_a_k;
        drv_c_sv(k_pack, imode) = dryvol_c_k;
        num_c_sv(k_pack, imode) = num_c_k;
      }
    }

    if (do_aitacc_transfer) {
      aitken_accum_exchange(nlevels, dt, q_i, q_c, n_i, n_c, didt, dcdt, dnidt,
                            dncdt);
    }
  }

 private:
  /**
   * \brief Set initial defaults for the dry diameter, volume to num and dry
   * volume.
   */
  KOKKOS_INLINE_FUNCTION
  void set_initial_sz_and_volumes_(std::size_t imode, PackType &dgncur,
                                   PackType &v2ncur) const {
    dgncur = dgnnom_nmodes(imode);
    v2ncur = v2nnom_nmodes(imode);
  }
  KOKKOS_INLINE_FUNCTION
  void aitken_accum_exchange(
      const int nlevs, const Real dt, const SpeciesColumnView q_i,
      const SpeciesColumnView q_c, const SpeciesColumnView n_i,
      const SpeciesColumnView n_c, const SpeciesColumnView didt,
      const SpeciesColumnView dcdt, const ModeColumnView dnidt,
      const ModeColumnView dncdt) const {
    // foo
  }

  /*----------------------------------------------------------------------------
   * Compute particle diameter and volume to number ratios using dry bulk volume
   * (drv)
   *--------------------------------------------------------------------------*/
  KOKKOS_INLINE_FUNCTION
  void update_diameter_and_vol2num(std::size_t klev, std::size_t imode,
                                   PackType drv, PackType num, Real v2nmin,
                                   Real v2nmax, Real dgnmin, Real dgnmax,
                                   Real cmn_factor, SpeciesColumnView dgncur,
                                   SpeciesColumnView v2ncur) const {
    const auto drv_gt_0 = drv > 0.0;
    if (!drv_gt_0.any()) return;

    const auto drv_mul_v2nmin = drv * v2nmin;
    const auto drv_mul_v2nmax = drv * v2nmax;

    auto &dgncur_k_i = dgncur(klev, imode);
    auto &v2ncur_k_i = v2ncur(klev, imode);

    dgncur_k_i.set(num <= drv_mul_v2nmin, dgnmin);
    dgncur_k_i.set(num >= drv_mul_v2nmax, dgnmax);
    dgncur_k_i.set(num > drv_mul_v2nmin and num < drv_mul_v2nmax,
                   pow((drv / (cmn_factor * num)), (1.0 / 3)));

    v2ncur_k_i.set(num <= drv_mul_v2nmin, v2nmin);
    v2ncur_k_i.set(num >= drv_mul_v2nmax, v2nmax);
    v2ncur_k_i.set(num > drv_mul_v2nmin and num < drv_mul_v2nmax, num / drv);
  }

  /*
   * \brief Find mapping between the species of two different modes
   * \note Global variables declared at modeul level will be updated to store
   * the mapping generated by the following call
   */
  void find_species_mapping(ModalAerosolConfig const &config);

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
        dryvol_a(k_pack) += ekat::max(0.0, q_i(ispec, k_pack)) * inv_density;
        dryvol_c(k_pack) += ekat::max(0.0, q_i(ispec, k_pack)) * inv_density;
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
      dryvol += ekat::max(0.0, q_k) * inv_density;
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
    const auto adj_tscale = fmax(seconds_in_a_day, dt);

    // inverse of the adjustment time scale
    const auto adj_tscale_inv = 1.0 / (adj_tscale * close_to_one);

    // fraction of adj_tscale covered in the current time step "dt"
    const auto frac_adj_in_dt = fmax(0.0, fmin(1.0, dt * adj_tscale_inv));

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
  /*
   * num_levels_upper_bound is a naive way of allocating enough memory in work
   * arrays to operate on data that depends on the number of levels currently
   * being used.
   *
   * For example, dgncur_a has shape (num_levels, num_modes). The number of
   * levels is not known until the run member function is called. All allocation
   * must be performed *before* the run method is called to keep the codebase
   * portable.
   *
   * This is a depiction of the circular dependency that necessitates this
   * assumption about the number of levels:
   *
   *                    ┌────────┐      ┌──────────┐
   *                    │dgncur_a├─────►│num levels│
   *                    └────────┘      └─────┬────┘
   *                         ▲                │
   *                         │  ┌──────────┐  │
   *                         └──┤run method│◄─┘
   *                            └──────────┘
   *
   * To remove this assumption, we will have to be able to pass the
   * number of levels to the init member function, calculate an upper bound for
   * the number of levels in the init member function, or perform device
   * allocation in a portable manner.
   */
  static constexpr std::size_t num_levels_upper_bound = 128;

  static constexpr int top_level = 0;
  static constexpr bool do_adjust = true;
  static constexpr bool do_aitacc_transfer = true;
  static constexpr auto max_real = std::numeric_limits<Real>::max();

  // Max number of radiation diagnostics FIXME: it should be set somewhere in
  // the "init" codes for the model
  static constexpr std::size_t n_diag = 0;

  // "list_idx=0" is reservered for the prognostic call
  // FIXME: We are currently supporting only list_idx=0, generalize and test the
  // code with other list_idx values
  static constexpr std::size_t list_idx = 0;

  // Maximum number of aitken-accumulation pairs we can have (one for each
  // diagnostic list). NOTE: "0" is reserved for the prognostic call, so
  // maxpair_csizxf represents just the diagnostic pairs.
  static constexpr std::size_t maxpair_csizxf = n_diag;

  // Total number of pairs of aitken-accumulation modes
  // ---------------------------------------------------------------------------------
  // Note: For diagnostic calls, users can ask for any diagnostics like
  // rad_diag_1
  //       and rad_diag_3 (e.g. skipping rad_diag_2). Therefore the arrays
  //       should have a length of N_DIAG (unless we define another array which
  //       maps info such as to include info about the missing diagnostics)
  // ---------------------------------------------------------------------------------

  // total number of possible diagnostic calls
  static constexpr std::size_t npair_csizxf = n_diag;

  // Magic number taken from fortran process
  static constexpr std::size_t index_mapping_extent_1 = 10;

  // "srcmode_csizxf" stores source mode number from  which species will be
  // moved to a mode stored in "destmode_csizxf". [E.g. if srcmode_csizxf(3)=2
  // and destmode_csizxf(3)=1, for rad_diag_3 (notice that these arrays are
  // indexed "3" as rad_diag is rad_diag_3), species will be moved from 2nd mode
  // to the 1st mode.]
  IntView srcmode_csizxf;
  IntView destmode_csizxf;
  IntView ait_mode_inter;
  IntView acc_mode_inter;
  IntView ait_mode_cldbrn;
  IntView acc_mode_cldbrn;
  IntView do_aitacc_transfer_allowed;

  IntView2D ait_spec_in_acc_inter;
  IntView2D acc_spec_in_ait_inter;
  IntView2D ait_spec_in_acc_cldbrn;
  IntView2D acc_spec_in_ait_cldbrn;

  // interstitial particle diameter[m](FIXME: This should be diagnostic
  // variable)
  SpeciesColumnView dgncur_a;
  // cldborne particle diameter [m](FIXME: This should be diagnostic variable)
  SpeciesColumnView dgncur_c;
  // interstitial vol2num ratio [FIXME: units????]
  SpeciesColumnView v2ncur_a;
  // cldborne vol2num ratio [FIXME:units???]
  SpeciesColumnView v2ncur_c;

  // dry volume of a particle[m3/kg(of air)]
  RealView dryvol_a;
  RealView dryvol_c;

  // Work variables for aitken<-->accumulation transfer sub process
  ColumnView drv_a_aitsv;

  // saves aitken and accumulation interstitial modes dryvolume
  ColumnView drv_a_accsv;
  ColumnView drv_c_aitsv;

  // saves aitken and accumulation cloudborne modes dryvolume
  ColumnView drv_c_accsv;
  ColumnView num_a_aitsv;

  // saves aitken and accumulation interstitial modes num concentrations
  ColumnView num_a_accsv;
  ColumnView num_c_aitsv;

  // saves aitken and accumulation cloudborne modes num concentrations
  ColumnView num_c_accsv;

  SpeciesColumnView drv_a_sv;

  // saves dryvolume for each mode and level
  SpeciesColumnView drv_c_sv;
  SpeciesColumnView num_a_sv;

  // saves num conc. for each mode and level
  SpeciesColumnView num_c_sv;

  // number of species found common between aitken and accumulation modes
  IntView nspec_common;

  std::size_t nmodes;
  std::size_t max_nspec;
  std::size_t num_populations;
  std::size_t aitken_idx;
  std::size_t accum_idx;

  IntView population_offsets;
  IntView num_mode_species;

  // NOTE: this has been flattened just like the aero species/modes held in
  // the modal_aerosol_config.
  RealView spec_density;

  RealView v2nmin_nmodes;
  RealView v2nmax_nmodes;
  RealView v2nnom_nmodes;
  RealView dgnmin_nmodes;
  RealView dgnmax_nmodes;
  RealView dgnnom_nmodes;

  // There is a common factor calculated over and over in the core loop of this
  // process. This factor has been pulled out so the calculation only has to be
  // performed once.
  RealView common_factor_nmodes;

  // The following are only used in run
  mutable RealView density;
};

}  // namespace haero

#endif
