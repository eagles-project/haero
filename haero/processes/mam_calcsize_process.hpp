#ifndef HAERO_MAM_CALCSIZE_PROCESS_HPP
#define HAERO_MAM_CALCSIZE_PROCESS_HPP

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
      : DeviceAerosolProcess<MAMCalcsizeProcess>(pp) {}

  /// MAMCalcsizeProcess objects are not assignable.
  AerosolProcess &operator=(const MAMCalcsizeProcess &) = delete;

 protected:
  //------------------------------------------------------------------------
  //                                Overrides
  //------------------------------------------------------------------------

  void init_(const ModalAerosolConfig &modal_aerosol_config) override;

  KOKKOS_FUNCTION
  void run_(const TeamType &team, Real t, Real dt,
            const Prognostics &prognostics, const Atmosphere &atmosphere,
            const Diagnostics &diagnostics,
            const Tendencies &tendencies) const override;

 private:
  /**
   * \brief Set initial defaults for the dry diameter, volume to num and dry
   * volume.
   */
  KOKKOS_FUNCTION
  void set_initial_sz_and_volumes_(const int imode, const int top_lev,
                                   const int nlevs, SpeciesColumnView dgncur,
                                   SpeciesColumnView v2ncur, ColumnView dryvol,
                                   const std::size_t num_vert_packs) const;


  /**
   * \brief Compute initial dry volume based on bulk mass mixing ratio (mmr) and
   * specie density volume = mmr/density
   */
  KOKKOS_FUNCTION
  void compute_dry_volume(const int imode, const int top_lev, const int nlevs,
                          const int s_spec_ind, const int e_spec_ind,
                          const RealView &density,
                          const SpeciesColumnView q_i,
                          const SpeciesColumnView q_c,
                          ColumnView dryvol_a,
                          ColumnView dryvol_c,
                          const std::size_t num_vert_packs) const;

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
  KOKKOS_FUNCTION
  void get_relaxed_v2n_limits(const bool do_aitacc_transfer,
                              const bool is_aitken_mode,
                              const bool is_accum_mode, Real &v2nmin,
                              Real &v2nmax, Real &v2nminrl,
                              Real &v2nmaxrl) const;

  /*
   * \brief number adjustment routine. See the implementation for more detailed
   * comments.
   */
  void adjust_num_sizes(const PackType &drv_a, const PackType &drv_c,
                        const PackType &init_num_a, const PackType &init_num_c,
                        const Real &dt, const Real &v2nmin, const Real &v2nmax,
                        const Real &v2nminrl, const Real &v2nmaxrl,
                        PackType &num_a, PackType &num_c, PackType &dqdt,
                        PackType &dqqcwdt) const;

  static inline PackType update_number_mixing_ratio_tendencies(
      const PackType &num, const PackType &num0, const PackType &dt_inverse) {
    return (num - num0) * dt_inverse;
  }

  static inline PackType min_max_bounded(const PackType &drv,
                                         const PackType &v2nmin,
                                         const PackType &v2nmax,
                                         const PackType &num) {
    return ekat::max(drv * v2nmin, ekat::min(drv * v2nmax, num));
  }

 private:
  static constexpr int top_level = 0;
  static constexpr bool do_adjust = true;
  static constexpr bool do_aitacc_transfer = true;

  int nmodes;
  int max_nspec;
  int num_populations;
  int aitken_idx;
  int accum_idx;

  IntView population_offsets;
  IntView num_mode_species;

  // NOTE: this has been linearized just like the aero species/modes held in
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
};

}  // namespace haero

#endif
