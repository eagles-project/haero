#ifndef HAERO_MAM_CALCSIZE_HOSTCXX_PROCESS_HPP
#define HAERO_MAM_CALCSIZE_HOSTCXX_PROCESS_HPP

#include <ekat/ekat.hpp>
#include <ekat/logging/ekat_logger.hpp>
#include <iomanip>

#include "haero/aerosol_process.hpp"

#ifdef HAERO_USE_CUDA
#error "MAMCalcsizeHostCXXProcess is host-only!"
#endif

namespace haero {

/// @class MAMCalcsizeHostCXXProcess
///
/// @remark This processes is an intermediate step to safely transition from the
/// fortran to the Kokkos process.
class MAMCalcsizeHostCXXProcess final
    : public DeviceAerosolProcess<MAMCalcsizeHostCXXProcess> {
 public:
  using HostColumnVec = std::vector<PackType>;
  using HostModeColumnVec = std::vector<std::vector<PackType>>;
  using HostSpeciesColumnVec = std::vector<std::vector<PackType>>;

  MAMCalcsizeHostCXXProcess();

  MAMCalcsizeHostCXXProcess(const std::string &name,
                            const ModalAerosolConfig &config)
      : DeviceAerosolProcess<MAMCalcsizeHostCXXProcess>(name) {}

  /// Destructor.
  KOKKOS_INLINE_FUNCTION
  ~MAMCalcsizeHostCXXProcess() {}

  /// Default copy constructor. For use in moving host instance to device.
  KOKKOS_INLINE_FUNCTION
  MAMCalcsizeHostCXXProcess(const MAMCalcsizeHostCXXProcess &pp)
      : DeviceAerosolProcess<MAMCalcsizeHostCXXProcess>(pp) {}

  /// MAMCalcsizeHostCXXProcess objects are not assignable.
  AerosolProcess &operator=(const MAMCalcsizeHostCXXProcess &) = delete;

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

 public:
  /**
   * \brief Set initial defaults for the dry diameter, volume to num and dry
   * volume.
   */
  void set_initial_sz_and_volumes_(const int imode, const int top_lev,
                                   const int nlevs,
                                   std::vector<std::vector<PackType>> &dgncur,
                                   std::vector<std::vector<PackType>> &v2ncur,
                                   std::vector<PackType> &dryvol,
                                   const std::size_t num_vert_packs) const;

  /**
   * \brief Compute initial dry volume based on bulk mass mixing ratio (mmr) and
   * specie density volume = mmr/density
   */
  void compute_dry_volume(const int imode, const int top_lev, const int nlevs,
                          const int s_spec_ind, const int e_spec_ind,
                          const std::vector<Real> &density,
                          const std::vector<std::vector<PackType>> &q_i,
                          const std::vector<std::vector<PackType>> &q_c,
                          std::vector<PackType> &dryvol_a,
                          std::vector<PackType> &dryvol_c,
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
  void get_relaxed_v2n_limits(const bool do_aitacc_transfer,
                              const bool is_aitken_mode,
                              const bool is_accum_mode, Real &v2nmin,
                              Real &v2nmax, Real &v2nminrl,
                              Real &v2nmaxrl) const;

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

 public:
  static constexpr int top_level = 0;
  static constexpr bool do_adjust = true;
  static constexpr bool do_aitacc_transfer = true;

  int nmodes;
  int max_nspec;
  int num_populations;
  int aitken_idx;
  int accum_idx;

  std::vector<int> population_offsets;
  std::vector<int> num_mode_species;

  // NOTE: this has been linearized just like the aero species/modes held in
  // the modal_aerosol_config.
  std::vector<int> spec_density;

  std::vector<Real> v2nmin_nmodes;
  std::vector<Real> v2nmax_nmodes;
  std::vector<Real> dgnmin_nmodes;
  std::vector<Real> dgnmax_nmodes;

  // There is a common factor calculated over and over in the core loop of this
  // process. This factor has been pulled out so the calculation only has to be
  // performed once.
  std::vector<Real> common_factor_nmodes;

 private:
  std::shared_ptr<ekat::logger::Log::logger> logger =
      ekat::logger::Log::stdout_color_mt("console");
};

}  // namespace haero

#endif
