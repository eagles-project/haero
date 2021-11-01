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

  MAMCalcsizeHostCXXProcess(const AerosolProcessType type,
                            const std::string &name,
                            const ModalAerosolConfig &config)
      : DeviceAerosolProcess<MAMCalcsizeHostCXXProcess>(type, name) {}

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

 public:
  static constexpr int top_level = 0;

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
  std::vector<Real> common_factor;

 private:
  std::shared_ptr<ekat::logger::Log::logger> logger =
      ekat::logger::Log::stdout_color_mt("console");
};

}  // namespace haero

#endif
