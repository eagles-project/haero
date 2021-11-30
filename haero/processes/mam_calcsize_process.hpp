#ifndef HAERO_MAM_CALCSIZE_PROCESS_HPP
#define HAERO_MAM_CALCSIZE_PROCESS_HPP

#include <iomanip>

#include "haero/aerosol_process.hpp"

namespace haero {

/// @class MAMCalcsizeProcess
class MAMCalcsizeProcess final
    : public DeviceAerosolProcess<MAMCalcsizeProcess> {
 public:
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
            const Tendencies &tendencies) const override{};

private:
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
};

}  // namespace haero

#endif
