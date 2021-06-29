#ifndef HAERO_MAM_CALCSIZE_PROCESS_HPP
#define HAERO_MAM_CALCSIZE_PROCESS_HPP

#include <iomanip>

#include "haero/aerosol_process.hpp"

namespace haero {

/// @class MAMCalcsizeProcess
class MAMCalcsizeProcess : public AerosolProcess {

 public:
  MAMCalcsizeProcess();

  MAMCalcsizeProcess(const AerosolProcessType type, const std::string &name,
                       const ModalAerosolConfig &config)
    : AerosolProcess(type, name){ }

  /// Destructor.
  KOKKOS_INLINE_FUNCTION
  virtual ~MAMCalcsizeProcess() {}

  /// Default copy constructor. For use in moving host instance to device.
  KOKKOS_INLINE_FUNCTION
  MAMCalcsizeProcess(const MAMCalcsizeProcess &pp)
    : AerosolProcess(pp) {}

  /// MAMCalcsizeProcess objects are not assignable.
  AerosolProcess &operator=(const MAMCalcsizeProcess &) = delete;

  //------------------------------------------------------------------------
  //                                Accessors
  //------------------------------------------------------------------------

  KOKKOS_FUNCTION
  virtual void run(const ModalAerosolConfig &modal_aerosol_config, Real t,
                   Real dt, const Prognostics &prognostics,
                   const Atmosphere &atmosphere, const Diagnostics &diagnostics,
                   Tendencies &tendencies) const override{};
};

}  // namespace haero

#endif
