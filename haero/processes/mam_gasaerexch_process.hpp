#ifndef HAERO_MAM_GASAEROSOLEXCHANGE_PROCESS_HPP
#define HAERO_MAM_GASAEROSOLEXCHANGE_PROCESS_HPP

#include <iomanip>

#include "haero/aerosol_process.hpp"

namespace haero {

class MAMGasAerosolExchangeProcess : public DeviceAerosolProcess<MAMGasAerosolExchangeProcess> {

 public:
  MAMGasAerosolExchangeProcess();

  MAMGasAerosolExchangeProcess(const AerosolProcessType type, const std::string &name,
                       const ModalAerosolConfig &config,
                       const HostDiagnostics &diagnostics)
      : DeviceAerosolProcess(type, name) { }

  /// Destructor.
  KOKKOS_INLINE_FUNCTION
  virtual ~MAMGasAerosolExchangeProcess() {}

  /// Default copy constructor. For use in moving host instance to device.
  KOKKOS_INLINE_FUNCTION
  MAMGasAerosolExchangeProcess(const MAMGasAerosolExchangeProcess &pp)
      : DeviceAerosolProcess(pp) {}

  /// MAMNucleationProcess objects are not assignable.
  MAMGasAerosolExchangeProcess &operator=(const MAMGasAerosolExchangeProcess &) = delete;

  //------------------------------------------------------------------------
  //                                Accessors
  //------------------------------------------------------------------------

  void init(const ModalAerosolConfig& config) {
     init_(config);
  }

  KOKKOS_FUNCTION
  void run(Real t,
           Real dt, const Prognostics &prognostics,
           const Atmosphere &atmosphere, 
           const Diagnostics &diagnostics,
           Tendencies &tendencies) const {
     run_(t, dt, prognostics, atmosphere, diagnostics, tendencies);
  }

  /// Set the named parameter to the given value.
  /// It is a fatal error to pass an unknown name.
  virtual void set_param(const std::string &name, Real value) final {}

protected :

  virtual void init_(const ModalAerosolConfig& config) override;

  KOKKOS_FUNCTION
  virtual void run_(Real t,
                    Real dt, 
                    const Prognostics &prognostics,
                    const Atmosphere &atmosphere, 
                    const Diagnostics &diagnostics,
                    const Tendencies &tendencies) const override {}

};

}  // namespace haero

#endif
