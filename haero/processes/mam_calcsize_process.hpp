#ifndef HAERO_MAM_CALCSIZE_PROCESS_HPP
#define HAERO_MAM_CALCSIZE_PROCESS_HPP

#include <iomanip>

#include "haero/aerosol_process.hpp"
#include "haero/physical_constants.hpp"

namespace haero {

/// @class MAMCalcsizeProcess
/// This type is an interface (base class) to an aerosol process quantified by
/// a parametrization that computes a set of tendencies for prognostic variables
/// within an aerosol system. Each subclass of this type implements a particular
/// **implementation** of a specific **parametrization** for a particular
/// **process**.
///
/// To make these ideas more complete, consider the following examples of
/// important **physical processes** in the AerosolProcessType above.
///
/// Each of these processes has one or more **parametrizations**--mathematical
/// models that quantify the outcomes of these processes in specific
/// circumstances. For example, **surface emissions** may be parametrized by:
/// * time series input from a file
/// * a low-order polynomial in the time variable
/// * estimated from, e.g., an agent-based representation of human activity.
///
/// Finally, each of these parametrizations can have one or more
/// **implementations**. For example, every parametrization can have a Fortran
/// implementation that runs only on CPUs, as well as a C++ implementation that
/// can run on CPUs and GPUs.
///
/// The AerosolProcess class provides an interface for all
/// implementations of all parametrizations for all physical processes that
/// compute tendencies for aerosol systems.
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

  virtual void init(const ModalAerosolConfig &modal_aerosol_config) override;

  KOKKOS_FUNCTION
  virtual void run(const ModalAerosolConfig &modal_aerosol_config, Real t,
                   Real dt, const Prognostics &prognostics,
                   const Atmosphere &atmosphere, const Diagnostics &diagnostics,
                   Tendencies &tendencies) const override{};

  /// Set the named parameter to the given value.
  /// It is a fatal error to pass an unknown name.
  virtual void set_param(const std::string &name, Real value) final {
  }
};

}  // namespace haero

#endif
