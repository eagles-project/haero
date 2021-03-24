#ifndef HAERO_MAM_NUCLEATION_PROCESS_HPP
#define HAERO_MAM_NUCLEATION_PROCESS_HPP

#include "haero/process.hpp"


namespace haero {


/// @class MAMNucleationProcess
/// This type is an interface (base class) to an aerosol process quantified by
/// a parametrization that computes a set of tendencies for prognostic variables
/// within an aerosol system. Each subclass of this type implements a particular
/// **implementation** of a specific **parametrization** for a particular **process**.
///
/// To make these ideas more complete, consider the following examples of
/// important **physical processes** in the ProcessType above.
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
/// The PrognosticProcess class provides an interface for all
/// implementations of all parametrizations for all physical processes that
/// compute tendencies for aerosol systems.
class MAMNucleationProcess : public PrognosticProcess {
  public:

  /// Constructor, called by all PrognosticProcess subclasses.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  MAMNucleationProcess();

  /// Destructor.
  KOKKOS_INLINE_FUNCTION
  virtual ~MAMNucleationProcess() {}

  /// Default copy constructor. For use in moving host instance to device.
  KOKKOS_INLINE_FUNCTION
  MAMNucleationProcess(const MAMNucleationProcess& pp) :
    PrognosticProcess(pp) {}

  /// MAMNucleationProcess objects are not assignable.
  PrognosticProcess& operator=(const MAMNucleationProcess&) = delete;

  //------------------------------------------------------------------------
  //                                Accessors
  //------------------------------------------------------------------------


  //------------------------------------------------------------------------
  //                Methods to be overridden by subclasses
  //------------------------------------------------------------------------

  /// Override this method if your aerosol process needs to be initialized
  /// with information about the model. The default implementation does nothing.
  /// @param [in] modal_aerosol_config The aerosol configuration describing the
  ///                                  aerosol system to which this process
  ///                                  belongs.
  virtual void init(const ModalAerosolConfig& modal_aerosol_config);

  /// Override this method to implement the aerosol process using the specific
  /// parameterization for the subclass.
  /// @param [in] modal_aerosol_config The aerosol configuration describing the
  ///                                  aerosol system to which this process
  ///                                  belongs.
  /// @param [in] t The simulation time at which this process is being invoked
  ///               (in seconds).
  /// @param [in] dt The simulation time interval ("timestep size") over which
  ///                this process occurs.
  /// @param [in] prognostics The prognostic variables used by and affected by
  ///                         this process.
  /// @param [in] atmosphere The atmosphere state variables used by this process.
  /// @param [in] diagnostics The prognostic variables used by and affected by
  ///                         this process.
  /// @param [out] tendencies A container that stores time derivatives for
  ///                         prognostic variables evolved by this process.
  KOKKOS_FUNCTION
  virtual void run(const ModalAerosolConfig& modal_aerosol_config,
                   Real t, Real dt,
                   const Prognostics& prognostics,
                   const Atmosphere& atmosphere,
                   const Diagnostics& diagnostics,
                   Tendencies& tendencies) const;


  KOKKOS_FUNCTION 
  static void ternary_nuc_merik2007(const double t, 
                                    const double rh, 
                                    const double c2,
                                    const double c3, 
                                    double &j_log, 
                                    double &ntot, 
                                    double &nacid, 
                                    double &namm, 
                                    double &r);

  
};


}

#endif
