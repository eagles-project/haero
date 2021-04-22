#ifndef HAERO_AEROSOL_PROCESS_HPP
#define HAERO_AEROSOL_PROCESS_HPP

#include "haero/modal_aerosol_config.hpp"
#include "haero/prognostics.hpp"
#include "haero/atmosphere.hpp"
#include "haero/diagnostics.hpp"
#include "haero/tendencies.hpp"

namespace haero {

/// @enum AerosolProcessType
/// This enumerated type lists all relevant physical aerosol processes.
enum AerosolProcessType {
  ActivationProcess,
  CloudBorneWetRemovalProcess,
  CoagulationProcess,
  CondensationProcess,
  DryDepositionProcess,
  EmissionsProcess,
  InterstitialWetRemovalProcess,
  NucleationProcess,
  ResuspensionProcess,
  Terminator // This isn't a real process type--it only terminates the enum!
};

/// @class AerosolProcess
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
/// The AerosolProcess class provides an interface for all
/// implementations of all parametrizations for all physical processes that
/// compute tendencies for aerosol systems.
class AerosolProcess {
  public:

  /// Constructor, called by all AerosolProcess subclasses.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  AerosolProcess(AerosolProcessType type, const std::string& name):
    type_(type), name_(name) {}

  /// Destructor.
  KOKKOS_INLINE_FUNCTION
  virtual ~AerosolProcess() {}

  /// Default constructor is disabled.
  AerosolProcess() = delete;

  /// Default copy constructor. For use in moving host instance to device.
  KOKKOS_INLINE_FUNCTION
  AerosolProcess(const AerosolProcess& pp) :
    type_(pp.type_), name_(pp.name_) {}

  /// AerosolProcess objects are not assignable.
  AerosolProcess& operator=(const AerosolProcess&) = delete;

  //------------------------------------------------------------------------
  //                                Accessors
  //------------------------------------------------------------------------

  /// Returns the type of physical process modeled by this process.
  AerosolProcessType type() const { return type_; }

  /// Returns the name of this process/parametrization/realization.
  std::string name() const { return name_.label(); }

  //------------------------------------------------------------------------
  //                Methods to be overridden by subclasses
  //------------------------------------------------------------------------

  /// Override this method if your aerosol process needs to be initialized
  /// with information about the model. The default implementation does nothing.
  /// @param [in] modal_aerosol_config The aerosol configuration describing the
  ///                                  aerosol system to which this process
  ///                                  belongs.
  virtual void init(const ModalAerosolConfig& modal_aerosol_config) {}

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
                   Tendencies& tendencies) const = 0;

  private:

  const AerosolProcessType type_;
  // Use View as a struct to store a string and allows copy to device.
  // Since std::string can not be used, it was either this or a char *
  const Kokkos::View<int>  name_;
};

/// @class NullAerosolProcess
/// This AerosolProcess represents a null process in which no tendencies
/// are computed. It can be used for all prognostic processes that have been
/// disabled.
class NullAerosolProcess: public AerosolProcess {
  public:

  /// Constructor: constructs a null aerosol process of the given type.
  /// @param [in] type The type of aerosol process.
  explicit NullAerosolProcess(AerosolProcessType type):
    AerosolProcess(type, "Null prognostic aerosol process") {}

  // Overrides
  KOKKOS_FUNCTION
  void run(const ModalAerosolConfig& modal_aerosol_config,
           Real t, Real dt,
           const Prognostics& prognostics,
           const Atmosphere& atmosphere,
           const Diagnostics& diagnostics,
           Tendencies& tendencies) const override {}

};

#if HAERO_FORTRAN

/// @class FAerosolProcess
/// This AerosolProcess makes it easier to implement an aerosol process in
/// Fortran by wrapping the run method in a Fortran call that creates the
/// proper Fortran proxies for the model, prognostics, diagnostics, and
/// tendencies.
class FAerosolProcess: public AerosolProcess
{
  public:

  /// This type is a pointer to a Fortran function that initializes a process.
  /// Since the model is already available to the Fortran process via Haero's
  /// Fortran helper module, this type of function takes no arguments.
  typedef void (*InitProcessFunction)(void);

  /// This type is a pointer to a Fortran function that runs a process.
  typedef void (*RunProcessFunction)(Real t, Real dt, void* progs, void* atm,
                                     void* diags, void* tends);

  /// This type is a pointer to a Fortran function that finalizes a process,
  /// deallocating any resources allocated in the process's init function.
  typedef void (*FinalizeProcessFunction)(void);

  /// Constructor.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  /// @param [in] create_process A pointer to an interoperable Fortran function
  ///                            that returns a C pointer to a newly allocated
  ///                            Fortran-backed prognostic process.
  FAerosolProcess(AerosolProcessType type,
                  const std::string& name,
                  InitProcessFunction init_process,
                  RunProcessFunction run_process,
                  FinalizeProcessFunction finalize_process):
    AerosolProcess(type, name), init_process_(init_process),
    run_process_(run_process), finalize_process_(finalize_process),
    initialized_(false) {}

  /// Copy constructor.
  FAerosolProcess(const FAerosolProcess& pp) :
    AerosolProcess(pp), init_process_(pp.init_process_),
    run_process_(pp.run_process_), finalize_process_(pp.finalize_process_),
    initialized_(false) {}

  /// Destructor.
  ~FAerosolProcess() {
    if (initialized_) {
      finalize_process_();
      initialized_ = false;
    }
  }

  // Overrides.
  void init(const ModalAerosolConfig& modal_aerosol_config) override {
    if (not initialized_) {
      init_process_();
      initialized_ = true;
    }
  }

  void run(const ModalAerosolConfig& modal_aerosol_config,
           Real t, Real dt,
           const Prognostics& prognostics,
           const Atmosphere& atmosphere,
           const Diagnostics& diagnostics,
           Tendencies& tendencies) const override {
    // Set tendencies to zero.
    tendencies.scale(0.0);

    // Call the Fortran subroutine for the process.
    run_process_(t, dt, (void*)&prognostics, (void*)&atmosphere,
                 (void*)&diagnostics, (void*)&tendencies);
  }

  private:

  // Pointers to Fortran subroutines.
  InitProcessFunction init_process_;
  RunProcessFunction run_process_;
  FinalizeProcessFunction finalize_process_;

  // Has the process been initialized?
  bool initialized_;
};

#endif // HAERO_FORTRAN

} // namespace haero

#endif
