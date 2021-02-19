#ifndef HAERO_PROCESS_HPP
#define HAERO_PROCESS_HPP

#include "haero/prognostics.hpp"
#include "haero/atmosphere.hpp"
#include "haero/diagnostics.hpp"
#include "haero/tendencies.hpp"

namespace haero {

/// This is a forward declaration, since we refer to Model below.
class Model;

/// @enum ProcessType
/// This enumerated type lists all relevant physical aerosol processes.
enum ProcessType {
  ActivationProcess,
  CloudBorneWetRemovalProcess,
  CoagulationProcess,
  CondensationProcess,
  DryDepositionProcess,
  EmissionsProcess,
  InterstitialWetRemovalProcess,
  NucleationProcess,
  ResuspensionProcess,
  WaterUptakeProcess
};

/// @class PrognosticProcess
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
class PrognosticProcess {
  public:

  /// Constructor, called by all PrognosticProcess subclasses.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  PrognosticProcess(ProcessType type, const std::string& name):
    type_(type), name_(name) {}

  /// Destructor.
  virtual ~PrognosticProcess() {}

  /// Default constructor is disabled.
  PrognosticProcess() = delete;

  /// PrognosticProcess objects are not deep-copyable. They should be passed
  /// by reference or as pointers.
  PrognosticProcess(const PrognosticProcess&) = delete;

  /// PrognosticProcess objects are not assignable either.
  PrognosticProcess& operator=(const PrognosticProcess&) = delete;

  //------------------------------------------------------------------------
  //                                Accessors
  //------------------------------------------------------------------------

  /// Returns the type of physical process modeled by this process.
  ProcessType type() const { return type_; }

  /// Returns the name of this process/parametrization/realization.
  const std::string& name() const { return name_; }

  //------------------------------------------------------------------------
  //                Methods to be overridden by subclasses
  //------------------------------------------------------------------------

  /// Override this method if your aerosol process needs to be initialized
  /// with information about the model. The default implementation does nothing.
  /// @param [in] model The aerosol model describing the aerosol system to
  ///                   which this process belongs.
  virtual void init(const Model& model) {}

  /// Override this method to implement the aerosol process using the specific
  /// parameterization for the subclass.
  /// @param [in] model The aerosol model describing the aerosol system to
  ///                   which this process belongs.
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
  virtual void run(const Model& model,
                   Real t, Real dt,
                   const Prognostics& prognostics,
                   const Atmosphere& atmosphere,
                   const Diagnostics& diagnostics,
                   Tendencies& tendencies) const = 0;

  private:

  ProcessType type_;
  std::string name_;
};

/// @class NullPrognosticProcess
/// This PrognosticProcess represents a null process in which no tendencies
/// are computed. It can be used for all prognostic processes that have been
/// disabled.
class NullPrognosticProcess: public PrognosticProcess {
  public:

  /// Constructor: constructs a null aerosol process of the given type.
  /// @param [in] type The type of aerosol process.
  explicit NullPrognosticProcess(ProcessType type):
    PrognosticProcess(type, "Null prognostic aerosol process") {}

  // Overrides
  void run(const Model& model,
           Real t, Real dt,
           const Prognostics& prognostics,
           const Atmosphere& atmosphere,
           const Diagnostics& diagnostics,
           Tendencies& tendencies) const override {}

};

#if HAERO_FORTRAN

/// @class FPrognosticProcess
/// This PrognosticProcess makes it easier to implement an aerosol process in
/// Fortran by wrapping the run method in a Fortran call that creates the
/// proper Fortran proxies for the model, prognostics, diagnostics, and
/// tendencies.
class FPrognosticProcess: public PrognosticProcess
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
  FPrognosticProcess(ProcessType type,
                     const std::string& name,
                     InitProcessFunction init_process,
                     RunProcessFunction run_process,
                     FinalizeProcessFunction finalize_process):
    PrognosticProcess(type, name), init_process_(init_process),
    run_process_(run_process), finalize_process_(finalize_process),
    initialized_(false) {}

  /// Destructor.
  ~FPrognosticProcess() {
    if (initialized_) {
      finalize_process_();
    }
  }

  // Overrides.
  void init(const Model& model) override {
    init_process_();
    initialized_ = true;
  }

  void run(const Model& model,
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

/// @def HAERO_CREATE_PROGNOSTIC_FORTRAN_BRIDGE(module_name)
/// This macro generates C declarations for three interoperable Fortran
/// subroutines that correspond to the init, run, and finalize behaviors
/// for an FPrognosticProcess.
#define CREATE_PROGNOSTIC_FORTRAN_BRIDGE(module_name) \
extern "C" { \
extern void module_name##_bridge_init(void); \
extern void module_name##_bridge_run(Real, Real, void*, void*, void*, void*); \
extern void module_name##_bridge_finalize(void); \
}

/// @def PROGNOSTIC_FORTRAN_BRIDGE(module_name)
/// Use this in place of the last three arguments to the FPrognosticProcess
/// base class constructor. You must use @ref CREATE_PROGNOSTIC_FORTRAN_BRIDGE
/// to declare the interfaces for the interoperable bridge functions first.
#define PROGNOSTIC_FORTRAN_BRIDGE(module_name) \
  module_name##_bridge_init, \
  module_name##_bridge_run, \
  module_name##_bridge_finalize

#endif // HAERO_FORTRAN

/// @class DiagnosticProcess
/// This type is an interface (base class) to an aerosol process quantified by
/// a parametrization that updates the **diagnostic** variables in an aerosol
/// system.
///
class DiagnosticProcess {
  public:

  /// Constructor, called by all DiagnosticProcess subclasses.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  DiagnosticProcess(ProcessType type, const std::string& name):
    type_(type), name_(name) {}

  /// Constructor that accepts names of all modal and non-modal diagnostic
  /// variables needed by this process.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  /// @param [in] variables The set of (non-modal) diagnostic variables
  ///                       required by this process.
  /// @param [in] aero_variables The set of per-aerosol-species diagnostic
  ///                            variables required by this process.
  /// @param [in] gas_variables The set of per-gas-species diagnostic
  ///                            variables required by this process.
  /// @param [in] modal_variables The set of modal diagnostic variables
  ///                             required by this process.
  DiagnosticProcess(ProcessType type, const std::string& name,
                    const std::vector<std::string>& variables,
                    const std::vector<std::string>& aero_variables,
                    const std::vector<std::string>& gas_variables,
                    const std::vector<std::string>& modal_variables):
    type_(type), name_(name) {
    set_diag_vars(variables, aero_variables, gas_variables, modal_variables);
  }

  /// Destructor.
  virtual ~DiagnosticProcess() {}

  /// Default constructor is disabled.
  DiagnosticProcess() = delete;

  /// DiagnosticProcess objects are not deep-copyable. They should be
  /// passed by reference or as pointers.
  DiagnosticProcess(const DiagnosticProcess&) = delete;

  /// DiagnosticProcess objects are not assignable either.
  DiagnosticProcess& operator=(const DiagnosticProcess&) = delete;

  //------------------------------------------------------------------------
  //                                Accessors
  //------------------------------------------------------------------------

  /// Returns the type of physical process modeled by this process.
  ProcessType type() const { return type_; }

  /// Returns the name of this process/parametrization/realization.
  const std::string& name() const { return name_; }

  //------------------------------------------------------------------------
  //                Methods called during initialization
  //------------------------------------------------------------------------

  /// This method creates all necessary diagnostic variables within the
  /// given Diagnostics object. It's called during the model initialization
  /// process at the beginning of a simulation.
  /// @param [inout] diagnostics The Diagnostics object in which the needed
  ///                            variables are created.
  void prepare(Diagnostics& diagnostics) const {
    for (const auto& var: required_vars_) {
      if (not diagnostics.has_var(var)) {
        diagnostics.create_var(var);
      }
    }
    for (const auto& var: required_aero_vars_) {
      if (not diagnostics.has_aerosol_var(var)) {
        diagnostics.create_aerosol_var(var);
      }
    }
    for (const auto& var: required_gas_vars_) {
      if (not diagnostics.has_gas_var(var)) {
        diagnostics.create_gas_var(var);
      }
    }
    for (const auto& var: required_modal_vars_) {
      if (not diagnostics.has_modal_var(var)) {
        diagnostics.create_modal_var(var);
      }
    }
  }

  //------------------------------------------------------------------------
  //                Methods to be overridden by subclasses
  //------------------------------------------------------------------------

  /// Override this method if your aerosol process needs to be initialized
  /// with information about the model. The default implementation does nothing.
  /// @param [in] model The aerosol model describing the aerosol system to
  ///                   which this process belongs.
  virtual void init(const Model& model) {}

  /// Override this method to update diagnostic variables at the given time.
  /// @param [in] model The aerosol model describing the aerosol system to
  ///                   which this process belongs.
  /// @param [in] t The simulation time at which this process is being invoked
  ///               (in seconds).
  /// @param [in] prognostics The prognostic variables used by this process.
  /// @param [in] atmosphere The atmosphere state variables used by this process.
  /// @param [inout] diagnostics The diagnostic variables used by and updated by
  ///                            this process.
  virtual void update(const Model& model, Real t,
                      const Prognostics& prognostics,
                      const Atmosphere& atmosphere,
                      Diagnostics& diagnostics) const = 0;

  protected:

  //------------------------------------------------------------------------
  //                Methods for constructing diagnostic processes
  //------------------------------------------------------------------------

  /// This method accepts the names of diagnostic variables required by this
  /// process. These variables include those variables used by the process, as
  /// well as variables computed by the process. The process ensures that any
  /// Diagnostics object passed to it has storage for these variables. This
  /// is usually called inside the constructor of a DiagnosticProcess subclass.
  /// @param [in] vars A list of names of (non-mode-ѕpecific variables required
  ///                  by the process
  /// @param [in] aero_vars The set of per-aerosol-species diagnostic
  ///                       variables required by this process.
  /// @param [in] gas_vars The set of per-gas-species diagnostic
  ///                      variables required by this process.
  /// @param [in] modal_vars A list of names of mode-ѕpecific variables required
  ///                        by the process
  void set_diag_vars(const std::vector<std::string>& vars,
                     const std::vector<std::string>& aero_vars,
                     const std::vector<std::string>& gas_vars,
                     const std::vector<std::string>& modal_vars) {
    required_vars_ = vars;
    required_aero_vars_ = aero_vars;
    required_gas_vars_ = gas_vars;
    required_modal_vars_ = modal_vars;
  }

  private:

  // The specific type of aerosol process.
  ProcessType type_;
  // The name of this diagnostic process.
  std::string name_;
  // The non-modal variables required by this process.
  std::vector<std::string> required_vars_;
  // The per-aerosol-species variables required by this process.
  std::vector<std::string> required_aero_vars_;
  // The per-gas-species variables required by this process.
  std::vector<std::string> required_gas_vars_;
  // The modal variables required by this process.
  std::vector<std::string> required_modal_vars_;
};

/// @class NullDiagnosticProcess
/// This DiagnosticProcess represents a null process in which no updates
/// are performed. It can be used for all processes that
/// have been disabled.
class NullDiagnosticProcess: public DiagnosticProcess {
  public:

  /// Constructor: constructs a null aerosol process of the given type.
  /// @param [in] type The type of aerosol process.
  explicit NullDiagnosticProcess(ProcessType type):
    DiagnosticProcess(type, "Null diagnostic aerosol process") {}

  // Overrides
  void update(const Model& model, Real t,
              const Prognostics& prognostics,
              const Atmosphere& atmosphere,
              Diagnostics& diagnostics) const override {}
};

#if HAERO_FORTRAN
extern "C" {

/// Given a C pointer to a Fortran-backed diagnostic process, this function
/// initializes the process with the given aerosol model.
void fortran_diagnostic_process_init(void* process, void* model);

/// Given a C pointer to a Fortran-backed diagnostic process, this function
/// invokes the process to update diagnostic variables with the given arguments.
void fortran_diagnostic_process_update(void* process, void* model, Real t,
                                       void* prognostics, void* atmosphere,
                                       void* diagnostics);

/// Given a C pointer to a Fortran-backed diagnostic process, this function
/// frees all resources allocated within the process
void fortran_diagnostic_process_destroy(void* process);

} // extern "C"

/// @class FDiagnosticProcess
/// This DiagnosticProcess makes it easier to implement an aerosol process in
/// Fortran by wrapping the update method in a Fortran call that creates the
/// proper Fortran proxies for the model, prognostics, and diagnostics.
class FDiagnosticProcess: public DiagnosticProcess {
  public:

  /// This type is a pointer to a Fortran function that initializes a process.
  /// Since the model is already available to the Fortran process via Haero's
  /// Fortran helper module, this type of function takes no arguments.
  typedef void (*InitProcessFunction)(void);

  /// This type is a pointer to a Fortran function that runs a process update.
  typedef void (*UpdateProcessFunction)(Real t, void* progs, void* atm,
                                        void* diags);

  /// This type is a pointer to a Fortran function that finalizes a process,
  /// deallocating any resources allocated in the process's init function.
  typedef void (*FinalizeProcessFunction)(void);

  /// Constructor.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  /// @param [in] variables The set of (non-modal) diagnostic variables
  ///                       required by this process.
  /// @param [in] aero_variables The set of per-aerosol-species diagnostic
  ///                            variables required by this process.
  /// @param [in] gas_variables The set of per-gas-species diagnostic
  ///                            variables required by this process.
  /// @param [in] modal_variables The set of modal diagnostic variables
  ///                             required by this process.
  /// @param [in] create_process A pointer to an interoperable Fortran function
  ///                            that returns a C pointer to a newly allocated
  ///                            Fortran-backed diagnostic process.
  FDiagnosticProcess(ProcessType type,
                     const std::string& name,
                     const std::vector<std::string>& variables,
                     const std::vector<std::string>& aero_variables,
                     const std::vector<std::string>& gas_variables,
                     const std::vector<std::string>& modal_variables,
                     InitProcessFunction init_process,
                     UpdateProcessFunction update_process,
                     FinalizeProcessFunction finalize_process):
    DiagnosticProcess(type, name, variables, aero_variables, gas_variables,
                      modal_variables),
    init_process_(init_process), update_process_(update_process),
    finalize_process_(finalize_process), initialized_(false) {}

  /// Destructor.
  ~FDiagnosticProcess() {
    if (initialized_) {
      finalize_process_();
    }
  }

  // Overrides.
  void init(const Model& model) override {
    init_process_();
    initialized_ = true;
  }

  void update(const Model& model, Real t,
              const Prognostics& prognostics,
              const Atmosphere& atmosphere,
              Diagnostics& diagnostics) const override {
    update_process_(t, (void*)&prognostics, (void*)&atmosphere,
                    (void*)&diagnostics);
  }

  private:

  // Pointers to Fortran subroutines.
  InitProcessFunction init_process_;
  UpdateProcessFunction update_process_;
  FinalizeProcessFunction finalize_process_;

  // Has the process been initialized?
  bool initialized_;
};

/// @def HAERO_CREATE_DIAGNOSTIC_FORTRAN_BRIDGE(module_name)
/// This macro generates C declarations for three interoperable Fortran
/// subroutines that correspond to the init, update, and finalize behaviors
/// for an FDiagnosticProcess.
#define CREATE_DIAGNOSTIC_FORTRAN_BRIDGE(module_name) \
extern "C" { \
extern void module_name##_bridge_init(void); \
extern void module_name##_bridge_update(Real, void*, void*, void*); \
extern void module_name##_bridge_finalize(void); \
}

/// @def DIAGNOSTIC_FORTRAN_BRIDGE(module_name)
/// Use this in place of the last three arguments to the FDiagnosticProcess
/// base class constructor. You must use @ref CREATE_DIAGNOSTIC_FORTRAN_BRIDGE
/// to declare the interfaces for the interoperable bridge functions first.
#define DIAGNOSTIC_FORTRAN_BRIDGE(module_name) \
  module_name##_bridge_init, \
  module_name##_bridge_update, \
  module_name##_bridge_finalize

#endif // HAERO_FORTRAN

}

#endif
