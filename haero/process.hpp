#ifndef HAERO_PROCESS_HPP
#define HAERO_PROCESS_HPP

#include "haero/prognostics.hpp"
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
  virtual ~PrognosticProcess();

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
  /// @param [in] diagnostics The prognostic variables used by and affected by
  ///                         this process.
  /// @param [out] tendencies A container that stores time derivatives for
  ///                         prognostic variables evolved by this process.
  virtual void run(const Model& model,
                   Real t, Real dt,
                   const Prognostics& prognostics,
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
           const Diagnostics& diagnostics,
           Tendencies& tendencies) const override {}

};

/// @class FPrognosticProcess
/// This PrognosticProcess makes it easier to implement an aerosol process in
/// Fortran by wrapping the run method in a Fortran call that creates the
/// proper Fortran proxies for the model, prognostics, diagnostics, and
/// tendencies.
class FPrognosticProcess: public PrognosticProcess
{
  public:

  /// This type represents an interoperable Fortran subroutine that initializes
  /// a process using information from an aerosol model. It accepts an aerosol
  /// model and returns a pointer to an initialized Fortran process data type.
  /// Can return a null pointer if the process doesn't store any data.
  typedef void* (*init_subroutine_t)(void* model);

  /// This type represents an interoperable Fortran subroutine that runs
  /// a process with interoperable arguments.
  typedef void (*run_subroutine_t)(void* process, void* model,
                                   Real t, Real dt,
                                   void* prognostics, void* diagnostics,
                                   void* tendencies);

  /// This type represents an interoperable Fortran subroutine that destroys
  /// a process using information from an aerosol model. It accepts a pointer to
  /// a Fortran process data type.
  typedef void (*destroy_subroutine_t)(void* process);

  /// Constructor.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  /// @param [in] init_subroutine A pointer to an interoperable Fortran
  ///                             subroutine that initializes the process.
  /// @param [in] run_subroutine A pointer to an interoperable Fortran
  ///                            subroutine that runs the process.
  /// @param [in] destroy_subroutine A pointer to an interoperable Fortran
  ///                                subroutine that destroys the process.
  FPrognosticProcess(ProcessType type,
                     const std::string& name,
                     init_subroutine_t init_subroutine,
                     run_subroutine_t run_subroutine,
                     destroy_subroutine_t destroy_subroutine):
    PrognosticProcess(type, name), process_(nullptr),
    init_(init_subroutine), run_(run_subroutine), destroy_(destroy_subroutine)
  {
    EKAT_ASSERT_MSG((init_ != nullptr) && (destroy_ != nullptr),
      "init and destroy subroutines must both be specified or omitted.");
    EKAT_ASSERT(run_ != nullptr);
  }

  /// Destructor.
  ~FPrognosticProcess() {
    if (process_ != nullptr) {
      destroy_(process_);
    }
  }

  // Overrides.
  void init(const Model& model) override {
    EKAT_ASSERT(process_ == nullptr);
    if (init_ != nullptr) {
      process_ = init_((void*)&model);
    }
  }

  void run(const Model& model,
           Real t, Real dt,
           const Prognostics& prognostics,
           const Diagnostics& diagnostics,
           Tendencies& tendencies) const override {
    run_(process_, (void*)&model, t, dt,
         (void*)&prognostics, (void*)&diagnostics, (void*)&tendencies);
  }

  private:

  // Pointer to Fortran process data type. Can be null if the process stores
  // no data.
  void* process_;

  // Interoperable Fortran subroutines for initializing and running a process.
  init_subroutine_t init_;
  run_subroutine_t run_;
  destroy_subroutine_t destroy_;

};

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

  /// Destructor.
  virtual ~DiagnosticProcess();

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

  /// Override this method to update the given aerosol state at the given time.
  /// @param [in] model The aerosol model describing the aerosol system to
  ///                   which this process belongs.
  /// @param [in] t The simulation time at which this process is being invoked
  ///               (in seconds).
  /// @param [in] prognostics The prognostic variables used by this process.
  /// @param [inout] diagnostics The diagnostic variables used by and updated by
  ///                            this process.
  virtual void update(const Model& model, Real t,
                      const Prognostics& prognostics,
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
  /// @param [in] modal_vars A list of names of mode-ѕpecific variables required
  ///                        by the process
  void set_diag_vars(const std::vector<std::string>& vars,
                     const std::vector<std::string>& modal_vars) {
    required_vars_ = vars;
    required_modal_vars_ = modal_vars;
  }

  private:

  // The specific type of aerosol process.
  ProcessType type_;
  // The name of this diagnostic process.
  std::string name_;
  // The non-modal variables required by this process.
  std::vector<std::string> required_vars_;
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
              Diagnostics& diagnostics) const override {}
};

/// @class FDiagnosticProcess
/// This DiagnosticProcess makes it easier to implement an aerosol process in
/// Fortran by wrapping the update method in a Fortran call that creates the
/// proper Fortran proxies for the model, prognostics, and diagnostics.
class FDiagnosticProcess: public DiagnosticProcess {
  public:

  /// This type represents an interoperable Fortran subroutine that initializes
  /// a process using information from an aerosol model. It accepts an aerosol
  /// model and returns a pointer to an initialized Fortran process data type.
  /// Can return a null pointer if the process doesn't store any data.
  typedef void* (*init_subroutine_t)(void* model);

  /// This type represents an interoperable Fortran subroutine that updates
  /// diagnostic variables using interoperable arguments.
  typedef void (*update_subroutine_t)(void* process, void* model, Real t,
                                      void* prognostics, void* diagnostics);

  /// This type represents an interoperable Fortran subroutine that destroys
  /// a process using information from an aerosol model. It accepts a pointer to
  /// a Fortran process data type.
  typedef void (*destroy_subroutine_t)(void* process);

  /// Constructor.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  /// @param [in] init_subroutine A pointer to an interoperable Fortran
  ///                             subroutine that initializes the process.
  /// @param [in] update_subroutine A pointer to an interoperable Fortran
  ///                               subroutine that updates diagnostic variables.
  /// @param [in] destroy_subroutine A pointer to an interoperable Fortran
  ///                                subroutine that destroys the process.
  FDiagnosticProcess(ProcessType type,
                     const std::string& name,
                     init_subroutine_t init_subroutine,
                     update_subroutine_t update_subroutine,
                     destroy_subroutine_t destroy_subroutine):
    DiagnosticProcess(type, name), process_(nullptr),
    init_(init_subroutine), update_(update_subroutine),
    destroy_(destroy_subroutine)
  {
    EKAT_ASSERT_MSG((init_ != nullptr) && (destroy_ != nullptr),
      "init and destroy subroutines must both be specified or omitted.");
    EKAT_ASSERT(update_ != nullptr);
  }

  /// Destructor.
  ~FDiagnosticProcess() {
    if (process_ != nullptr) {
      destroy_(process_);
    }
  }

  // Overrides.
  void init(const Model& model) override {
    EKAT_ASSERT(process_ == nullptr);
    if (init_ != nullptr) {
      process_ = init_((void*)&model);
    }
  }

  void update(const Model& model, Real t,
              const Prognostics& prognostics,
              Diagnostics& diagnostics) const override {
    update_(process_, (void*)&model, t, (void*)&prognostics,
            (void*)&diagnostics);
  }

  private:

  // Pointer to Fortran process data type. Can be null if the process stores
  // no data.
  void* process_;

  // Interoperable Fortran subroutines for initializing and running a process.
  init_subroutine_t init_;
  update_subroutine_t update_;
  destroy_subroutine_t destroy_;

};

}

#endif
