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

/// @class PrognosticProceess
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

  private:

  ProcessType type_;
  std::string name_;
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

}

#endif
