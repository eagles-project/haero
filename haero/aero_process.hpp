#ifndef HAERO_AERO_PROCESS_HPP
#define HAERO_AERO_PROCESS_HPP

#include "haero/aero_state.hpp"
#include "haero/aero_tendencies.hpp"

namespace haero {

/// This is a forward declaration, since we refer to AeroModel below.
class AeroModel;

/// @enum AeroProcessType
/// This enumerated type lists all relevant physical aerosol processes.
enum AeroProcessType {
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

/// @class PrognosticAeroProceess
/// This type is an interface (base class) to an aerosol process quantified by
/// a parametrization that computes a set of tendencies for prognostic variables
/// within an aerosol system. Each subclass of this type implements a particular
/// **implementation** of a specific **parametrization** for a particular **process**.
///
/// To make these ideas more complete, consider the following examples of
/// important **physical processes** in the AeroProcessType above.
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
/// The PrognosticAeroProcess class provides an interface for all
/// implementations of all parametrizations for all physical processes that
/// compute tendencies for aerosol systems.
class PrognosticAeroProcess {
  public:

  /// Constructor, called by all PrognosticAeroProcess subclasses.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  PrognosticAeroProcess(AeroProcessType type, const std::string& name):
    type_(type), name_(name) {}

  /// Destructor.
  virtual ~PrognosticAeroProcess();

  /// Default constructor is disabled.
  PrognosticAeroProcess() = delete;

  /// PrognosticAeroProcess objects are not deep-copyable. They should be passed
  /// by reference or as pointers.
  PrognosticAeroProcess(const PrognosticAeroProcess&) = delete;

  /// PrognosticAeroProcess objects are not assignable either.
  PrognosticAeroProcess& operator=(const PrognosticAeroProcess&) = delete;

  //------------------------------------------------------------------------
  //                                Accessors
  //------------------------------------------------------------------------

  /// Returns the type of physical process modeled by this AeroProcess.
  AeroProcessType type() const { return type_; }

  /// Returns the name of this process/parametrization/realization.
  const std::string& name() const { return name_; }

  //------------------------------------------------------------------------
  //                Methods to be overridden by subclasses
  //------------------------------------------------------------------------

  /// Override this method to implement the aerosol process using the specific
  /// parameterization for the subclass.
  /// @param [in] model The aerosol model describing the aerosol system to
  ///                   which this process belongs.
  /// @param [in] t The simulation time at which this process is being invoked
  ///               (in seconds).
  /// @param [in] dt The simulation time interval ("timestep size") over which
  ///                this process occurs.
  /// @param [in] state The aerosol state containing the prognostic variables
  ///                   used by and affected by this process.
  /// @param [out] tendencies A container that stores time derivatives for
  ///                         prognostic variables evolved by this process.
  virtual void run(const AeroModel& model,
                   Real t, Real dt,
                   const AeroState& state,
                   AeroTendencies& tendencies) const = 0;

  private:

  AeroProcessType type_;
  std::string name_;
};

/// @class NullPrognosticAeroProcess
/// This PrognosticAeroProcess represents a null process in which no tendencies
/// are computed. It can be used for all prognostic processes that have been
/// disabled.
class NullPrognosticAeroProcess: public PrognosticAeroProcess {
public:
  /// Constructor: constructs a null aerosol process of the given type.
  /// @param [in] type The type of aerosol process.
  explicit NullPrognosticAeroProcess(AeroProcessType type):
    PrognosticAeroProcess(type, "Null prognostic aerosol process") {}

  // Overrides
  void run(const AeroModel& model,
           Real t, Real dt,
           const AeroState& state,
           AeroTendencies& tendencies) const override {}
};

/// @class DiagnosticAeroProcess
/// This type is an interface (base class) to an aerosol process quantified by
/// a parametrization that updates the **diagnostic** variables in an aerosol
/// system.
///
class DiagnosticAeroProcess {
  public:

  /// Constructor, called by all DiagnosticAeroProcess subclasses.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  DiagnosticAeroProcess(AeroProcessType type, const std::string& name):
    type_(type), name_(name) {}

  /// Destructor.
  virtual ~DiagnosticAeroProcess();

  /// Default constructor is disabled.
  DiagnosticAeroProcess() = delete;

  /// DiagnosticAeroProcess objects  are not deep-copyable. They should be
  /// passed by reference or as pointers.
  DiagnosticAeroProcess(const DiagnosticAeroProcess&) = delete;

  /// DiagnosticAeroProcess objects are not assignable either.
  DiagnosticAeroProcess& operator=(const DiagnosticAeroProcess&) = delete;

  //------------------------------------------------------------------------
  //                                Accessors
  //------------------------------------------------------------------------

  /// Returns the type of physical process modeled by this AeroProcess.
  AeroProcessType type() const { return type_; }

  /// Returns the name of this process/parametrization/realization.
  const std::string& name() const { return name_; }

  //------------------------------------------------------------------------
  //                Methods to be overridden by subclasses
  //------------------------------------------------------------------------

  /// Override this method to update the given aerosol state at the given time.
  /// @param [in] model The aerosol model describing the aerosol system to
  ///                   which this process belongs.
  /// @param [in] t The simulation time at which this process is being invoked
  ///               (in seconds).
  /// @param [in] state The aerosol state containing the prognostic variables
  ///                   used by and affected by this process.
  virtual void update(const AeroModel& model, Real t, AeroState& state) const = 0;

  private:

  AeroProcessType type_;
  std::string name_;
};

/// @class NullDiagnosticAeroProcess
/// This DiagnosticAeroProcess represents a null process in which no updates
/// are performed. It can be used for all processes that
/// have been disabled.
class NullDiagnosticAeroProcess: public DiagnosticAeroProcess {
public:
  /// Constructor: constructs a null aerosol process of the given type.
  /// @param [in] type The type of aerosol process.
  explicit NullDiagnosticAeroProcess(AeroProcessType type):
   DiagnosticAeroProcess(type, "Null diagnostic aerosol process") {}

  // Overrides
  void update(const AeroModel& model, Real t, AeroState& state) const override {}
};

}

#endif
