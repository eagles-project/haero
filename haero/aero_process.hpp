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
  NucleationProcess,
  ResuspensionProcess,
  WaterUptakeProcess
};

/// @class AeroProcess
/// This type is an abstract interface (base class) to an aerosol process
/// quantified by a parametrization. Each subclass of this process implements a
/// particular **realization** of a specific **parametrization** for a
/// particular **process**.
///
/// To make these ideas more complete, consider the following examples of
/// important **physical processes** in aerosol systems:
/// * activation
/// * condensation
/// * emissions
/// * nucleation
/// * resuspension
/// * water uptake
/// The \ref AeroProcessType enum above includes these and other important aerosol
/// processes.
///
/// Each of these processes has one or more **parametrizations**--mathematical
/// models that quantify the outcomes of these processes in specific
/// circumstances. For example, **surface emissions** may be parametrized by:
/// * time series input from a file
/// * a low-order polynomial in the time variable
/// * estimated from, e.g., an agent-based representation of human activity.
///
/// Finally, each of these parametrizations can have one or more
/// **realizations** (or **implementations**). For example, every
/// parametrization can have a Fortran implementation that runs only on CPUs, as
/// well as a C++ implementation that can run on CPUs and GPUs.
///
/// The AeroProcess class provides an interface for all realizations of all
/// parametrizations for all physical processes relevant to aerosols.
class AeroProcess {
  public:

  /// Constructor, called by all AeroProcess subclasses.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its realization.
  AeroProcess(AeroProcessType type, const std::string& name):
    type_(type), name_(name) {}

  /// Destructor.
  virtual ~AeroProcess();

  /// Default constructor is disabled.
  AeroProcess() = delete;

  /// AeroProcesses are not deep-copyable. They should be passed by reference
  /// or as pointers.
  AeroProcess(const AeroProcess&) = delete;

  /// AeroProcesses are not assignable either.
  AeroProcess& operator=(const AeroProcess&) = delete;

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
  /// @param [in] state The aerosol state containing the prognostic variables
  ///                   used by and affected by this process.
  /// @param [in] t The simulation time at which this process is being invoked
  ///               (in seconds).
  /// @param [in] dt The simulation time interval ("timestep size") over which
  ///                this process occurs.
  /// @param [out] tendencies A container that stores time derivatives for
  ///                         prognostic variables evolved by this process.
  virtual void compute(const AeroModel& model,
                       const AeroState& state,
                       Real t,
                       Real dt,
                       AeroTendencies& tendencies) = 0;

  private:

  AeroProcessType type_;
  std::string name_;
};

}

#endif
