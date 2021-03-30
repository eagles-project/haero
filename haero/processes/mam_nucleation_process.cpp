
#include <cmath>

#include "haero/processes/mam_nucleation_process.hpp"


namespace haero {


  MAMNucleationProcess::MAMNucleationProcess():
    PrognosticProcess(NucleationProcess, "MAMNucleationProcess") {}


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
  void MAMNucleationProcess::init(const ModalAerosolConfig& modal_aerosol_config) {}

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
  void MAMNucleationProcess::run(const ModalAerosolConfig& modal_aerosol_config,
                   Real t, Real dt,
                   const Prognostics& prognostics,
                   const Atmosphere& atmosphere,
                   const Diagnostics& diagnostics,
                   Tendencies& tendencies) const {}

}
