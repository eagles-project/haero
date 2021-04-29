#ifndef HAERO_MODEL_HPP
#define HAERO_MODEL_HPP

#include "haero/modal_aerosol_config.hpp"
#include "haero/selected_processes.hpp"
#include "haero/process.hpp"
#include "haero/prognostics.hpp"
#include "haero/atmosphere.hpp"
#include "haero/diagnostics.hpp"
#include "haero/tendencies.hpp"
#include <map>

namespace haero {

/// @class Model
/// This type represents an aerosol system to be simulated, including all
/// information about modes, species, chemistry, and selected processes.
class Model final {
  public:

  /// Creates an aerosol model that supports the selected processes.
  /// @param [in] modal_aerosol_config The configuration of aerosols and gas
  ///                                  supported by the model.
  /// @param [in] selected_processes the set of aerosol processes supported by
  ///                                the model
  /// @param [in] num_levels The number of vertical levels in each column within
  ///                        the Context's computational domain
  Model(const ModalAerosolConfig& modal_aerosol_config,
        const SelectedProcesses& selected_processes,
        int num_levels);

  /// This factory function creates a model for use with unit tests. It
  /// initializes the Fortran subsystem, but doesn't perform any process
  /// initialization, and omits the selection of processes.
  static Model* ForUnitTests(const ModalAerosolConfig& modal_aerosol_config,
                             int num_levels);

  /// Copy constructor.
  Model(const Model&);

  /// Destructor.
  ~Model();

  /// Models are not assignable.
  Model& operator=(const Model&) = delete;

  /// Creates a new Prognostics object that can be used with this Model, given
  /// a set of views representing aerosol data managed by a host model. See
  /// the Prognostics class constructor for details on these views.
  Prognostics* create_prognostics(SpeciesColumnView int_aerosols,
                                  SpeciesColumnView cld_aerosols,
                                  SpeciesColumnView gases,
                                  ModalColumnView   modal_num_concs) const;

  /// Creates a new empty HostDiagnostics object that can be used with this Model.
  /// All fields within this new HostDiagnostics are owned and managed by it.
  HostDiagnostics* create_diagnostics() const;

  // Processes

  /// Runs the (prognostic) aerosol process of the given type, computing any
  /// relevant tendencies.
  /// @param [in] type The type of aerosol state to be updated.
  /// @param [in] t The time at which the process runs.
  /// @param [in] dt The time interval over which the process runs.
  /// @param [in] prognostics The prognostic variables used by this process.
  /// @param [in] atmosphere The atmospheric state variables used by this process.
  /// @param [in] diagnostics The diagnostic variables used by this process.
  /// @param [out] tendencies The aerosol tendencies computed.
  void run_process(ProcessType type,
                   Real t, Real dt,
                   const Prognostics& prognostics,
                   const Atmosphere& atmosphere,
                   const Diagnostics& diagnostics,
                   Tendencies& tendencies);

  /// Updates the state with the (diagnostic) aerosol process of the given type.
  /// @param [in] type The type of aerosol state to be updated.
  /// @param [in] t The time at which the process runs.
  /// @param [in] prognostics The prognostic variables used by this process.
  /// @param [in] atmosphere The atmospheric state variables used by this process.
  /// @param [inout] diagnostics The diagnostic variables used by and updated by
  ///                            this process.
  void update_diagnostics(ProcessType type, Real t,
                          const Prognostics& prognostics,
                          const Atmosphere& atmosphere,
                          Diagnostics& diagnostics);

  // Accessors

  /// Returns the modal aerosol configuration associated with this aerosol modeļ.
  KOKKOS_INLINE_FUNCTION
  const ModalAerosolConfig& modal_aerosol_config() const {
    return modal_aerosol_config_;
  }

  /// Returns the selected set of processes associated with this aerosol model.
  const SelectedProcesses& selected_processes() const;

  /// Returns the number of vertical levels in the model.
  int num_levels() const { return num_levels_; }

  /// Returns the number of modes in the model
  KOKKOS_INLINE_FUNCTION
  int num_modes() const {return modal_aerosol_config_.num_modes();}

  KOKKOS_INLINE_FUNCTION
  int num_gases() const {return modal_aerosol_config_.num_gases();}

  /// Returns the total number of distinct aerosol species populations
  /// (mode-species pairs).
  int num_aerosol_populations() const {
    return modal_aerosol_config_.num_aerosol_populations;
  }

  private:

  // Default constructor--used only internally
  Model();

  // This validates the Model's given parameters, throwing an exception if an
  // issue is encountered.
  void validate();

  // This gathers the set of processes based on selections, returning true if
  // any are backed by Fortran, and false if not.
  bool gather_processes();

  // This initializes the Model's Fortran subsystem.
  void init_fortran();

  // The modal aerosol configuration.
  ModalAerosolConfig modal_aerosol_config_;

  // Process selections.
  SelectedProcesses selected_processes_;

  // Grid parameters.
  int num_levels_;

  // Selected implementations of prognostic processes used by this model.
  std::map<ProcessType, PrognosticProcess*> prog_processes_;

  // Selected implementations of diagnostic processes used by this model.
  std::map<ProcessType, DiagnosticProcess*> diag_processes_;

  // This flag is set if this model initializes the haero Fortran model.
  bool uses_fortran_;
};

}

#endif
