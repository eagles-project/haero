#ifndef HAERO_MODEL_HPP
#define HAERO_MODEL_HPP

#include "haero/mode.hpp"
#include "haero/species.hpp"
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

  /// This is the device on which the Prognostics stores its data.
  using DeviceType = ekat::KokkosTypes<ekat::DefaultDevice>;

  /// This type represents vectorizable packs of Reals of length HAERO_PACK_SIZE.
  using PackType = ekat::Pack<Real, HAERO_PACK_SIZE>;

  /// Creates an aerosol model that supports the selected processes.
  /// @param [in] selected_processes the set of processes (including
  ///                                implementations) supported by the resulting
  ///                                Context
  /// @param [in] aerosol_modes a list of aerosol modes supported by the Context
  /// @param [in] aerosol_species a list of aerosol species supported by the
  ///                             Context
  /// @param [in] mode_species a map that defines the association of aerosol
  ///                          species with an aerosol mode. The keys in this
  ///                          map are names of aerosol modes (corresponding to
  ///                          those found in `modes`), and the values are lists
  ///                          of symbolic names (symbols) of aerosol species
  ///                          (supplied in `aerosol_species`) that belong to
  ///                          those modes.
  /// @param [in] gas_species a list of gas species supported by the Context
  /// @param [in] num_levels The number of vertical levels in each column within
  ///                        the Context's computational domain
  Model(const SelectedProcesses& selected_processes,
        const std::vector<Mode>& aerosol_modes,
        const std::vector<Species>& aerosol_species,
        const std::map<std::string, std::vector<std::string> >& mode_species,
        const std::vector<Species>& gas_species,
        int num_levels);

  /// This factory function creates a model for use with unit tests. It
  /// initializes the Fortran subsystem, but doesn't perform any process
  /// initialization, and omits the selection of processes.
  static Model* ForUnitTests(const std::vector<Mode>& aerosol_modes,
                             const std::vector<Species>& aerosol_species,
                             const std::map<std::string, std::vector<std::string> >& mode_species,
                             const std::vector<Species>& gas_species,
                             int num_levels);

  Model(const Model&);

  /// Destructor.
  ~Model();

  /// Models are not assignable either.
  Model& operator=(const Model&) = delete;

  /// Creates a new Prognostics object that can be used with this Model, given
  /// a set of views representing aerosol data managed by a host model. See
  /// the Prognostics class constructor for details on these views.
  using kokkos_device_type = ekat::KokkosTypes<ekat::DefaultDevice>;
  using SpeciesColumnView  = kokkos_device_type::view_2d<PackType>;
  using ModalColumnView    = kokkos_device_type::view_2d<PackType>;
  Prognostics* create_prognostics(SpeciesColumnView int_aerosols,
                                  SpeciesColumnView cld_aerosols,
                                  SpeciesColumnView gases,
                                  ModalColumnView   modal_num_concs) const;

  /// Creates a new empty DiagnosticsRegister object that can be used with this Model.
  /// All fields within this new DiagnosticsRegister are owned and managed by it.
  DiagnosticsRegister* create_diagnostics() const;

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

  /// Returns the selected set of processes associated with this aerosol model.
  const SelectedProcesses& selected_processes() const;

  /// Returns the list of aerosol modes associated with this aerosol model.
  const std::vector<Mode>& modes() const;

  /// Returns the list of all aerosol species associated with this aerosol
  /// model.
  const std::vector<Species>& aerosol_species() const;

  /// Returns the list of aerosol species associated with the model with the
  /// given mode index.
  /// @param [in]mode_index An integer index identifying the mode in question. This
  ///                       This index goes from 0 to num_modes-1.
  std::vector<Species> aerosol_species_for_mode(int mode_index) const;

  /// Returns the total number of distinct aerosol species populations in the
  /// model, counting appearances of one species in different modes separately.
  int num_aerosol_populations() const { return num_aero_populations_; }

  /// Returns the list of gas species associated with this aerosol model.
  const std::vector<Species>& gas_species() const;

  /// Returns the number of vertical levels in the model.
  int num_levels() const { return num_levels_; }

  private:

  // Default constructor--used only internally
  Model();

  // This validates the Model's given parameters, throwing an exception if an
  // issue is encountered.
  void validate();

  // This sets mode->species indexing. Throws an exception if the mode_species
  // mapping produces an inconsistent configuration.
  void index_modal_species(const std::map<std::string, std::vector<std::string> >& mode_species);

  // This gathers the set of processes based on selections, returning true if
  // any are backed by Fortran, and false if not.
  bool gather_processes();

  // This initializes the Model's Fortran subsystem.
  void init_fortran();

  // Process selections, modes, species...
  SelectedProcesses selected_processes_;
  std::vector<Mode> modes_;
  std::vector<Species> aero_species_;
  std::vector<Species> gas_species_;

  // The association of aerosol species with modes.
  // species_for_modes_[mode_name] = vector of species names
  std::vector<std::vector<int> > species_for_mode_;

  // The total number of modal aerosol populations.
  int num_aero_populations_;

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
