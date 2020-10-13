#ifndef HAERO_AERO_MODEL_HPP
#define HAERO_AERO_MODEL_HPP

#include "haero/mode.hpp"
#include "haero/species.hpp"
#include "haero/chemistry.hpp"
#include "haero/parameterizations.hpp"
#include "haero/aero_process.hpp"
#include "haero/aero_state.hpp"
#include "haero/aero_tendencies.hpp"
#include <map>

namespace haero {

/// @class AeroModel
/// This type represents an aerosol system to be simulated, including all
/// information about modes, species, chemistry, and selected parametrizations.
class AeroModel final {
  public:

  /// Creates a AeroModel that supports the specified parameterizations.
  /// @param [in] parameterizations the set of parameterizations (including
  ///                               implementations) supported by the resulting
  ///                               Context
  /// @param [in] aerosol_modes a list of aerosol modes supported by the Context
  /// @param [in] aerosol_species a list of aerosol species supported by the
  ///                             Context
  /// @param [in] mode_species a map that defines the association of aerosol
  ///             species with an aerosol mode. The keys in this map are names
  ///             of aerosol modes (corresponding to those found in `modes`),
  ///             and the values are lists of symbolic names of aerosol species
  ///             (supplied in `aerosol_species`) that belong to those modes.
  /// @param [in] gas_species a list of gas species supported by the Context
  /// @param [in] gas_chemistry a ChemicalMechanism representing gas chemistry
  /// @param [in] aqueous_chemistry a ChemicalMechanism representing aqueous
  ///                               chemistry
  /// @param [in] num_columns The number of columns in the Context's
  ///                         computational domain
  /// @param [in] num_levels The number of vertical levels in each column within
  ///                        the Context's computational domain
  AeroModel(const Parameterizations& parameterizations,
            const std::vector<Mode>& aerosol_modes,
            const std::vector<Species>& aerosol_species,
            const std::map<std::string, std::vector<std::string> >& mode_species,
            const std::vector<Species>& gas_species,
            ChemicalMechanism* gas_chemistry,
            ChemicalMechanism* aqueous_chemistry,
            int num_columns,
            int num_levels);

  /// AeroModels are not deep-copyable. They should be passed by reference.
  AeroModel(const AeroModel&) = delete;

  /// Destructor.
  ~AeroModel();

  /// AeroModels are not assignable either.
  AeroModel& operator=(const AeroModel&) = delete;

  /// Creates a new AeroState object that can be used with this AeroModel.
  /// All fields within this new AeroState are owned and managed by it. In
  /// the general case, this might not be what you want--in particular, a host
  /// model may demand to manage all of its own state information. Nevertheless,
  /// this simplified state creation can be useful for testing.
  AeroState* create_state() const;

  // Parameterizations

  /// Runs the (prognostic) aerosol process of the given type, computing any
  /// relevant tendencies.
  /// @param [in] type The type of aerosol state to be updated.
  /// @param [in] t The time at which the process runs.
  /// @param [in] dt The time interval over which the process runs.
  /// @param [in] state The aerosol state to be updated.
  /// @param [out] tendencies The aerosol tendencies computed.
  void run_process(AeroProcessType type,
                   Real t, Real dt,
                   const AeroState& state,
                   AeroTendencies& tendencies);

  /// Updates the state with the (diagnostic) aerosol process of the given type.
  /// @param [in] type The type of aerosol state to be updated.
  /// @param [in] t The time at which the process runs.
  /// @param [out] state The aerosol state to be updated.
  void update_state(AeroProcessType type, Real t, AeroState& state);

  // Accessors

  /// Returns the parameterizations associated with this aerosol model.
  const Parameterizations& parameterizations() const;

  /// Returns the list of aerosol modes associated with this aerosol model.
  const std::vector<Mode>& modes() const;

  /// Returns the list of aerosol species associated with this aerosol model.
  const std::vector<Species>& aerosol_species() const;

  /// Returns the list of gas species associated with this aerosol model.
  const std::vector<Species>& gas_species() const;

  /// Returns the chemical mechanism for gas chemistry associated with this
  /// aerosol model.
  const ChemicalMechanism& gas_chemistry() const;

  /// Returns the chemical mechanism for aqueous chemistry associated with this
  /// aerosol model.
  const ChemicalMechanism& aqueous_chemistry() const;

  /// Returns the number of columns in the model.
  int num_columns() const { return num_columns_; }

  /// Returns the number of vertical levels in the model.
  int num_levels() const { return num_levels_; }

  private:

  // Parameterizations, modes, species...
  Parameterizations parameterizations_;
  std::vector<Mode> modes_;
  std::vector<Species> aero_species_;
  std::vector<Species> gas_species_;

  // The association of aerosol species with modes.
  // species_for_modes_[mode_name] = vector of species names
  std::vector<std::vector<int> > species_for_modes_;

  // Chemical mechanisms.
  ChemicalMechanism* gas_chem_;
  ChemicalMechanism* aqueous_chem_;

  // Grid parameters.
  int num_columns_;
  int num_levels_;

  // Selected implementations of prognostic processes used by this model.
  std::map<AeroProcessType, PrognosticAeroProcess*> prog_processes_;

  // Selected implementations of diagnostic processes used by this model.
  std::map<AeroProcessType, DiagnosticAeroProcess*> diag_processes_;
};

}

#endif
