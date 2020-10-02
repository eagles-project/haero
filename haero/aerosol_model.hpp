#ifndef HAERO_AEROSOL_MODEL_HPP
#define HAERO_AEROSOL_MODEL_HPP

#include "haero/mode.hpp"
#include "haero/species.hpp"
#include "haero/chemistry.hpp"
#include "haero/parameterizations.hpp"
#include "haero/aerosol_state.hpp"
#include "haero/aerosol_tendencies.hpp"
#include <map>

namespace haero {

/// This type represents an aerosol system to be simulated, including all
/// information about modes, species, chemistry, and selected parametrizations.
class AerosolModel final {
  public:

  /// Creates a AerosolModel that supports the specified parameterizations.
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
  AerosolModel(const Parameterizations& parameterizations,
               const std::vector<Mode>& aerosol_modes,
               const std::vector<Species>& aerosol_species,
               const std::map<std::string, std::vector<std::string> >& mode_species,
               const std::vector<Species>& gas_species,
               const ChemicalMechanism& gas_chemistry,
               const ChemicalMechanism& aqueous_chemistry,
               int num_columns,
               int num_levels);

  /// AerosolModels are not deep-copyable. They should be passed by reference.
  AerosolModel(const AerosolModel&) = delete;

  /// Destructor.
  ~AerosolModel();

  /// AerosolModels are not assignable either.
  AerosolModel& operator=(const AerosolModel&) = delete;

  /// Creates a new AerosolState object that can be used with this AerosolModel.
  /// All fields within this new AerosolState are owned and managed by it. In
  /// the general case, this might not be what you want--in particular, a host
  /// model may demand to manage all of its own state information. Nevertheless,
  /// this simplified state creation can be useful for testing.
  AerosolState* create_state() const;

  // Parameterizations

  /// Runs the selected water uptake parametrization, computing related
  /// tendencies (time derivatives).
  /// @param [in] state The aerosol state to be updated.
  /// @param [out] tendencies The aerosol tendencies computed.
  void run_water_uptake(const AerosolState& state,
                        AerosolTendencies& tendencies);

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
  ChemicalMechanism gas_chem_;
  ChemicalMechanism aqueous_chem_;

  // Grid parameters.
  int num_columns_;
  int num_levels_;
};

}

#endif
