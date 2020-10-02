#include "haero/aerosol_model.hpp"
#include "ekat/util/ekat_units.hpp"

namespace haero {

AerosolModel::AerosolModel(
  const Parameterizations& parameterizations,
  const std::vector<Mode>& aerosol_modes,
  const std::vector<Species>& aerosol_species,
  const std::map<std::string, std::vector<std::string> >& mode_species,
  const std::vector<Species>& gas_species,
  const ChemicalMechanism& gas_chemistry,
  const ChemicalMechanism& aqueous_chemistry,
  int num_columns,
  int num_levels):
  parameterizations_(parameterizations),
  modes_(aerosol_modes),
  aero_species_(aerosol_species),
  gas_species_(gas_species),
  species_for_modes_(),
  gas_chem_(gas_chemistry),
  aqueous_chem_(aqueous_chemistry),
  num_columns_(num_columns),
  num_levels_(num_levels)
{
  // Set up mode/species indexing.
  // TODO
}

AerosolModel::~AerosolModel() {
}

AerosolState* AerosolModel::create_state() const {

  auto state = new AerosolState(num_columns_, num_levels_);

  // Add aerosol modes/species data.
  for (size_t i = 0; i < modes_.size(); ++i) {
    std::vector<Species> species;
    for (size_t j = 0; j < species_for_modes_[i].size(); ++j) {
      species.push_back(aero_species_[species_for_modes_[i][j]]);
    }
    state->add_aerosol_mode(modes_[i], species);
  }

  // Add gas species data.
  state->add_gas_species(gas_species_);

  return state;
}

void AerosolModel::run_water_uptake(const AerosolState& state,
                                    AerosolTendencies& tendencies) {
}

const Parameterizations& AerosolModel::parameterizations() const {
  return parameterizations_;
}

const std::vector<Mode>& AerosolModel::modes() const {
  return modes_;
}

const std::vector<Species>& AerosolModel::aerosol_species() const {
  return aero_species_;
}

const std::vector<Species>& AerosolModel::gas_species() const {
  return gas_species_;
}

const ChemicalMechanism& AerosolModel::gas_chemistry() const {
  return gas_chem_;
}

const ChemicalMechanism& AerosolModel::aqueous_chemistry() const {
  return aqueous_chem_;
}

}
