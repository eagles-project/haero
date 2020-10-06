#include "haero/aero_model.hpp"
#include "ekat/util/ekat_units.hpp"

namespace haero {

AeroModel::AeroModel(
  const Parameterizations& parameterizations,
  const std::vector<Mode>& aerosol_modes,
  const std::vector<Species>& aerosol_species,
  const std::map<std::string, std::vector<std::string> >& mode_species,
  const std::vector<Species>& gas_species,
  ChemicalMechanism* gas_chemistry,
  ChemicalMechanism* aqueous_chemistry,
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
  EKAT_ASSERT_MSG(gas_chemistry != nullptr,
                  "A gas chemistry mechanism must be provided!");
  EKAT_ASSERT_MSG(aqueous_chemistry != nullptr,
                  "An aqueous chemistry mechanism must be provided!");

  // Set up mode/species indexing.
  // TODO
}

AeroModel::~AeroModel() {
}

AeroState* AeroModel::create_state() const {

  auto state = new AeroState(num_columns_, num_levels_);

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

void AeroModel::run_water_uptake(const AeroState& state,
                                 AeroTendencies& tendencies) {
}

const Parameterizations& AeroModel::parameterizations() const {
  return parameterizations_;
}

const std::vector<Mode>& AeroModel::modes() const {
  return modes_;
}

const std::vector<Species>& AeroModel::aerosol_species() const {
  return aero_species_;
}

const std::vector<Species>& AeroModel::gas_species() const {
  return gas_species_;
}

const ChemicalMechanism& AeroModel::gas_chemistry() const {
  return *gas_chem_;
}

const ChemicalMechanism& AeroModel::aqueous_chemistry() const {
  return *aqueous_chem_;
}

}
