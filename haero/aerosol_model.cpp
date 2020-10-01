#include "haero/context.hpp"

namespace haero {

AerosolModel::Context(
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
  species_for_modes_(mode_species),
  gas_chem_(gas_chemistry),
  aqueous_chem_(aqueous_chemistry),
  num_columns_(num_columns),
  num_levels_(num_levels)
{
}

AerosolModel::~Context() {
}

AerosolState* AerosolModel::create_state() const {
  std::vector<std::string> mode_names(modes_.size());
  for (size_t i = 0; i < modes_.size(); ++i) {
    mode_names[i] = modes_[i].name();
  }
  auto state = new AerosolState(modes);

  // Add fields for modal number densities.
  auto num_density_units = pow(ekat::units::m, -3);
  for (size_t i = 0; i < modes_.size(); ++i) {
    auto field_name = modes_[i].name + std::string("_num_density");
    state.create_field(field_name, num_density_units);
  }

  // Add fields for modal aerosol species mix fractions.
  auto mix_frac_units = ekat::units::Units(1);
  for (size_t i = 0; i < modes_.size(); ++i) {
    auto field_name = gas_species_[i].symbol + std::string("_mix_frac");
    state.create_modal_field(modes_[i].name(), field_name, mix_frac_units);
  }

  // Add fields for gas species mole fractions.
  auto mole_frac_units = ekat::units::Units(1);
  for (size_t i = 0; i < modes_.size(); ++i) {
    auto field_name = gas_species_[i].symbol + std::string("_mole_frac");
    state.create_field(field_name, mole_frac_units);
  }

  return state;
}

void AerosolModel::run_water_uptake(AerosolState& state) {
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
