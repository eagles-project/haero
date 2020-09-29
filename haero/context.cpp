#include "haero/context.hpp"

namespace haero {

Context::Context(
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
  gas_chem_(gas_chemistry),
  aqueous_chem_(aqueous_chemistry),
  num_columns_(num_columns),
  num_levels_(num_levels)
{

}

Context::~Context() {
}

const Parameterizations& Context::parameterizations() const {
  return parameterizations_;
}

const std::vector<Mode>& Context::modes() const {
  return modes_;
}

const std::vector<Species>& Context::aerosol_species() const {
  return aero_species_;
}

const std::vector<Species>& Context::gas_species() const {
  return gas_species_;
}

const ChemicalMechanism& Context::gas_chemistry() const {
  return gas_chem_;
}

const ChemicalMechanism& Context::aqueous_chemistry() const {
  return aqueous_chem_;
}

State Context::create_state() const {
  std::vector<int> num_modal_species; // TODO
  return State(modes_.size(), num_modal_species, gas_species_.size(),
               num_columns_, num_levels_);
}

}
