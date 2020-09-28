#include "haero/context.hpp"

namespace haero {

Context::Context(
  const Parametrizations& parametrizations,
  const std::vector<Mode>& aerosol_modes,
  const std::vector<Species>& aerosol_species,
  const std::map<std::string, std::vector<std::string> >& mode_species,
  const std::vector<Species>& gas_species,
  const ChemicalMechanism& gas_chemistry):
  parametrizations_(parametrizations),
  modes_(aerosol_modes),
  aero_species_(aerosol_species),
  gas_species_(gas_species),
  gas_chem_(gas_chemistry)
{

}

Context::~Context() {
}

const Parametrizations& Context::parametrizations() const {
  return parametrizations_;
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

State Context::create_state() const {
  return State();
}

}
