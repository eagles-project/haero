#include "haero/model.hpp"
#include "ekat/util/ekat_units.hpp"

namespace haero {

Model::Model(
  const Parameterizations& parameterizations,
  const std::vector<Mode>& aerosol_modes,
  const std::vector<Species>& aerosol_species,
  const std::map<std::string, std::vector<std::string> >& mode_species,
  const std::vector<Species>& gas_species,
  int num_columns,
  int num_levels):
  parameterizations_(parameterizations),
  modes_(aerosol_modes),
  aero_species_(aerosol_species),
  gas_species_(gas_species),
  species_for_modes_(),
  num_columns_(num_columns),
  num_levels_(num_levels),
  prog_processes_(),
  diag_processes_()
{
  // Set up mode/species indexing.
  // TODO

  // Set up prognostic and diagnostic processes.
  ProcessType progProcessTypes[] = {
    ActivationProcess,
    CloudBorneWetRemovalProcess,
    CoagulationProcess,
    CondensationProcess,
    DryDepositionProcess,
    EmissionsProcess,
    NucleationProcess,
    ResuspensionProcess
  };
  for (auto p: progProcessTypes) {
    prog_processes_[p] = select_prognostic_process(p, parameterizations);
    prog_processes_[p]->init(*this);
  }

  ProcessType diagProcessTypes[] = {
    WaterUptakeProcess
  };
  for (auto p: diagProcessTypes) {
    diag_processes_[p] = select_diagnostic_process(p, parameterizations);
    diag_processes_[p]->init(*this);
  }
}

Model::~Model() {
  for (auto p: prog_processes_) {
    delete p.second;
  }
  for (auto p: diag_processes_) {
    delete p.second;
  }
}

Prognostics* Model::create_prognostics() const {

  auto progs = new Prognostics(num_columns_, num_levels_);

  // Add aerosol modes/species data.
  for (size_t i = 0; i < modes_.size(); ++i) {
    std::vector<Species> species;
    for (size_t j = 0; j < species_for_modes_[i].size(); ++j) {
      species.push_back(aero_species_[species_for_modes_[i][j]]);
    }
    progs->add_aerosol_mode(modes_[i], species);
  }

  // Add gas species data.
  progs->add_gas_species(gas_species_);

  return progs;
}

Diagnostics* Model::create_diagnostics() const {
  // Create an empty Diagnostics object.
  std::vector<int> num_aero_species(modes_.size());
  for (size_t m = 0; m < modes_.size(); ++m) {
    num_aero_species[m] = static_cast<int>(species_for_modes_[m].size());
  }
  auto diags =  new Diagnostics(num_columns_, num_levels_,
                                num_aero_species, gas_species_.size());

  // Make sure that all diagnostic variables needed by the model's processes
  // are present.
  for (auto iter = diag_processes_.begin(); iter != diag_processes_.end(); ++iter) {
    iter->second->prepare(*diags);
  }

  return diags;
}

void Model::run_process(ProcessType type,
                        Real t, Real dt,
                        const Prognostics& prognostics,
                        const Diagnostics& diagnostics,
                        Tendencies& tendencies) {
  auto iter = prog_processes_.find(type);
  EKAT_REQUIRE_MSG(iter != prog_processes_.end(),
                   "No process of the selected type is available!");
  EKAT_REQUIRE_MSG(iter->second != nullptr,
                   "Null process pointer encountered!");
  EKAT_REQUIRE_MSG(iter->second->type() == type,
                   "Invalid process type encountered!");
  iter->second->run(*this, t, dt, prognostics, diagnostics, tendencies);
}

void Model::update_state(ProcessType type,
                         Real t,
                         const Prognostics& prognostics,
                         Diagnostics& diagnostics) {
  auto iter = diag_processes_.find(type);
  EKAT_REQUIRE_MSG(iter != diag_processes_.end(),
                   "No process of the selected type is available!");
  EKAT_REQUIRE_MSG(iter->second != nullptr,
                   "Null process pointer encountered!");
  EKAT_REQUIRE_MSG(iter->second->type() == type,
                   "Invalid process type encountered!");
  iter->second->update(*this, t, prognostics, diagnostics);
}

const Parameterizations& Model::parameterizations() const {
  return parameterizations_;
}

const std::vector<Mode>& Model::modes() const {
  return modes_;
}

const std::vector<Species>& Model::aerosol_species() const {
  return aero_species_;
}

const std::vector<Species>& Model::gas_species() const {
  return gas_species_;
}

}
