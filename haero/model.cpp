#include "haero/model.hpp"
#include "ekat/util/ekat_units.hpp"

#include <set>

// Functions for intializing the haero Fortran module.
extern "C" {

using haero::Real;

void haerotran_begin_init();
void haerotran_set_num_modes(int);
void haerotran_set_max_mode_species(int);
void haerotran_set_mode(int, const char*, Real, Real, Real);
void haerotran_set_aero_species(int, int, const char*, const char*,
                                Real, Real, Real);
void haerotran_set_num_gas_species(int);
void haerotran_set_gas_species(int, const char*, const char*,
                               Real, Real, Real);
void haerotran_set_num_columns(int);
void haerotran_set_num_levels(int);
void haerotran_end_init();
void haerotran_finalize();

}

namespace haero {

namespace {

// Indicates whether we have initialized the Haero Fortran helper module. Only
// one model gets to do this, so if more than one model instance has
// Fortran-backed processes, we encounter a fatal error.
bool initialized_fortran_ = false;

// Here are C strings that have been constructed from Fortran strings.
std::set<std::string>* fortran_strings_ = nullptr;

// Prognostic process types.
const ProcessType progProcessTypes[] = {
  ActivationProcess,
  CloudBorneWetRemovalProcess,
  CoagulationProcess,
  CondensationProcess,
  DryDepositionProcess,
  EmissionsProcess,
  NucleationProcess,
  ResuspensionProcess
};

// Diagnostic process types.
const ProcessType diagProcessTypes[] = {
  WaterUptakeProcess
};

} // anonymous namespace

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
  species_for_mode_(),
  num_columns_(num_columns),
  num_levels_(num_levels),
  prog_processes_(),
  diag_processes_(),
  uses_fortran_(false)
{
  // Validate model parameters.
  validate();

  // Set up mode/species indexing.
  index_modal_species(mode_species);

  // Gather processes, and determine whether any of them are backed by Fortran.
  bool have_fortran_processes = gather_processes();

  // If we encountered a Fortran-backed process, attempt to initialize the
  // Haero Fortran helper module with our model data.
  if (have_fortran_processes) {
    init_fortran();
  }

  // Now we can initialize the processes.
  for (auto p: prog_processes_) {
    p.second->init(*this);
  }
  for (auto p: diag_processes_) {
    p.second->init(*this);
  }
}

Model::Model() {
}

Model* Model::ForUnitTests(
  const std::vector<Mode>& aerosol_modes,
  const std::vector<Species>& aerosol_species,
  const std::map<std::string, std::vector<std::string> >& mode_species,
  const std::vector<Species>& gas_species,
  int num_columns,
  int num_levels) {
  Model* model = new Model();
  model->modes_ = aerosol_modes;
  model->aero_species_ = aerosol_species;
  model->gas_species_ = gas_species;
  model->num_columns_ = num_columns;
  model->num_levels_ = num_levels;
  model->uses_fortran_ = false;

  // Validate model parameters.
  model->validate();

  // Set up mode/species indexing.
  model->index_modal_species(mode_species);

  // Initialize Fortran
  model->init_fortran();

  return model;
}

Model::~Model() {
  for (auto p: prog_processes_) {
    delete p.second;
  }
  for (auto p: diag_processes_) {
    delete p.second;
  }

  // If we initialized the Haero Fortran module, we must finalize it.
  if (uses_fortran_) {
    haerotran_finalize();
    initialized_fortran_ = false;

    // Delete all relevant Fortran string instances.
    if (fortran_strings_ != nullptr) {
      delete fortran_strings_;
      fortran_strings_ = nullptr;
    }
  }
}

Prognostics* Model::create_prognostics() const {

  auto progs = new Prognostics(num_columns_, num_levels_);

  // Add aerosol modes/species data.
  for (size_t i = 0; i < modes_.size(); ++i) {
    std::vector<Species> species;
    for (size_t j = 0; j < species_for_mode_[i].size(); ++j) {
      species.push_back(aero_species_[species_for_mode_[i][j]]);
    }
    progs->add_aerosol_mode(modes_[i], species);
  }

  // Add gas species data.
  progs->add_gas_species(gas_species_);

  // Assemble the prognostics and return them.
  progs->assemble();
  return progs;
}

Diagnostics* Model::create_diagnostics() const {
  // Create an empty Diagnostics object.
  std::vector<int> num_aero_species(modes_.size());
  for (size_t m = 0; m < modes_.size(); ++m) {
    num_aero_species[m] = static_cast<int>(species_for_mode_[m].size());
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

void Model::update_diagnostics(ProcessType type,
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

std::vector<Species> Model::aerosol_species_for_mode(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < species_for_mode_.size());
  // Construct this vector from our association data.
  std::vector<Species> species;
  for (int s = 0; s < species_for_mode_[mode_index].size(); ++s) {
    species.push_back(aero_species_[species_for_mode_[mode_index][s]]);
  }
  return species;
}

const std::vector<Species>& Model::gas_species() const {
  return gas_species_;
}

void Model::index_modal_species(const std::map<std::string, std::vector<std::string> >& mode_species) {
  species_for_mode_.resize(modes_.size());
  for (auto iter = mode_species.begin(); iter != mode_species.end(); ++iter) {
    const auto& mode_name = iter->first;
    const auto& mode_species = iter->second;

    auto m_iter = std::find_if(modes_.begin(), modes_.end(),
      [&] (const Mode& mode) { return mode.name == mode_name; });
    int mode_index = m_iter - modes_.begin();

    for (int s = 0; s < mode_species.size(); ++s) {
      auto s_iter = std::find_if(aero_species_.begin(), aero_species_.end(),
        [&] (const Species& species) { return species.symbol == mode_species[s]; });
      int species_index = s_iter - aero_species_.begin();
      species_for_mode_[mode_index].push_back(species_index);
    }
  }

  // Make sure each mode contains at least one species.
  for (int m = 0; m < species_for_mode_.size(); ++m) {
    EKAT_REQUIRE_MSG(not species_for_mode_[m].empty(),
      modes_[m].name.c_str() << " mode contains no aerosol species!");
  }
}

void Model::init_fortran() {
  // If we've already initialized the Fortran module, this means that this is
  // not the first C++ model instance that has Fortran-backed processes. We
  // đon't allow this, since we've made assumptions in order to simplify
  // the process of implementing Fortran processes.
  EKAT_REQUIRE_MSG(not initialized_fortran_,
      "More than one C++ model includes Fortran-backed processes. This is not allowed!");

  // This series of calls sets things up in Haero's Fortran module.
  haerotran_begin_init();
  int num_modes = modes_.size();
  haerotran_set_num_modes(num_modes);
  int max_species = 0;
  for (int i = 0; i < num_modes; ++i) {
    max_species = std::max(max_species, (int)species_for_mode_[i].size());
  }
  haerotran_set_max_mode_species(max_species);
  for (int i = 0; i < num_modes; ++i) {
    // Set the properties of mode i+1 (as indexed in Fortran).
    const auto& mode = modes_[i];
    haerotran_set_mode(i+1, mode.name.c_str(), mode.min_diameter,
        mode.max_diameter, mode.mean_std_dev);

    // Set up aerosol species for this mode.
    int num_species = species_for_mode_[i].size();
    for (int j = 0; j < num_species; ++j) {
      const auto& species = aero_species_[species_for_mode_[i][j]];
      haerotran_set_aero_species(i+1, j+1, species.name.c_str(),
          species.symbol.c_str(), species.molecular_weight,
          species.crystalization_point, species.deliquescence_point);
    }
  }

  // Set up gas species.
  int num_gas_species = gas_species_.size();
  haerotran_set_num_gas_species(num_gas_species);
  for (int i = 0; i < num_gas_species; ++i) {
    const auto& species = gas_species_[i];
    haerotran_set_gas_species(i+1, species.name.c_str(),
        species.symbol.c_str(), species.molecular_weight,
        species.crystalization_point, species.deliquescence_point);
  }

  // Set dimensions.
  haerotran_set_num_columns(num_columns_);
  haerotran_set_num_levels(num_levels_);
  haerotran_end_init();

  // Okay, the Fortran module is initialized.
  initialized_fortran_ = true;
  uses_fortran_ = true;
}

bool Model::gather_processes() {

  // We determine whether we have Fortran-backed processes, and map the
  // aerosol process types to these processes. We don't initialize the processes
  // yet, because that requires that the Fortran representation of the Model
  // needs to be set up before that happens.
  bool have_fortran_processes = false;
  for (auto p: progProcessTypes) {
    PrognosticProcess* process = select_prognostic_process(p, parameterizations_);
    if (dynamic_cast<FPrognosticProcess*>(process) != nullptr) { // Fortran-backed!
      have_fortran_processes = true;
    }
    prog_processes_[p] = process;
  }

  for (auto p: diagProcessTypes) {
    DiagnosticProcess* process = select_diagnostic_process(p, parameterizations_);
    if (dynamic_cast<FPrognosticProcess*>(process) != nullptr) { // Fortran-backed!
      have_fortran_processes = true;
    }
    diag_processes_[p] = process;
  }

  return have_fortran_processes;
}

void Model::validate() {
  EKAT_REQUIRE_MSG(not modes_.empty(), "Model: No modes were defined!");
  EKAT_REQUIRE_MSG(not aero_species_.empty(), "Model: No aerosol species were given!");
  EKAT_REQUIRE_MSG(not gas_species_.empty(), "Model: No gas species were given!");
  EKAT_REQUIRE_MSG((num_columns_ > 0), "Model: No columns were specified!");
  EKAT_REQUIRE_MSG((num_levels_ > 0), "Model: No vertical levels were specified!");
}

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
extern "C" {

// Returns a newly-allocated C string for the given Fortran string pointer with
// the given length. Resources for this string are managed by the running model.
const char* new_c_string(char* f_str_ptr, int f_str_len) {
  if (fortran_strings_ == nullptr) {
    fortran_strings_ = new std::set<std::string>();
  }

  // Before we allocate any more strings, check to see whether we've already
  // got this one.
  std::string f_str(f_str_ptr, f_str_len);
  auto iter = fortran_strings_->find(f_str);
  if (iter != fortran_strings_->end()) { // got it!
    return iter->c_str();
  } else { // we need a new one
    auto inserted = fortran_strings_->insert(f_str);
    return inserted.first->c_str();
  }
}

} // extern "C"

}
