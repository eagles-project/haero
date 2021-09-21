#include "haero/model.hpp"

#include <set>

#include "ekat/util/ekat_units.hpp"

#if HAERO_FORTRAN
// Functions for intializing the haero Fortran module.
extern "C" {

using haero::Real;

void haerotran_begin_init();
void haerotran_set_num_modes(int);
void haerotran_set_max_mode_species(int);
void haerotran_set_mode(int, const char*, Real, Real, Real, Real, Real);
void haerotran_set_aero_species(int, int, const char*, const char*, Real, Real,
                                Real);
void haerotran_set_num_gas_species(int);
void haerotran_set_gas_species(int, const char*, const char*, Real);
void haerotran_set_num_levels(int);
void haerotran_end_init();
void haerotran_finalize();
}
#endif  // HAERO_FORTRAN

namespace haero {

namespace {

#if HAERO_FORTRAN
// Indicates whether we have initialized the Haero Fortran helper module. Only
// one model gets to do this, so if more than one model instance has
// Fortran-backed processes, we encounter a fatal error.
bool initialized_fortran_ = false;

// Here are C strings that have been constructed from Fortran strings.
std::set<std::string>* fortran_strings_ = nullptr;

#endif  // HAERO_FORTRAN

}  // anonymous namespace

Model::Model(const ModalAerosolConfig& modal_aerosol_config,
             const SelectedProcesses& selected_processes, int num_levels)
    : config_(modal_aerosol_config),
      selected_processes_(selected_processes),
      num_levels_(num_levels),
      aero_processes_(),
      uses_fortran_(false) {
  // Validate model parameters.
  validate();

  // Gather processes, and determine whether any of them are backed by Fortran.
  bool have_fortran_processes = gather_processes();

  // If we encountered a Fortran-backed process, attempt to initialize the
  // Haero Fortran helper module with our model data.
  if (have_fortran_processes) {
#if HAERO_FORTRAN
    init_fortran();
#endif  // HAERO_FORTRAN
  }

  // Now we can initialize the processes.
  for (auto p : aero_processes_) {
    p.second->init(config_);
  }
}

Model::Model() {}

Model::Model(const Model& model)
    : config_(model.config_),
      selected_processes_(model.selected_processes_),
      num_levels_(model.num_levels_),
      aero_processes_(model.aero_processes_),
      uses_fortran_(model.uses_fortran_) {}

Model* Model::ForUnitTests(const ModalAerosolConfig& modal_aerosol_config,
                           int num_levels) {
  Model* model = new Model();
  model->config_ = modal_aerosol_config;
  model->num_levels_ = num_levels;
  model->uses_fortran_ = false;

  // Validate model parameters.
  model->validate();

  // Initialize Fortran
  model->init_fortran();

  return model;
}

Model::~Model() {
  for (auto p : aero_processes_) {
    delete p.second;
  }

#if HAERO_FORTRAN
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
#endif  // HAERO_FORTRAN
}

Prognostics* Model::create_prognostics(
    SpeciesColumnView int_aerosols, SpeciesColumnView cld_aerosols,
    SpeciesColumnView gases, ModeColumnView interstitial_num_mix_ratios,
    ModeColumnView cloudborne_num_mix_ratios) const {
  std::vector<int> num_aero_species(config_.aerosol_modes.size());
  for (size_t m = 0; m < config_.aerosol_modes.size(); ++m) {
    const auto mode_species = config_.aerosol_species_for_mode(m);
    num_aero_species[m] = static_cast<int>(mode_species.size());
  }
  return new Prognostics(num_aero_species.size(), num_aero_species,
                         config_.gas_species.size(), num_levels_, int_aerosols,
                         cld_aerosols, gases, interstitial_num_mix_ratios,
                         cloudborne_num_mix_ratios);
}

HostDiagnostics* Model::create_diagnostics() const {
  // Create an empty Diagnostics object.
  std::vector<int> num_aero_species(config_.num_modes());
  for (size_t m = 0; m < num_aero_species.size(); ++m) {
    const auto mode_species = config_.aerosol_species_for_mode(m);
    num_aero_species[m] = static_cast<int>(mode_species.size());
  }
  auto diags = new HostDiagnostics(num_aero_species.size(), num_aero_species,
                                   config_.num_gases(), num_levels_);

  //  // Make sure that all diagnostic variables needed by the model's processes
  //  // are present.
  //  for (auto iter = diag_processes_.begin(); iter != diag_processes_.end();
  //  ++iter) {
  //    iter->second->prepare(*diags);
  //  }

  return diags;
}

AerosolProcess* Model::process_on_device(AerosolProcessType type) const {
  auto iter = aero_processes_.find(type);
  EKAT_REQUIRE_MSG(iter != aero_processes_.end(),
                   "No process of the selected type is available!");
  EKAT_REQUIRE_MSG(iter->second != nullptr,
                   "Null process pointer encountered!");
  EKAT_REQUIRE_MSG(iter->second->type() == type,
                   "Invalid process type encountered!");
  return iter->second->copy_to_device();
}

const SelectedProcesses& Model::selected_processes() const {
  return selected_processes_;
}

void Model::init_fortran() {
#if HAERO_FORTRAN
  // If we've already initialized the Fortran module, this means that this is
  // not the first C++ model instance that has Fortran-backed processes. We
  // đon't allow this, since we've made assumptions in order to simplify
  // the process of implementing Fortran processes.
  EKAT_REQUIRE_MSG(not initialized_fortran_,
                   "More than one C++ model includes Fortran-backed processes. "
                   "This is not allowed!");

  // This series of calls sets things up in Haero's Fortran module.
  haerotran_begin_init();
  int num_modes = config_.aerosol_modes.size();
  haerotran_set_num_modes(num_modes);
  size_t max_species = 0;
  for (int m = 0; m < num_modes; ++m) {
    const auto mode_species = config_.aerosol_species_for_mode(m);
    max_species = std::max(max_species, mode_species.size());
  }
  haerotran_set_max_mode_species(max_species);
  for (int m = 0; m < num_modes; ++m) {
    const auto mode_species = config_.aerosol_species_for_mode(m);

    // Set the properties of mode i+1 (as indexed in Fortran).
    const auto& mode = config_.aerosol_modes[m];
    haerotran_set_mode(m + 1, mode.name().c_str(), mode.min_diameter,
                       mode.max_diameter, mode.mean_std_dev,
                       mode.deliquescence_pt, mode.crystallization_pt);

    // Set up aerosol species for this mode.
    int num_species = mode_species.size();
    for (int s = 0; s < num_species; ++s) {
      const auto species = config_.aerosol_species[s];
      haerotran_set_aero_species(
          m + 1, s + 1, species.name().c_str(), species.symbol().c_str(),
          species.molecular_weight, species.density, species.hygroscopicity);
    }
  }

  // Set up gas species.
  int num_gas_species = config_.gas_species.size();
  haerotran_set_num_gas_species(num_gas_species);
  for (int i = 0; i < num_gas_species; ++i) {
    const auto species = config_.gas_species[i];
    haerotran_set_gas_species(i + 1, species.name().c_str(),
                              species.symbol().c_str(),
                              species.molecular_weight);
  }

  // Set dimensions.
  haerotran_set_num_levels(num_levels_);
  haerotran_end_init();

  // Okay, the Fortran module is initialized.
  initialized_fortran_ = true;
  uses_fortran_ = true;
#endif  // HAERO_FORTRAN
}

bool Model::gather_processes() {
  // We determine whether we have Fortran-backed processes, and map the
  // aerosol process types to these processes. We don't initialize the processes
  // yet, because that requires that the Fortran representation of the Model
  // needs to be set up before that happens.
  bool have_fortran_processes = false;
  for (int p = ActivationProcess; p != Terminator; ++p) {
    auto pType = static_cast<AerosolProcessType>(p);
    AerosolProcess* process =
        select_aerosol_process(pType, selected_processes_);
#if HAERO_FORTRAN
    if (dynamic_cast<FAerosolProcess*>(process) !=
        nullptr) {  // Fortran-backed!
      have_fortran_processes = true;
    }
#endif  // HAERO_FORTRAN
    aero_processes_[pType] = process;
  }

  return have_fortran_processes;
}

void Model::validate() {
  EKAT_REQUIRE_MSG(config_.aerosol_modes.size(),
                   "Model: No modes were defined!");
  EKAT_REQUIRE_MSG(config_.aerosol_species.size(),
                   "Model: No aerosol species were given!");
  EKAT_REQUIRE_MSG(config_.gas_species.size(),
                   "Model: No gas species were given!");
  EKAT_REQUIRE_MSG((num_levels_ > 0),
                   "Model: No vertical levels were specified!");
}

#if HAERO_FORTRAN
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
  if (iter != fortran_strings_->end()) {  // got it!
    return iter->c_str();
  } else {  // we need a new one
    auto inserted = fortran_strings_->insert(f_str);
    return inserted.first->c_str();
  }
}

}  // extern "C"
#endif  // HAERO_FORTRAN

}  // namespace haero
