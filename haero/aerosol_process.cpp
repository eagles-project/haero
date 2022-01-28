#include "haero/aerosol_process.hpp"

#include <set>

#include "haero/utils.hpp"

namespace {

#if HAERO_FORTRAN

// Indicates whether we have initialized the Haero Fortran helper module. Only
// one model gets to do this, so if more than one model instance has
// Fortran-backed processes, we encounter a fatal error.
bool initialized_fortran_ = false;

// Stored aerosol configuration for Fortran.
haero::ModalAerosolConfig* config_ = nullptr;

// Here are C strings that have been constructed from Fortran strings.
std::set<std::string>* fortran_strings_ = nullptr;

// Functions for intializing the haero Fortran module.
extern "C" {

using haero::Real;

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
void haerotran_begin_init();
void haerotran_set_num_aerosol_modes(int);
void haerotran_set_max_mode_species(int);
void haerotran_set_aerosol_mode(int, const char*, Real, Real, Real, Real, Real,
                                Real);
void haerotran_set_aerosol_species(int, int, const char*, const char*, Real,
                                   Real, Real);
void haerotran_set_num_gas_species(int);
void haerotran_set_gas_species(int, const char*, const char*, Real);
void haerotran_end_init();
void haerotran_finalize();

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

// This function frees all resources associated with the Fortran subsystem.
void finalize_fortran() {
  // If we initialized the Haero Fortran module, we must finalize it.
  if (initialized_fortran_) {
    haerotran_finalize();
    delete config_;
    config_ = nullptr;
    initialized_fortran_ = false;

    // Delete all relevant Fortran string instances.
    if (fortran_strings_ != nullptr) {
      delete fortran_strings_;
      fortran_strings_ = nullptr;
    }
  }
}

}  // extern "C"

#endif  // HAERO_FORTRAN

}  // anonymous namespace

namespace haero {

#if HAERO_FORTRAN
void AerosolProcess::init_fortran_(const ModalAerosolConfig& config) {
  if (initialized_fortran_) {
    // Check to make sure our aerosol configurations are compatible!
    EKAT_REQUIRE_MSG(
        config == *config_,
        "Two incompatible Fortran aerosol processes have been created!\n"
        "All Fortran aerosol processes must use the same configuration.");
    return;
  }

  // This series of calls sets things up in Haero's Fortran module.
  config_ = new ModalAerosolConfig(config);
  haerotran_begin_init();
  int num_modes = config_->aerosol_modes.size();
  haerotran_set_num_aerosol_modes(num_modes);
  size_t max_species = 0;
  for (int m = 0; m < num_modes; ++m) {
    const auto mode_species = config_->aerosol_species_for_mode(m);
    max_species = std::max(max_species, mode_species.size());
  }
  haerotran_set_max_mode_species(max_species);
  for (int m = 0; m < num_modes; ++m) {
    const auto mode_species = config_->aerosol_species_for_mode(m);

    // Set the properties of mode i+1 (as indexed in Fortran).
    const auto& mode = config_->aerosol_modes[m];
    haerotran_set_aerosol_mode(m + 1, mode.name().c_str(), mode.min_diameter,
                               mode.nom_diameter, mode.max_diameter,
                               mode.mean_std_dev, mode.deliquescence_pt,
                               mode.crystallization_pt);

    // Set up aerosol species for this mode.
    int num_aero_species = mode_species.size();
    for (int s = 0; s < num_aero_species; ++s) {
      const auto species = mode_species[s];
      haerotran_set_aerosol_species(
          m + 1, s + 1, species.name().c_str(), species.symbol().c_str(),
          species.molecular_weight, species.density, species.hygroscopicity);
    }
  }

  // Set up gas species.
  int num_gas_species = config_->gas_species.size();
  haerotran_set_num_gas_species(num_gas_species);
  for (int i = 0; i < num_gas_species; ++i) {
    const auto species = config_->gas_species[i];
    haerotran_set_gas_species(i + 1, species.name().c_str(),
                              species.symbol().c_str(),
                              species.molecular_weight);
  }

  // Set dimensions.
  haerotran_end_init();

  // Okay, the Fortran module is initialized.
  initialized_fortran_ = true;
  // Register our Fortran finalization function to be called when the process
  // exits.
  atexit(finalize_fortran);
}

#endif

void AerosolProcess::interpret_and_set_param(const std::string& name,
                                             const std::string& value) {
  // Is this an integer?
  try {
    int int_value = std::stoi(value);
    set_param(name, int_value);
  } catch (const std::invalid_argument&) {
    // Okay, it's not an integer. Is it a real number?
    try {
#if HAERO_DOUBLE_PRECISION
      Real real_value = std::stod(value);
#else
      Real real_value = std::stof(value);
#endif
      set_param(name, real_value);
    } catch (const std::invalid_argument&) {
      // Boolean?
      if (is_boolean(value)) {
        set_param(name, as_boolean(value));
      } else {
        // Okay, we can only interpret this value as a string. String parameters
        // aren't supported for Fortran processes, so make sure we're not one of
        // those.
#if HAERO_FORTRAN
        if (dynamic_cast<haero::FAerosolProcess*>(this) != nullptr) {
          fprintf(
              stderr,
              "Parameter '%s' with string value '%s' given for\n"
              "Fortran aerosol process %s. Fortran aerosol processes cannot\n"
              "accept string parameter values.",
              name.c_str(), value.c_str(), name_.label().c_str());
          return;
        }
#endif
        set_param(name, value);
      }
    }
  }
}

}  // namespace haero
