#include "haero/diagnostics.hpp"

namespace haero {

Diagnostics::Diagnostics(int num_aerosol_modes,
                         const std::vector<int>& num_aerosol_species,
                         int num_gases,
                         int num_levels):
  num_aero_species_(num_aerosol_species), num_gases_(num_gases),
  num_levels_(num_levels) {
  EKAT_ASSERT_MSG(num_aerosol_modes == num_aerosol_species.size(),
                  "num_aerosol_species must be a vector of length " << num_aerosol_modes);
}

Diagnostics::~Diagnostics() {
}

int Diagnostics::num_aerosol_modes() const {
  return num_aero_species_.size();
}

int Diagnostics::num_aerosol_species(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < num_aero_species_.size());
  return num_aero_species_[mode_index];
}

int Diagnostics::num_gases() const {
  return num_gases_;
}

int Diagnostics::num_levels() const {
  return num_levels_;
}

bool Diagnostics::has_var(const std::string& name) const {
  return (vars_.find(name) != vars_.end());
}

void Diagnostics::create_var(const std::string& name) {
  auto iter = vars_.find(name);
  EKAT_REQUIRE_MSG(iter == vars_.end(),
    "Diagnostic variable " << name << " already exists!");
  vars_[name] = ColumnView(name, num_levels_);
}

Diagnostics::ColumnView&
Diagnostics::var(const std::string& name) {
  auto iter = vars_.find(name);
  EKAT_REQUIRE_MSG(iter != vars_.end(),
    "Diagnostic variable " << name << " not found!");
  return iter->second;
}

const Diagnostics::ColumnView&
Diagnostics::var(const std::string& name) const {
  auto iter = vars_.find(name);
  EKAT_REQUIRE_MSG(iter != vars_.end(),
    "Diagnostic variable " << name << " not found!");
  return iter->second;
}

bool Diagnostics::has_aerosol_var(const std::string& name) const {
  return (aero_vars_.find(name) != aero_vars_.end());
}

void Diagnostics::create_aerosol_var(const std::string& name) {
  auto iter = aero_vars_.find(name);
  EKAT_REQUIRE_MSG(iter == aero_vars_.end(),
    "Aerosol diagnostic variable " << name << " already exists!");

  // Figure out the dimensions of the view.
  int num_mode_species_pairs = 0;
  for (int m = 0; m < num_aero_species_.size(); ++m) {
    num_mode_species_pairs += num_aero_species_[m];
  }
  int num_vert_packs = num_levels_/HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels_) {
    num_vert_packs++;
  }
  aero_vars_[name] = SpeciesColumnView(name, num_mode_species_pairs,
                                       num_vert_packs);
}

Diagnostics::SpeciesColumnView&
Diagnostics::aerosol_var(const std::string& name) {
  auto iter = aero_vars_.find(name);
  EKAT_REQUIRE_MSG(iter != aero_vars_.end(),
    "Aerosol diagnostic variable " << name << " not found!");
  return iter->second;
}

const Diagnostics::SpeciesColumnView&
Diagnostics::aerosol_var(const std::string& name) const {
  auto iter = aero_vars_.find(name);
  EKAT_REQUIRE_MSG(iter != aero_vars_.end(),
    "Aerosol diagnostic variable " << name << " not found!");
  return iter->second;
}

bool Diagnostics::has_gas_var(const std::string& name) const {
  return (gas_vars_.find(name) != gas_vars_.end());
}

void Diagnostics::create_gas_var(const std::string& name) {
  auto iter = gas_vars_.find(name);
  EKAT_REQUIRE_MSG(iter == gas_vars_.end(),
    "Gas diagnostic variable " << name << " already exists!");
  int num_vert_packs = num_levels_/HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels_) {
    num_vert_packs++;
  }
  gas_vars_[name] = SpeciesColumnView(name, num_gases_, num_vert_packs);
}

Diagnostics::SpeciesColumnView&
Diagnostics::gas_var(const std::string& name) {
  auto iter = gas_vars_.find(name);
  EKAT_REQUIRE_MSG(iter != gas_vars_.end(),
    "Gas diagnostic variable " << name << " not found!");
  return iter->second;
}

const Diagnostics::SpeciesColumnView&
Diagnostics::gas_var(const std::string& name) const {
  auto iter = gas_vars_.find(name);
  EKAT_REQUIRE_MSG(iter != gas_vars_.end(),
    "Gas diagnostic variable " << name << " not found!");
  return iter->second;
}

bool Diagnostics::has_modal_var(const std::string& name) const {
  return (modal_vars_.find(name) != modal_vars_.end());
}

void Diagnostics::create_modal_var(const std::string& name) {
  auto iter = modal_vars_.find(name);
  EKAT_REQUIRE_MSG(iter == modal_vars_.end(),
    "Modal diagnostic variable " << name << " already exists!");
  int num_vert_packs = num_levels_/HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels_) {
    num_vert_packs++;
  }
  modal_vars_[name] = ModalColumnView(name, num_aero_species_.size(),
                                      num_vert_packs);
}

Diagnostics::ModalColumnView&
Diagnostics::modal_var(const std::string& name) {
  auto iter = modal_vars_.find(name);
  EKAT_REQUIRE_MSG(iter != modal_vars_.end(),
    "Modal diagnostic variable " << name << " not found!");
  return iter->second;
}

const Diagnostics::ModalColumnView&
Diagnostics::modal_var(const std::string& name) const {
  auto iter = modal_vars_.find(name);
  EKAT_REQUIRE_MSG(iter != modal_vars_.end(),
    "Modal diagnostic variable " << name << " not found!");
  return iter->second;
}

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
extern "C" {

bool d_has_var_c(void* d, const char* name)
{
  auto* diags = static_cast<Diagnostics*>(d);
  return diags->has_var(name);
}

void* d_var_c(void* d, const char* name)
{
  auto* diags = static_cast<Diagnostics*>(d);
  auto& var = diags->var(name);
  return (void*)var.data();
}

bool d_has_aerosol_var_c(void* d, const char* name)
{
  auto* diags = static_cast<Diagnostics*>(d);
  return diags->has_aerosol_var(name);
}

void* d_aerosol_var_c(void* d, const char* name)
{
  auto* diags = static_cast<Diagnostics*>(d);
  auto& var = diags->aerosol_var(name);
  return (void*)var.data();
}

bool d_has_gas_var_c(void* d, const char* name)
{
  auto* diags = static_cast<Diagnostics*>(d);
  return diags->has_gas_var(name);
}

void* d_gas_var_c(void* d, const char* name)
{
  auto* diags = static_cast<Diagnostics*>(d);
  auto& var = diags->gas_var(name);
  return (void*)var.data();
}

bool d_has_modal_var_c(void* d, const char* name)
{
  auto* diags = static_cast<Diagnostics*>(d);
  return diags->has_modal_var(name);
}

void* d_modal_var_c(void* d, const char* name)
{
  auto* diags = static_cast<Diagnostics*>(d);
  auto& var = diags->modal_var(name);
  return (void*)var.data();
}

} // extern "C"

} // haero

