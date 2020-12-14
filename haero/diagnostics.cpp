#include "haero/diagnostics.hpp"

namespace haero {

Diagnostics::Diagnostics(int num_columns, int num_levels,
                         const std::vector<int>& num_aero_species,
                         int num_gas_species):
  num_columns_(num_columns), num_levels_(num_levels),
  num_aero_species_(num_aero_species), num_gas_species_(num_gas_species) {
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

int Diagnostics::num_gas_species() const {
  return num_gas_species_;
}

int Diagnostics::num_columns() const {
  return num_columns_;
}

int Diagnostics::num_levels() const {
  return num_levels_;
}

bool Diagnostics::has_var(const std::string& name) const {
  return (vars_.find(name) != vars_.end());
}

void Diagnostics::create_var(const std::string& name) {
  auto iter = vars_.find(name);
  EKAT_REQUIRE_MSG(iter == vars_.end(), "Diagnostic variable already exists!");
  vars_[name] = ColumnView(name, num_columns_, num_levels_);
}

Diagnostics::ColumnView&
Diagnostics::var(const std::string& name) {
  auto iter = vars_.find(name);
  EKAT_REQUIRE_MSG(iter != vars_.end(), "Diagnostic variable not found!");
  return iter->second;
}

const Diagnostics::ColumnView&
Diagnostics::var(const std::string& name) const {
  auto iter = vars_.find(name);
  EKAT_REQUIRE_MSG(iter != vars_.end(), "Diagnostic variable not found!");
  return iter->second;
}

bool Diagnostics::has_aerosol_var(const std::string& name) const {
  return (aero_vars_.find(name) != aero_vars_.end());
}

void Diagnostics::create_aerosol_var(const std::string& name) {
  auto iter = aero_vars_.find(name);
  EKAT_REQUIRE_MSG(iter == aero_vars_.end(), "Aerosol diagnostic variable already exists!");
  aero_vars_[name] = std::vector<ColumnSpeciesView>();
  for (int m = 0; m < num_aerosol_modes(); ++m) {
    aero_vars_[name].push_back(ColumnSpeciesView(name, num_columns_,
                                                 num_levels_,
                                                 num_aero_species_[m]));
  }
}

Diagnostics::ColumnSpeciesView&
Diagnostics::aerosol_var(const std::string& name, int mode_index) {
  auto iter = aero_vars_.find(name);
  EKAT_REQUIRE_MSG(iter != aero_vars_.end(), "Aerosol diagnostic variable not found!");
  return iter->second[mode_index];
}

const Diagnostics::ColumnSpeciesView&
Diagnostics::aerosol_var(const std::string& name, int mode_index) const {
  auto iter = aero_vars_.find(name);
  EKAT_REQUIRE_MSG(iter != aero_vars_.end(), "Aerosol diagnostic variable not found!");
  return iter->second[mode_index];
}

bool Diagnostics::has_gas_var(const std::string& name) const {
  return (gas_vars_.find(name) != gas_vars_.end());
}

void Diagnostics::create_gas_var(const std::string& name) {
  auto iter = gas_vars_.find(name);
  EKAT_REQUIRE_MSG(iter == gas_vars_.end(), "Gas diagnostic variable already exists!");
  gas_vars_[name] = ColumnSpeciesView(name, num_columns_, num_levels_, num_gas_species_);
}

Diagnostics::ColumnSpeciesView&
Diagnostics::gas_var(const std::string& name) {
  auto iter = gas_vars_.find(name);
  EKAT_REQUIRE_MSG(iter != gas_vars_.end(), "Gas diagnostic variable not found!");
  return iter->second;
}

const Diagnostics::ColumnSpeciesView&
Diagnostics::gas_var(const std::string& name) const {
  auto iter = gas_vars_.find(name);
  EKAT_REQUIRE_MSG(iter != gas_vars_.end(), "Gas diagnostic variable not found!");
  return iter->second;
}

bool Diagnostics::has_modal_var(const std::string& name) const {
  return (modal_vars_.find(name) != modal_vars_.end());
}

void Diagnostics::create_modal_var(const std::string& name) {
  auto iter = modal_vars_.find(name);
  EKAT_REQUIRE_MSG(iter == modal_vars_.end(), "Modal diagnostic variable already exists!");
  modal_vars_[name] = ModalColumnView(name, num_aero_species_.size(),
                                      num_columns_, num_levels_);
}

Diagnostics::ModalColumnView&
Diagnostics::modal_var(const std::string& name) {
  auto iter = modal_vars_.find(name);
  EKAT_REQUIRE_MSG(iter != modal_vars_.end(), "Modal diagnostic variable not found!");
  return iter->second;
}

const Diagnostics::ModalColumnView&
Diagnostics::modal_var(const std::string& name) const {
  auto iter = modal_vars_.find(name);
  EKAT_REQUIRE_MSG(iter != modal_vars_.end(), "Modal diagnostic variable not found!");
  return iter->second;
}

// Interoperable C functions for providing data to Fortran.
// See haero.F90.in for details on how these functions are used.
extern "C" {

bool d_has_var_c(void* d, char name[32])
{
  Diagnostics* diags = (Diagnostics*)d;
  return diags->has_var(name);
}

void* d_var_c(void* d, char name[32])
{
  Diagnostics* diags = (Diagnostics*)d;
  auto var = diags->var(name);
  return (void*)var.data();
}

bool d_has_aerosol_var_c(void* d, char name[32])
{
  Diagnostics* diags = (Diagnostics*)d;
  return diags->has_aerosol_var(name);
}

void* d_aerosol_var_c(void* d, char name[32], int mode)
{
  Diagnostics* diags = (Diagnostics*)d;
  auto var = diags->aerosol_var(name, mode);
  return (void*)var.data();
}

bool d_has_gas_var_c(void* d, char name[32])
{
  Diagnostics* diags = (Diagnostics*)d;
  return diags->has_gas_var(name);
}

void* d_gas_var_c(void* d, char name[32])
{
  Diagnostics* diags = (Diagnostics*)d;
  auto var = diags->gas_var(name);
  return (void*)var.data();
}

bool d_has_modal_var_c(void* d, char name[32])
{
  Diagnostics* diags = (Diagnostics*)d;
  return diags->has_modal_var(name);
}

void* d_modal_var_c(void* d, char name[32])
{
  Diagnostics* diags = (Diagnostics*)d;
  auto var = diags->modal_var(name);
  return (void*)var.data();
}

} // extern "C"

} // haero

