#include "haero/diagnostics.hpp"

namespace haero {

Diagnostics::Diagnostics(int num_aerosol_modes,
                         const std::vector<int>& num_aerosol_species,
                         int num_gases,
                         int num_levels):
  num_aero_species_(vector_to_basic_1dview(num_aerosol_species, "Diagnostics::num_aerosol_species")), 
  num_aero_populations_(0),
  num_gases_(num_gases), num_levels_(num_levels) {
  EKAT_ASSERT_MSG(num_aerosol_modes == num_aerosol_species.size(),
                  "num_aerosol_species must be a vector of length " << num_aerosol_modes);
  for (int m = 0; m < num_aerosol_modes; ++m) {
    num_aero_populations_ += num_aerosol_species[m];
  }
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

int Diagnostics::num_aerosol_populations() const {
  return num_aero_populations_;
}

int Diagnostics::num_gases() const {
  return num_gases_;
}

int Diagnostics::num_levels() const {
  return num_levels_;
}

int Diagnostics::has_var(const std::string& name) const {
  return get_string_to_token_vars (name);
}

int Diagnostics::create_var(const std::string& name) {
  int return_val = has_var (name);
  EKAT_REQUIRE_MSG(NOT_FOUND == return_val,
    "Diagnostic variable " << name << " already exists!");
  return_val = vars_.extent(0);
  set_string_to_token_vars (name, return_val);
  Kokkos::resize(vars_, return_val+1, num_levels_);
  return return_val;
}

Diagnostics::ColumnView
Diagnostics::var(const int token) {
  EKAT_REQUIRE_MSG(token < vars_.extent(0),
    "Diagnostic variable " << token << " not found!");
  ColumnView vars = Kokkos::subview(vars_, token, Kokkos::ALL);
  return vars;
}

const Diagnostics::ColumnView
Diagnostics::var(const int token) const {
  EKAT_REQUIRE_MSG(token < vars_.extent(0),
    "Diagnostic variable " << token << " not found!");
  const ColumnView vars = Kokkos::subview(vars_, token, Kokkos::ALL);
  return vars;
}

int Diagnostics::has_aerosol_var(const std::string& name) const {
  return get_string_to_token_aero (name);
}

void Diagnostics::create_aerosol_var(const std::string& name) {
  auto iter = aero_vars_.find(name);
  EKAT_REQUIRE_MSG(iter == aero_vars_.end(),
    "Aerosol diagnostic variable " << name << " already exists!");

  // Figure out the dimensions of the view.
  int num_populations = 0;
  for (int m = 0; m < num_aero_species_.size(); ++m) {
    num_populations += num_aero_species_[m];
  }
  int num_vert_packs = num_levels_/HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels_) {
    num_vert_packs++;
  }
  aero_vars_[name] = SpeciesColumnView(name, num_populations, num_vert_packs);
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

int Diagnostics::has_gas_var(const std::string& name) const {
  return get_string_to_token_gas  (name);
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

int Diagnostics::has_modal_var(const std::string& name) const {
  return get_string_to_token_modal(name);
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

int  Diagnostics::set_string_to_token(std::map<std::string,int> &registered_strings,
                                      const std::string &name, 
                                      const int token) {
  int return_val=get_string_to_token(registered_strings, name);

  if (NOT_FOUND != token && NOT_FOUND==return_val) {
    return_val = registered_strings.size();
    const std::pair<std::string,int> value(name, return_val);
    registered_strings.insert( value );
  }
  return return_val;
}

int  Diagnostics::get_string_to_token(const std::map<std::string,int> &registered_strings,
                                      const std::string &name) {
  int return_val=NOT_FOUND;
  const auto iter = registered_strings.find(name);
  if (registered_strings.end() != iter) 
    return_val = iter->second;
  return return_val;
}

int Diagnostics::set_string_to_token_vars(const std::string &name, const int token) {
  return set_string_to_token(registered_strings_vars, name, token);
}
int Diagnostics::set_string_to_token_aero(const std::string &name, const int token) {
  return set_string_to_token(registered_strings_aero, name, token);
}
int Diagnostics::set_string_to_token_gas (const std::string &name, const int token) {
  return set_string_to_token(registered_strings_gas , name, token);
}
int Diagnostics::set_string_to_token_modal(const std::string &name, const int token) {
  return set_string_to_token(registered_strings_modal, name, token);
}

int Diagnostics::get_string_to_token_vars(const std::string &name) const {
  return get_string_to_token(registered_strings_vars, name);
}
int Diagnostics::get_string_to_token_aero(const std::string &name) const {
  return get_string_to_token(registered_strings_aero, name);
}
int Diagnostics::get_string_to_token_gas (const std::string &name) const {
  return get_string_to_token(registered_strings_gas , name);
}
int Diagnostics::get_string_to_token_modal(const std::string &name) const {
  return get_string_to_token(registered_strings_modal, name);
}

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
extern "C" {

int d_has_var_c(void* d, const char* name)
{
  auto* diags = static_cast<Diagnostics*>(d);
  return diags->has_var(name);
}

void* d_var_c(void* d, const int token)
{
  auto* diags = static_cast<Diagnostics*>(d);
  auto var = diags->var(token);
  return (void*)var.data();
}

int d_has_aerosol_var_c(void* d, const char* name)
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

int d_has_gas_var_c(void* d, const char* name)
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

int d_has_modal_var_c(void* d, const char* name)
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

