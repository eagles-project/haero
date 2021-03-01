#include "haero/diagnostics.hpp"

namespace haero {

Diagnostics::Diagnostics(int num_aerosol_modes,
                         const std::vector<int>& num_aerosol_species,
                         int num_gases,
                         int num_levels):
  num_aero_species_(vector_to_basic_1dview(num_aerosol_species, "Diagnostics::num_aerosol_species")), 
  num_aero_populations_(0),
  num_gases_(num_gases), 
  num_levels_(num_levels) {

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
  auto T = Kokkos::create_mirror_view(num_aero_species_);
  Kokkos::deep_copy(T,num_aero_species_);
  return T[mode_index];
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

Diagnostics::TOKEN Diagnostics::find_var(const std::string& name) const {
  return get_string_to_token_vars (name);
}

Diagnostics::TOKEN Diagnostics::create_var(const std::string& name) {
  TOKEN return_val = find_var (name);
  EKAT_REQUIRE_MSG(NOT_FOUND == return_val,
    "Diagnostic variable " << name << " already exists!");
  return_val = vars_.extent(0);
  set_string_to_token_vars (name, return_val);

  int num_vert_packs = num_levels_/HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels_) {
    num_vert_packs++;
  }
  Kokkos::resize(vars_, return_val+1, num_vert_packs);
  return return_val;
}

Diagnostics::ColumnView
Diagnostics::var(const TOKEN token) {
  EKAT_KERNEL_REQUIRE_MSG(token < vars_.extent(0),
    "Diagnostic variable token not found!");
  ColumnView vars = Kokkos::subview(vars_, token, Kokkos::ALL);
  return vars;
}

//const Diagnostics::ColumnView
//Diagnostics::var(const TOKEN token) const {
//  EKAT_KERNEL_REQUIRE_MSG(token < vars_.extent(0),
//    "Diagnostic variable token not found!");
//  const ColumnView vars = Kokkos::subview(vars_, token, Kokkos::ALL);
//  return vars;
//}

Diagnostics::TOKEN Diagnostics::find_aerosol_var(const std::string& name) const {
  return get_string_to_token_aero (name);
}

Diagnostics::TOKEN Diagnostics::create_aerosol_var(const std::string& name) {
  TOKEN return_val = find_aerosol_var (name);
  EKAT_REQUIRE_MSG(NOT_FOUND == return_val,
    "Aerosol diagnostic variable " << name << " already exists!");

  // Figure out the dimensions of the view.
  int num_populations = 0;
  auto T = Kokkos::create_mirror_view(num_aero_species_);
  Kokkos::deep_copy(T,num_aero_species_);
  for (int m = 0; m < T.size(); ++m) {
    num_populations += T[m];
  }
  int num_vert_packs = num_levels_/HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels_) {
    num_vert_packs++;
  }
  return_val = aero_vars_.extent(0);
  set_string_to_token_aero(name, return_val);
  Kokkos::resize(aero_vars_, return_val+1, num_populations, num_vert_packs);
  return return_val;
}

Diagnostics::SpeciesColumnView
Diagnostics::aerosol_var(const TOKEN token) {
  EKAT_KERNEL_REQUIRE_MSG(token < aero_vars_.extent(0),
    "Aerosol diagnostic variable token not found!");
  SpeciesColumnView vars = Kokkos::subview(aero_vars_, token, Kokkos::ALL, Kokkos::ALL);
  return vars;
}

//KOKKOS_FUNCTION
//const Diagnostics::SpeciesColumnView
//Diagnostics::aerosol_var(const TOKEN token) const {
//  EKAT_KERNEL_REQUIRE_MSG(token < aero_vars_.extent(0),
//    "Aerosol diagnostic variable token  not found!");
//  const SpeciesColumnView vars = Kokkos::subview(aero_vars_, token, Kokkos::ALL, Kokkos::ALL);
//  return vars;
//}

Diagnostics::TOKEN Diagnostics::find_gas_var(const std::string& name) const {
  return get_string_to_token_gas  (name);
}

Diagnostics::TOKEN Diagnostics::create_gas_var(const std::string& name) {
  TOKEN return_val = find_gas_var (name);
  EKAT_REQUIRE_MSG(NOT_FOUND == return_val,
    "Gas diagnostic variable " << name << " already exists!");
  int num_vert_packs = num_levels_/HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels_) {
    num_vert_packs++;
  }
  return_val = gas_vars_.extent(0);
  set_string_to_token_gas(name, return_val);
  Kokkos::resize(gas_vars_, return_val+1, num_gases_, num_vert_packs);
  return return_val;
}

Diagnostics::SpeciesColumnView
Diagnostics::gas_var(const TOKEN token) {
  EKAT_KERNEL_REQUIRE_MSG(token < gas_vars_.extent(0),
    "Gas diagnostic variable token not found!");
  SpeciesColumnView vars = Kokkos::subview(gas_vars_, token, Kokkos::ALL, Kokkos::ALL);
  return vars;
}

const Diagnostics::SpeciesColumnView
Diagnostics::gas_var(const TOKEN token) const {
  EKAT_KERNEL_REQUIRE_MSG(token < gas_vars_.extent(0),
    "Gas diagnostic variable token not found!");
  const SpeciesColumnView vars = Kokkos::subview(gas_vars_, token, Kokkos::ALL, Kokkos::ALL);
  return vars;
}

Diagnostics::TOKEN Diagnostics::find_modal_var(const std::string& name) const {
  return get_string_to_token_modal(name);
}

Diagnostics::TOKEN Diagnostics::create_modal_var(const std::string& name) {
  TOKEN return_val = find_modal_var (name);
  EKAT_REQUIRE_MSG(NOT_FOUND == return_val,
    "Modal diagnostic variable " << name << " already exists!");
  int num_vert_packs = num_levels_/HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels_) {
    num_vert_packs++;
  }
  return_val = modal_vars_.extent(0);
  set_string_to_token_modal(name, return_val);
  Kokkos::resize(modal_vars_, return_val+1, num_aero_species_.size(), num_vert_packs);
  return return_val;
}

Diagnostics::ModalColumnView
Diagnostics::modal_var(const TOKEN token) {
  EKAT_KERNEL_REQUIRE_MSG(token < modal_vars_.extent(0),
    "Modal diagnostic variable token not found!");
  SpeciesColumnView vars = Kokkos::subview(modal_vars_, token, Kokkos::ALL, Kokkos::ALL);
  return vars;
}

const Diagnostics::ModalColumnView
Diagnostics::modal_var(const TOKEN token) const {
  EKAT_KERNEL_REQUIRE_MSG(token < modal_vars_.extent(0),
    "Modal diagnostic variable token not found!");
  const SpeciesColumnView vars = Kokkos::subview(modal_vars_, token, Kokkos::ALL, Kokkos::ALL);
  return vars;
}

Diagnostics::TOKEN  Diagnostics::set_string_to_token(std::map<std::string,TOKEN> &registered_strings,
                                      const std::string &name, 
                                      const TOKEN token) {
  TOKEN return_val=get_string_to_token(registered_strings, name);

  if (NOT_FOUND != token && NOT_FOUND==return_val) {
    return_val = registered_strings.size();
    const std::pair<std::string,TOKEN> value(name, return_val);
    registered_strings.insert( value );
  }
  return return_val;
}

Diagnostics::TOKEN  Diagnostics::get_string_to_token(const std::map<std::string,TOKEN> &registered_strings,
                                      const std::string &name) {
  TOKEN return_val=NOT_FOUND;
  const auto iter = registered_strings.find(name);
  if (registered_strings.end() != iter) 
    return_val = iter->second;
  return return_val;
}

namespace {
  std::map<std::string,Diagnostics::TOKEN> &registered_strings_vars() {
    static std::map<std::string,Diagnostics::TOKEN> registered_strings;
    return registered_strings;
  }
  std::map<std::string,Diagnostics::TOKEN> &registered_strings_aero() {
    static std::map<std::string,Diagnostics::TOKEN> registered_strings;
    return registered_strings;
  }
  std::map<std::string,Diagnostics::TOKEN> &registered_strings_gas() {
    static std::map<std::string,Diagnostics::TOKEN> registered_strings;
    return registered_strings;
  }
  std::map<std::string,Diagnostics::TOKEN> &registered_strings_modal() {
    static std::map<std::string,Diagnostics::TOKEN> registered_strings;
    return registered_strings;
  }
}
Diagnostics::TOKEN Diagnostics::set_string_to_token_vars(const std::string &name, const TOKEN token) {
  return set_string_to_token(registered_strings_vars(), name, token);
}
Diagnostics::TOKEN Diagnostics::set_string_to_token_aero(const std::string &name, const TOKEN token) {
  return set_string_to_token(registered_strings_aero(), name, token);
}
Diagnostics::TOKEN Diagnostics::set_string_to_token_gas (const std::string &name, const TOKEN token) {
  return set_string_to_token(registered_strings_gas (), name, token);
}
Diagnostics::TOKEN Diagnostics::set_string_to_token_modal(const std::string &name, const TOKEN token) {
  return set_string_to_token(registered_strings_modal(), name, token);
}

Diagnostics::TOKEN Diagnostics::get_string_to_token_vars(const std::string &name) const {
  return get_string_to_token(registered_strings_vars(), name);
}
Diagnostics::TOKEN Diagnostics::get_string_to_token_aero(const std::string &name) const {
  return get_string_to_token(registered_strings_aero(), name);
}
Diagnostics::TOKEN Diagnostics::get_string_to_token_gas (const std::string &name) const {
  return get_string_to_token(registered_strings_gas (), name);
}
Diagnostics::TOKEN Diagnostics::get_string_to_token_modal(const std::string &name) const {
  return get_string_to_token(registered_strings_modal(), name);
}

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
extern "C" {

Diagnostics::TOKEN d_find_var_c(void* d, const char* name)
{
  auto* diags = static_cast<Diagnostics*>(d);
  return diags->find_var(name);
}

void* d_var_c(void* d, const Diagnostics::TOKEN token)
{
  auto* diags = static_cast<Diagnostics*>(d);
  auto var = diags->var(token);
  return (void*)var.data();
}

Diagnostics::TOKEN d_find_aerosol_var_c(void* d, const char* name)
{
  auto* diags = static_cast<Diagnostics*>(d);
  return diags->find_aerosol_var(name);
}

void* d_aerosol_var_c(void* d, const Diagnostics::TOKEN token)
{
  auto* diags = static_cast<Diagnostics*>(d);
  auto var = diags->aerosol_var(token);
  return (void*)var.data();
}

Diagnostics::TOKEN d_find_gas_var_c(void* d, const char* name)
{
  auto* diags = static_cast<Diagnostics*>(d);
  return diags->find_gas_var(name);
}

void* d_gas_var_c(void* d, const Diagnostics::TOKEN token)
{
  auto* diags = static_cast<Diagnostics*>(d);
  auto var = diags->gas_var(token);
  return (void*)var.data();
}

Diagnostics::TOKEN d_find_modal_var_c(void* d, const char* name)
{
  auto* diags = static_cast<Diagnostics*>(d);
  return diags->find_modal_var(name);
}

void* d_modal_var_c(void* d, const Diagnostics::TOKEN token)
{
  auto* diags = static_cast<Diagnostics*>(d);
  auto var = diags->modal_var(token);
  return (void*)var.data();
}

} // extern "C"

} // haero

