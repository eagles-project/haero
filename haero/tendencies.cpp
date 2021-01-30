#include "haero/tendencies.hpp"

namespace haero {

Tendencies::Tendencies(const Prognostics& prognostics) {
  int num_modes = prognostics.num_aerosol_modes();
  int num_mode_species_pairs = 0;
  for (int m = 0; m < num_modes; ++m) {
    num_mode_species_pairs += prognostics.num_aerosol_species(m);
  }
  int num_levels = prognostics.num_levels();
  int num_vert_packs = num_levels/HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels) {
    num_vert_packs++;
  }
  auto int_view_name = std::string("d/dt[") +
                       prognostics.interstitial_aerosols().label() +
                       std::string(")]");
  int_aero_species_ = SpeciesColumnView(int_view_name, num_mode_species_pairs,
                                        num_vert_packs);
  auto cld_view_name = std::string("d/dt[") +
                       prognostics.cloudborne_aerosols().label() +
                       std::string(")]");
  cld_aero_species_ = SpeciesColumnView(cld_view_name, num_mode_species_pairs,
                                        num_vert_packs);


  int num_gases = prognostics.num_gases();
  auto gas_view_name = std::string("d/dt[") +
                       prognostics.gases().label() +
                       std::string(")]");
  gases_ = SpeciesColumnView(gas_view_name, num_gases, num_vert_packs);

  auto n_view_name = std::string("d/dt[") +
                     prognostics.modal_num_concs().label() +
                     std::string(")]");
  modal_num_concs_ = ModalColumnView(n_view_name, num_modes, num_vert_packs);
}

Tendencies::~Tendencies() {
}

Tendencies::SpeciesColumnView&
Tendencies::interstitial_aerosols() {
  return int_aero_species_;
}

const Tendencies::SpeciesColumnView&
Tendencies::interstitial_aerosols() const {
  return int_aero_species_;
}

Tendencies::SpeciesColumnView&
Tendencies::cloudborne_aerosols() {
  return cld_aero_species_;
}

const Tendencies::SpeciesColumnView&
Tendencies::cloudborne_aerosols() const {
  return cld_aero_species_;
}

Tendencies::SpeciesColumnView& Tendencies::gases() {
  return gases_;
}

const Tendencies::SpeciesColumnView& Tendencies::gases() const {
  return gases_;
}

Tendencies::ModalColumnView& Tendencies::modal_num_concs() {
  return modal_num_concs_;
}

const Tendencies::ModalColumnView& Tendencies::modal_num_concs() const {
  return modal_num_concs_;
}

Tendencies& Tendencies::scale(Real factor) {
  EKAT_REQUIRE_MSG(false, "Tendencies::scale is not yet implemented!");
  return *this;
}

void Tendencies::accumulate(const Tendencies& tendencies) {
  EKAT_REQUIRE_MSG(false, "Tendencies::accumulate is not yet implemented!");
}

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
extern "C" {

void* t_int_aero_mix_frac_c(void* t)
{
  Tendencies* tends = (Tendencies*)t;
  auto& mix_fracs = tends->interstitial_aerosols();
  return (void*)mix_fracs.data();
}

void* t_cld_aero_mix_frac_c(void* t)
{
  Tendencies* tends = (Tendencies*)t;
  auto& mix_fracs = tends->cloudborne_aerosols();
  return (void*)mix_fracs.data();
}

void* t_gases_c(void* t)
{
  Tendencies* tends = (Tendencies*)t;
  auto& mole_fracs = tends->gases();
  return (void*)mole_fracs.data();
}

void* t_modal_num_concs_c(void* t)
{
  Tendencies* tends = (Tendencies*)t;
  auto& num_densities = tends->modal_num_concs();
  return (void*)num_densities.data();
}

} // extern "C"

} // haero

