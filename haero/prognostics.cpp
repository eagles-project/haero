#include "haero/prognostics.hpp"

namespace haero {

Prognostics::Prognostics(int num_aerosol_modes,
                         const std::vector<int>& num_aerosol_species,
                         int num_gases,
                         int num_levels,
                         Kokkos::View<PackType**>& int_aerosols,
                         Kokkos::View<PackType**>& cld_aerosols,
                         Kokkos::View<PackType**>& gases,
                         Kokkos::View<PackType**>& modal_num_concs):
  num_aero_species_(num_aerosol_species), num_gases_(num_gases),
  num_levels_(num_levels), int_aero_species_(int_aerosols),
  cld_aero_species_(cld_aerosols), gases_(gases),
  modal_num_concs_(modal_num_concs) {
#ifndef NDEBUG
  EKAT_ASSERT_MSG(num_aerosol_modes == num_aerosol_species.size(),
                  "num_aerosol_species must be a vector of length " << num_aerosol_modes);
  // Count up the mode/species combinations to validate the extents of the
  // aerosol views.
  int n = 0;
  for (int m = 0; m < num_aerosol_species.size(); ++m) {
    n += num_aerosol_species[m];
  }
  int num_vert_packs = num_levels_/HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels_) {
    num_vert_packs++;
  }
  EKAT_ASSERT_MSG(int_aerosols.extent(0) == n,
                  "int_aerosols must have extent(0) == " << n);
  EKAT_ASSERT_MSG(int_aerosols.extent(1) == num_vert_packs,
                  "int_aerosols must have extent(1) == " << num_vert_packs);
  EKAT_ASSERT_MSG(cld_aerosols.extent(0) == n,
                  "cld_aerosols must have extent(0) == " << n);
  EKAT_ASSERT_MSG(cld_aerosols.extent(1) == num_vert_packs,
                  "int_aerosols must have extent(1) == " << num_vert_packs);
  EKAT_ASSERT_MSG(gases.extent(0) == num_gases_,
                  "gases must have extent(0) == " << num_gases_);
  EKAT_ASSERT_MSG(gases.extent(0) == num_vert_packs,
                  "gases must have extent(1) == " << num_vert_packs);
  EKAT_ASSERT_MSG(modal_num_concs.extent(0) == num_aerosol_modes,
                  "modal_num_concs must have extent(0) == " << num_aerosol_modes);
  EKAT_ASSERT_MSG(modal_num_concs.extent(0) == num_vert_packs,
                  "modal_num_concs must have extent(1) == " << num_vert_packs);
#endif
}

Prognostics::~Prognostics() {
}

int Prognostics::num_aerosol_modes() const {
  return num_aero_species_.size();
}

int Prognostics::num_aerosol_species(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < num_aero_species_.size());
  return num_aero_species_[mode_index];
}

int Prognostics::num_gases() const {
  return num_gases_;
}

int Prognostics::num_levels() const {
  return num_levels_;
}

Prognostics::SpeciesColumnView&
Prognostics::interstitial_aerosols() {
  return int_aero_species_;
}

const Prognostics::SpeciesColumnView&
Prognostics::interstitial_aerosols() const {
  return int_aero_species_;
}

Prognostics::SpeciesColumnView& Prognostics::cloudborne_aerosols() {
  return cld_aero_species_;
}

const Prognostics::SpeciesColumnView&
Prognostics::cloudborne_aerosols() const {
  return cld_aero_species_;
}

Prognostics::SpeciesColumnView& Prognostics::gases() {
  return gases_;
}

const Prognostics::SpeciesColumnView& Prognostics::gases() const {
  return gases_;
}

Prognostics::ModalColumnView& Prognostics::modal_num_concs() {
  return modal_num_concs_;
}

const Prognostics::ModalColumnView& Prognostics::modal_num_concs() const {
  return modal_num_concs_;
}

void Prognostics::scale_and_add(Real scale_factor,
                                const Tendencies& tendencies) {
  EKAT_REQUIRE_MSG(false, "scale_and_add() is not yet implemented!");
}

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
extern "C" {

void* p_int_aero_mix_frac_c(void* p)
{
  auto* progs = static_cast<Prognostics*>(p);
  auto& mix_fracs = progs->interstitial_aerosols();
  return (void*)mix_fracs.data();
}

void* p_cld_aero_mix_frac_c(void* p)
{
  auto* progs = static_cast<Prognostics*>(p);
  auto& mix_fracs = progs->cloudborne_aerosols();
  return (void*)mix_fracs.data();
}

void* p_gases_c(void* p)
{
  auto* progs = static_cast<Prognostics*>(p);
  auto& mix_fracs = progs->gases();
  return (void*)mix_fracs.data();
}

void* p_modal_num_concs_c(void* p)
{
  auto* progs = static_cast<Prognostics*>(p);
  auto& num_concs = progs->modal_num_concs();
  return (void*)num_concs.data();
}

} // extern "C"

} // haero

