#include "haero/prognostics.hpp"

#include "haero/tendencies.hpp"

namespace haero {

Prognostics::Prognostics(int num_aerosol_modes,
                         const std::vector<int>& num_aerosol_species,
                         int num_gases, int num_levels,
                         SpeciesColumnView int_aerosols,
                         SpeciesColumnView cld_aerosols,
                         SpeciesColumnView gases_,
                         ModalColumnView interstitial_num_concs_,
                         ModalColumnView cloudborne_num_concs_)
    : interstitial_aerosols(int_aerosols),
      cloud_aerosols(cld_aerosols),
      gases(gases_),
      interstitial_num_concs(interstitial_num_concs_),
      cloudborne_num_concs(cloudborne_num_concs_),
      num_aero_species_(vector_to_basic_1dview(
          num_aerosol_species, "Prognostics::num_aerosol_species")),
      num_aero_populations_(0),
      num_gases_(num_gases),
      num_levels_(num_levels) {
  // Count up the mode/species combinations.
  for (int m = 0; m < num_aerosol_species.size(); ++m) {
    num_aero_populations_ += num_aerosol_species[m];
  }
  EKAT_REQUIRE_MSG(
      num_aerosol_modes == num_aerosol_species.size(),
      "num_aerosol_species must be a vector of length " << num_aerosol_modes);
  int num_vert_packs = PackInfo::num_packs(num_levels_);
  const int int_aerosols_extent_0 = int_aerosols.extent(0);
  const int int_aerosols_extent_1 = int_aerosols.extent(1);
  const int cld_aerosols_extent_0 = cld_aerosols.extent(0);
  const int cld_aerosols_extent_1 = cld_aerosols.extent(1);
  const int gases_extent_0 = gases.extent(0);
  const int gases_extent_1 = gases.extent(1);
  const int interstitial_num_concs_extent_0 = interstitial_num_concs.extent(0);
  const int interstitial_num_concs_extent_1 = interstitial_num_concs.extent(1);
  const int cloudborne_num_concs_extent_0 = cloudborne_num_concs.extent(0);
  const int cloudborne_num_concs_extent_1 = cloudborne_num_concs.extent(1);
  EKAT_REQUIRE_MSG(
      int_aerosols_extent_0 == num_aero_populations_,
      "int_aerosols must have extent(0) == " << num_aero_populations_);
  EKAT_REQUIRE_MSG(int_aerosols_extent_1 == num_vert_packs,
                   "int_aerosols must have extent(1) == " << num_vert_packs);
  EKAT_REQUIRE_MSG(
      cld_aerosols_extent_0 == num_aero_populations_,
      "cld_aerosols must have extent(0) == " << num_aero_populations_);
  EKAT_REQUIRE_MSG(cld_aerosols_extent_1 == num_vert_packs,
                   "int_aerosols must have extent(1) == " << num_vert_packs);
  EKAT_REQUIRE_MSG(gases_extent_0 == num_gases_,
                   "gases must have extent(0) == " << num_gases_);
  EKAT_REQUIRE_MSG(gases_extent_1 == num_vert_packs,
                   "gases must have extent(1) == " << num_vert_packs);
  EKAT_REQUIRE_MSG(
      interstitial_num_concs_extent_0 == num_aerosol_modes,
      "interstitial_num_concs must have extent(0) == " << num_aerosol_modes);
  EKAT_REQUIRE_MSG(
      interstitial_num_concs_extent_1 == num_vert_packs,
      "interstitial_num_concs must have extent(1) == " << num_vert_packs);
  EKAT_REQUIRE_MSG(
      cloudborne_num_concs_extent_0 == num_aerosol_modes,
      "cloudborne_num_concs must have extent(0) == " << num_aerosol_modes);
  EKAT_REQUIRE_MSG(
      cloudborne_num_concs_extent_1 == num_vert_packs,
      "cloudborne_num_concs must have extent(1) == " << num_vert_packs);
}

Prognostics::~Prognostics() {}

int Prognostics::num_aerosol_modes() const { return num_aero_species_.size(); }

int Prognostics::num_aerosol_species(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < num_aero_species_.size());
  return num_aero_species_[mode_index];
}

int Prognostics::num_aerosol_populations() const {
  return num_aero_populations_;
}

int Prognostics::num_gases() const { return num_gases_; }

int Prognostics::num_levels() const { return num_levels_; }

void Prognostics::scale_and_add(Real scale_factor,
                                const Tendencies& tendencies) {
  int num_vert_packs = PackInfo::num_packs(num_levels_);
  Kokkos::parallel_for(
      "Prognostics::scale_and_add (aerosols)", num_vert_packs,
      KOKKOS_LAMBDA(const int k) {
        for (int p = 0; p < num_aero_populations_; ++p) {
          interstitial_aerosols(p, k) +=
              scale_factor * tendencies.interstitial_aerosols(p, k);
          cloud_aerosols(p, k) +=
              scale_factor * tendencies.cloud_aerosols(p, k);
        }
      });

  Kokkos::parallel_for(
      "Prognostics::scale_and_add (gas species)", num_vert_packs,
      KOKKOS_LAMBDA(const int k) {
        for (int g = 0; g < num_gases_; ++g) {
          gases(g, k) += scale_factor * tendencies.gases(g, k);
        }
      });

  int num_modes = interstitial_num_concs.extent(0);
  Kokkos::parallel_for(
      "Prognostics::scale_and_add (modal num concs)", num_vert_packs,
      KOKKOS_LAMBDA(const int k) {
        for (int m = 0; m < num_modes; ++m) {
          interstitial_num_concs(m, k) +=
              scale_factor * tendencies.interstitial_num_concs(m, k);
          cloudborne_num_concs(m, k) +=
              scale_factor * tendencies.cloudborne_num_concs(m, k);
        }
      });
}

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
extern "C" {

void* p_int_aero_mix_frac_c(void* p) {
  auto* progs = static_cast<Prognostics*>(p);
  auto mix_fracs = progs->interstitial_aerosols;
  return (void*)mix_fracs.data();
}

void* p_cld_aero_mix_frac_c(void* p) {
  auto* progs = static_cast<Prognostics*>(p);
  auto mix_fracs = progs->cloud_aerosols;
  return (void*)mix_fracs.data();
}

void* p_gases_c(void* p) {
  auto* progs = static_cast<Prognostics*>(p);
  auto mix_fracs = progs->gases;
  return (void*)mix_fracs.data();
}

void* p_interstitial_num_concs_c(void* p) {
  auto* progs = static_cast<Prognostics*>(p);
  auto num_concs = progs->interstitial_num_concs;
  return (void*)num_concs.data();
}

void* p_cloudborne_num_concs_c(void* p) {
  auto* progs = static_cast<Prognostics*>(p);
  auto num_concs = progs->cloudborne_num_concs;
  return (void*)num_concs.data();
}

}  // extern "C"

}  // namespace haero
