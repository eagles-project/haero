#include "haero/prognostics.hpp"

#include "haero/tendencies.hpp"

namespace haero {

Prognostics::Prognostics(const ModalAerosolConfig& config, int num_levels)
    : interstitial_aerosols("interstitial aerosols",
                            config.num_aerosol_populations,
                            PackInfo::num_packs(num_levels)),
      cloud_aerosols("cloudborne aerosols", config.num_aerosol_populations,
                     PackInfo::num_packs(num_levels)),
      interstitial_num_mix_ratios("interstitial aerosol number mix ratios",
                                  config.num_aerosol_modes(),
                                  PackInfo::num_packs(num_levels)),
      cloud_num_mix_ratios("cloudborne aerosol number mix ratios",
                           config.num_aerosol_modes(),
                           PackInfo::num_packs(num_levels)),
      gases("gases", config.num_gases(), PackInfo::num_packs(num_levels)),
      num_aero_species_("Prognostics::num_aerosol_species",
                        config.aerosol_species.size()),
      num_aero_populations_(config.num_aerosol_populations),
      num_gases_(config.num_gases()),
      num_levels_(num_levels) {}

Prognostics::Prognostics(const ModalAerosolConfig& config, int num_levels,
                         SpeciesColumnView int_aerosols,
                         SpeciesColumnView cld_aerosols,
                         ModeColumnView int_mode_num_mix_ratios_,
                         ModeColumnView cld_mode_num_mix_ratios_,
                         SpeciesColumnView gases_)
    : interstitial_aerosols(int_aerosols),
      cloud_aerosols(cld_aerosols),
      interstitial_num_mix_ratios(int_mode_num_mix_ratios_),
      cloud_num_mix_ratios(cld_mode_num_mix_ratios_),
      gases(gases_),
      num_aero_species_("Prognostics::num_aerosol_species",
                        config.aerosol_species.size()),
      num_aero_populations_(config.num_aerosol_populations),
      num_gases_(config.num_gases()),
      num_levels_(num_levels) {
  int num_aerosol_modes = config.num_aerosol_modes();
  int num_vert_packs = PackInfo::num_packs(num_levels_);
  const int int_aerosols_extent_0 = int_aerosols.extent(0);
  const int int_aerosols_extent_1 = int_aerosols.extent(1);
  const int cld_aerosols_extent_0 = cld_aerosols.extent(0);
  const int cld_aerosols_extent_1 = cld_aerosols.extent(1);
  const int gases_extent_0 = gases.extent(0);
  const int gases_extent_1 = gases.extent(1);
  const int interstitial_num_mix_ratios_extent_0 =
      interstitial_num_mix_ratios.extent(0);
  const int interstitial_num_mix_ratios_extent_1 =
      interstitial_num_mix_ratios.extent(1);
  const int cloud_num_mix_ratios_extent_0 = cloud_num_mix_ratios.extent(0);
  const int cloud_num_mix_ratios_extent_1 = cloud_num_mix_ratios.extent(1);
  EKAT_REQUIRE_MSG(
      int_aerosols_extent_0 == num_aero_populations_,
      "int_aerosols must have extent(0) == " << num_aero_populations_);
  EKAT_REQUIRE_MSG(int_aerosols_extent_1 == num_vert_packs,
                   "int_aerosols must have extent(1) == " << num_vert_packs);
  EKAT_REQUIRE_MSG(
      cld_aerosols_extent_0 == num_aero_populations_,
      "cld_aerosols must have extent(0) == " << num_aero_populations_);
  EKAT_REQUIRE_MSG(cld_aerosols_extent_1 == num_vert_packs,
                   "cld_aerosols must have extent(1) == " << num_vert_packs);
  EKAT_REQUIRE_MSG(gases_extent_0 == num_gases_,
                   "gases must have extent(0) == " << num_gases_);
  EKAT_REQUIRE_MSG(gases_extent_1 == num_vert_packs,
                   "gases must have extent(1) == " << num_vert_packs);
  EKAT_REQUIRE_MSG(interstitial_num_mix_ratios_extent_0 == num_aerosol_modes,
                   "interstitial_num_mix_ratios must have extent(0) == "
                       << num_aerosol_modes);
  EKAT_REQUIRE_MSG(
      interstitial_num_mix_ratios_extent_1 == num_vert_packs,
      "interstitial_num_mix_ratios must have extent(1) == " << num_vert_packs);
  EKAT_REQUIRE_MSG(
      cloud_num_mix_ratios_extent_0 == num_aerosol_modes,
      "cloud_num_mix_ratios must have extent(0) == " << num_aerosol_modes);
  EKAT_REQUIRE_MSG(
      cloud_num_mix_ratios_extent_1 == num_vert_packs,
      "cloud_num_mix_ratios must have extent(1) == " << num_vert_packs);
}

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

  int num_modes = interstitial_num_mix_ratios.extent(0);
  Kokkos::parallel_for(
      "Prognostics::scale_and_add (modal num mix_ratios)", num_vert_packs,
      KOKKOS_LAMBDA(const int k) {
        for (int m = 0; m < num_modes; ++m) {
          interstitial_num_mix_ratios(m, k) +=
              scale_factor * tendencies.interstitial_num_mix_ratios(m, k);
          cloud_num_mix_ratios(m, k) +=
              scale_factor * tendencies.cloud_num_mix_ratios(m, k);
        }
      });
}

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
extern "C" {

int p_num_levels_c(void* p) {
  auto* progs = static_cast<Prognostics*>(p);
  return progs->num_levels();
}

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

void* p_interstitial_num_mix_ratios_c(void* p) {
  auto* progs = static_cast<Prognostics*>(p);
  auto num_mix_ratios = progs->interstitial_num_mix_ratios;
  return (void*)num_mix_ratios.data();
}

void* p_cloud_num_mix_ratios_c(void* p) {
  auto* progs = static_cast<Prognostics*>(p);
  auto num_mix_ratios = progs->cloud_num_mix_ratios;
  return (void*)num_mix_ratios.data();
}

}  // extern "C"

}  // namespace haero