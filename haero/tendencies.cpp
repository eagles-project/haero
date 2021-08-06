#include "haero/tendencies.hpp"

#include "haero/view_pack_helpers.hpp"
#include "kokkos/Kokkos_Core.hpp"

namespace haero {

Tendencies::Tendencies(const Prognostics& prognostics) {
  int num_aerosol_populations = prognostics.num_aerosol_populations();
  int num_levels = prognostics.num_levels();
  int num_vert_packs = num_levels / HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels) {
    num_vert_packs++;
  }
  auto int_view_name = std::string("d/dt[") +
                       prognostics.interstitial_aerosols.label() +
                       std::string(")]");
  interstitial_aerosols =
      SpeciesColumnView(int_view_name, num_aerosol_populations, num_vert_packs);
  Kokkos::deep_copy(interstitial_aerosols, PackType(0));
  auto cld_view_name = std::string("d/dt[") +
                       prognostics.cloud_aerosols.label() + std::string(")]");
  cloud_aerosols =
      SpeciesColumnView(cld_view_name, num_aerosol_populations, num_vert_packs);
  Kokkos::deep_copy(cloud_aerosols, PackType(0));

  int num_gases = prognostics.num_gases();
  auto gas_view_name =
      std::string("d/dt[") + prognostics.gases.label() + std::string(")]");
  gases = SpeciesColumnView(gas_view_name, num_gases, num_vert_packs);
  Kokkos::deep_copy(gases, PackType(0));

  auto n_view_name = std::string("d/dt[") +
                     prognostics.interstitial_num_mix_ratios.label() +
                     std::string(")]");
  int num_modes = prognostics.num_aerosol_modes();
  interstitial_num_mix_ratios =
      ModeColumnView(n_view_name, num_modes, num_vert_packs);
  Kokkos::deep_copy(interstitial_num_mix_ratios, PackType(0));

  cloud_num_mix_ratios = ModeColumnView(n_view_name, num_modes, num_vert_packs);
  Kokkos::deep_copy(cloud_num_mix_ratios, PackType(0));
}

Tendencies::~Tendencies() {}

const Tendencies& Tendencies::scale(Real factor) const {
  int num_populations = interstitial_aerosols.extent(0);
  int num_levels = interstitial_aerosols.extent(1);
  int num_vert_packs = PackInfo::num_packs(num_levels);
  int num_gases = gases.extent(0);
  int num_modes = interstitial_num_mix_ratios.extent(0);

  // Scale aerosol species mixing ratios
  Kokkos::parallel_for(
      "Tendencies::scale (aero mixіng ratios)", num_vert_packs,
      KOKKOS_LAMBDA(const int k) {
        for (int p = 0; p < num_populations; ++p) {
          interstitial_aerosols(p, k) *= factor;
          cloud_aerosols(p, k) *= factor;
        }
      });

  // Scale gas mole fractions.
  Kokkos::parallel_for(
      "tendencies::scale (gas species)", num_vert_packs,
      KOKKOS_LAMBDA(const int k) {
        for (int g = 0; g < num_gases; ++g) {
          gases(g, k) *= factor;
        }
      });

  // Scale interstitial number densities.
  Kokkos::parallel_for(
      "tendencies::scale (interstitial num mix_ratios)", num_vert_packs,
      KOKKOS_LAMBDA(const int k) {
        for (int m = 0; m < num_modes; ++m) {
          interstitial_num_mix_ratios(m, k) *= factor;
        }
      });

  // Scale cloud borne number densities.
  Kokkos::parallel_for(
      "tendencies::scale (cloud num mix_ratios)", num_vert_packs,
      KOKKOS_LAMBDA(const int k) {
        for (int m = 0; m < num_modes; ++m) {
          cloud_num_mix_ratios(m, k) *= factor;
        }
      });

  return *this;
}

void Tendencies::accumulate(const Tendencies& tendencies) {
  int num_populations = interstitial_aerosols.extent(0);
  int num_levels = interstitial_aerosols.extent(1);
  int num_vert_packs = PackInfo::num_packs(num_levels);
  int num_gases = gases.extent(0);
  int num_modes = interstitial_num_mix_ratios.extent(0);

  // Accumulate aerosol species mixing ratios
  Kokkos::parallel_for(
      "Tendencies::scale (aero mixіng ratios)", num_vert_packs,
      KOKKOS_LAMBDA(const int k) {
        for (int p = 0; p < num_populations; ++p) {
          interstitial_aerosols(p, k) += tendencies.interstitial_aerosols(p, k);
          cloud_aerosols(p, k) += tendencies.cloud_aerosols(p, k);
        }
      });

  // Scale gas mole fractions.
  Kokkos::parallel_for(
      "tendencies::scale (gas species)", num_vert_packs,
      KOKKOS_LAMBDA(const int k) {
        for (int g = 0; g < num_gases; ++g) {
          gases(g, k) += tendencies.gases(g, k);
        }
      });

  // Scale mode number densities.
  Kokkos::parallel_for(
      "tendencies::scale (interstitial num mix_ratios)", num_vert_packs,
      KOKKOS_LAMBDA(const int k) {
        for (int m = 0; m < num_modes; ++m) {
          interstitial_num_mix_ratios(m, k) +=
              tendencies.interstitial_num_mix_ratios(m, k);
        }
      });

  // Scale mode number densities.
  Kokkos::parallel_for(
      "tendencies::scale (cloud num mix_ratios)", num_vert_packs,
      KOKKOS_LAMBDA(const int k) {
        for (int m = 0; m < num_modes; ++m) {
          cloud_num_mix_ratios(m, k) += tendencies.cloud_num_mix_ratios(m, k);
        }
      });
}

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
extern "C" {

void* t_int_aero_mix_frac_c(void* t) {
  Tendencies* tends = (Tendencies*)t;
  auto mix_fracs = tends->interstitial_aerosols;
  return (void*)mix_fracs.data();
}

void* t_cld_aero_mix_frac_c(void* t) {
  Tendencies* tends = (Tendencies*)t;
  auto mix_fracs = tends->cloud_aerosols;
  return (void*)mix_fracs.data();
}

void* t_gases_c(void* t) {
  Tendencies* tends = (Tendencies*)t;
  auto mix_fracs = tends->gases;
  return (void*)mix_fracs.data();
}

void* t_interstitial_num_mix_ratios_c(void* t) {
  Tendencies* tends = (Tendencies*)t;
  auto num_mix_ratios = tends->interstitial_num_mix_ratios;
  return (void*)num_mix_ratios.data();
}

void* t_cloud_num_mix_ratios_c(void* t) {
  Tendencies* tends = (Tendencies*)t;
  auto num_mix_ratios = tends->cloud_num_mix_ratios;
  return (void*)num_mix_ratios.data();
}

}  // extern "C"

}  // namespace haero
