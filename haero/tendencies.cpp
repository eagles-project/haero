#include "haero/tendencies.hpp"

#include "haero/view_pack_helpers.hpp"
#include "kokkos/Kokkos_Core.hpp"

namespace haero {

Tendencies::Tendencies(const Prognostics& prognostics) {
  int num_aerosol_populations = prognostics.num_aerosol_populations();
  int num_levels = prognostics.num_levels();
  int num_vert_packs = num_levels/HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels) {
    num_vert_packs++;
  }
  auto int_view_name = std::string("d/dt[") +
                       prognostics.interstitial_aerosols().label() +
                       std::string(")]");
  int_aero_species_ = SpeciesColumnView(int_view_name, num_aerosol_populations,
                                        num_vert_packs);
  auto cld_view_name = std::string("d/dt[") +
                       prognostics.cloudborne_aerosols().label() +
                       std::string(")]");
  cld_aero_species_ = SpeciesColumnView(cld_view_name, num_aerosol_populations,
                                        num_vert_packs);

  int num_gases = prognostics.num_gases();
  auto gas_view_name = std::string("d/dt[") +
                       prognostics.gases().label() +
                       std::string(")]");
  gases_ = SpeciesColumnView(gas_view_name, num_gases, num_vert_packs);

  auto n_view_name = std::string("d/dt[") +
                     prognostics.modal_num_concs().label() +
                     std::string(")]");
  int num_modes = prognostics.num_aerosol_modes();
  modal_num_concs_ = ModalColumnView(n_view_name, num_modes, num_vert_packs);
}

Tendencies::~Tendencies() {
}

//Tendencies::SpeciesColumnView
//Tendencies::interstitial_aerosols() {
//  return int_aero_species_;
//}

const Tendencies::SpeciesColumnView
Tendencies::interstitial_aerosols() const {
  return int_aero_species_;
}

Tendencies::SpeciesColumnView
Tendencies::cloudborne_aerosols() {
  return cld_aero_species_;
}

const Tendencies::SpeciesColumnView
Tendencies::cloudborne_aerosols() const {
  return cld_aero_species_;
}

Tendencies::SpeciesColumnView Tendencies::gases() {
  return gases_;
}

const Tendencies::SpeciesColumnView Tendencies::gases() const {
  return gases_;
}

Tendencies::ModalColumnView Tendencies::modal_num_concs() {
  return modal_num_concs_;
}

const Tendencies::ModalColumnView Tendencies::modal_num_concs() const {
  return modal_num_concs_;
}

Tendencies& Tendencies::scale(Real factor) {
  int num_populations = int_aero_species_.extent(0);
  int num_levels = int_aero_species_.extent(1);
  int num_vert_packs = PackInfo::num_packs(num_levels);
  int num_gases = gases_.extent(0);
  int num_modes = modal_num_concs_.extent(0);

  // Scale aerosol species mixing ratios
  Kokkos::parallel_for("Tendencies::scale (aero mixіng ratios)", num_vert_packs,
    KOKKOS_LAMBDA (const int k) {
      for (int p = 0; p < num_populations; ++p) {
        int_aero_species_(p, k) *= factor;
        cld_aero_species_(p, k) *= factor;
      }
    });

  // Scale gas mole fractions.
  Kokkos::parallel_for("tendencies::scale (gas species)", num_vert_packs,
    KOKKOS_LAMBDA (const int k) {
      for (int g = 0; g < num_gases; ++g) {
        gases_(g, k) *= factor;
      }
    });

  // Scale modal number densities.
  Kokkos::parallel_for("tendencies::scale (modal num concs)", num_vert_packs,
    KOKKOS_LAMBDA (const int k) {
      for (int m = 0; m < num_modes; ++m) {
        modal_num_concs_(m, k) *= factor;
      }
    });

  return *this;
}

void Tendencies::accumulate(const Tendencies& tendencies) {
  int num_populations = int_aero_species_.extent(0);
  int num_levels = int_aero_species_.extent(1);
  int num_vert_packs = PackInfo::num_packs(num_levels);
  int num_gases = gases_.extent(0);
  int num_modes = modal_num_concs_.extent(0);

  // Accumulate aerosol species mixing ratios
  Kokkos::parallel_for("Tendencies::scale (aero mixіng ratios)", num_vert_packs,
    KOKKOS_LAMBDA (const int k) {
      for (int p = 0; p < num_populations; ++p) {
        int_aero_species_(p, k) += tendencies.int_aero_species_(p, k);
        cld_aero_species_(p, k) += tendencies.cld_aero_species_(p, k);
      }
    });

  // Scale gas mole fractions.
  Kokkos::parallel_for("tendencies::scale (gas species)", num_vert_packs,
    KOKKOS_LAMBDA (const int k) {
      for (int g = 0; g < num_gases; ++g) {
        gases_(g, k) += tendencies.gases_(g, k);
      }
    });

  // Scale modal number densities.
  Kokkos::parallel_for("tendencies::scale (modal num concs)", num_vert_packs,
    KOKKOS_LAMBDA (const int k) {
      for (int m = 0; m < num_modes; ++m) {
        modal_num_concs_(m, k) += tendencies.modal_num_concs_(m, k);
      }
    });

}

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
extern "C" {

void* t_int_aero_mix_frac_c(void* t)
{
  Tendencies* tends = (Tendencies*)t;
  auto mix_fracs = tends->interstitial_aerosols();
  return (void*)mix_fracs.data();
}

void* t_cld_aero_mix_frac_c(void* t)
{
  Tendencies* tends = (Tendencies*)t;
  auto mix_fracs = tends->cloudborne_aerosols();
  return (void*)mix_fracs.data();
}

void* t_gases_c(void* t)
{
  Tendencies* tends = (Tendencies*)t;
  auto mix_fracs = tends->gases();
  return (void*)mix_fracs.data();
}

void* t_modal_num_concs_c(void* t)
{
  Tendencies* tends = (Tendencies*)t;
  auto num_concs = tends->modal_num_concs();
  return (void*)num_concs.data();
}

} // extern "C"

} // haero

