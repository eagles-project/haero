#include "haero/atmosphere.hpp"

namespace haero {

Atmosphere::Atmosphere(int num_levels, Real pblh)
    : temperature("temperature", PackInfo::num_packs(num_levels)),
      pressure("pressure", PackInfo::num_packs(num_levels)),
      vapor_mixing_ratio("vapor mixing ratio", PackInfo::num_packs(num_levels)),
      height("height", PackInfo::num_packs(num_levels)),
      hydrostatic_dp("hydrostatic_dp", PackInfo::num_packs(num_levels)),
      planetary_boundary_height(pblh),
      num_levels_(num_levels) {
  EKAT_REQUIRE_MSG(num_levels > 0,
                   "Number of vertical levels must be positive");
  EKAT_REQUIRE_MSG(pblh >= 0.0,
                   "Planetary boundary height must be non-negative");
}

Atmosphere::Atmosphere(int num_levels, const ColumnView temp,
                       const ColumnView press, const ColumnView qv,
                       const ColumnView ht, const ColumnView pdel, Real pblh)
    : temperature(temp),
      pressure(press),
      vapor_mixing_ratio(qv),
      height(ht),
      hydrostatic_dp(pdel),
      planetary_boundary_height(pblh),
      num_levels_(num_levels) {
  EKAT_REQUIRE_MSG(num_levels > 0,
                   "Number of vertical levels must be positive");
  EKAT_REQUIRE_MSG(pblh >= 0.0,
                   "Planetary boundary height must be non-negative");

  // Make sure the views we're given are properly sized.
  int num_vert_packs = PackInfo::num_packs(num_levels);
  EKAT_REQUIRE_MSG(temp.extent(0) == num_vert_packs,
                   "Temperature view must have extent == " << num_vert_packs);
  EKAT_REQUIRE_MSG(press.extent(0) == num_vert_packs,
                   "Pressure view must have extent == " << num_vert_packs);
  EKAT_REQUIRE_MSG(
      qv.extent(0) == num_vert_packs,
      "Vapor mixing ratio view must have extent == " << num_vert_packs);
  EKAT_REQUIRE_MSG(
      pdel.extent(0) == num_vert_packs,
      "Hydrostatic pressure thickness must have extent == " << num_vert_packs);
  int num_iface_packs = (num_levels_ + 1) / HAERO_PACK_SIZE;
  if (num_iface_packs * HAERO_PACK_SIZE < (num_levels_ + 1)) {
    num_iface_packs++;
  }
  EKAT_REQUIRE_MSG(ht.extent(0) == num_iface_packs,
                   "Height view must have extent == " << num_iface_packs);
}

#if HAERO_FORTRAN

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
extern "C" {

int a_num_levels_c(void* a) {
  auto* atm = static_cast<Atmosphere*>(a);
  return atm->num_levels();
}

void* a_temperature_c(void* a) {
  auto* atm = static_cast<Atmosphere*>(a);
  auto& T = atm->temperature;
  return (void*)T.data();
}

void* a_pressure_c(void* a) {
  auto* atm = static_cast<Atmosphere*>(a);
  auto& p = atm->temperature;
  return (void*)p.data();
}

void* a_vapor_mixing_ratio_c(void* a) {
  auto* atm = static_cast<Atmosphere*>(a);
  auto& qv = atm->vapor_mixing_ratio;
  return (void*)qv.data();
}

void* a_height_c(void* a) {
  auto* atm = static_cast<Atmosphere*>(a);
  auto& h = atm->height;
  return (void*)h.data();
}

void* a_hydrostatic_dp_c(void* a) {
  auto* atm = static_cast<Atmosphere*>(a);
  auto& hdp = atm->hydrostatic_dp;
  return (void*)hdp.data();
}

Real a_pblh_c(void* a) {
  auto* atm = static_cast<Atmosphere*>(a);
  return atm->planetary_boundary_height;
}

}  // extern "C"

#endif  // HAERO_FORTRAN

}  // namespace haero
