#include "haero/atmosphere.hpp"

namespace haero {

Atmosphere::Atmosphere(int num_levels,
                       const ColumnView temp,
                       const ColumnView press,
                       const ColumnView rel_hum,
                       const ColumnView ht,
                       const ColumnView pdel,
                       Real pbl):
  temperature(temp),
  pressure(press),
  relative_humidity(rel_hum),
  height(ht),
  hydrostatic_dp(pdel),
  planetary_boundary_height(pbl),
  num_levels_(num_levels) {
  EKAT_REQUIRE_MSG(num_levels > 0,
                   "Number of vertical levels must be positive");
  EKAT_REQUIRE_MSG(pbl >= 0.0,
                   "Planetary boundary height must be non-negative");

  // Make sure the views we're given are properly sized.
  int num_vert_packs = num_levels_/HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels_) {
    num_vert_packs++;
  }
  EKAT_REQUIRE_MSG(temp.extent(0) == num_vert_packs,
                   "Temperature view must have extent == " << num_vert_packs);
  EKAT_REQUIRE_MSG(press.extent(0) == num_vert_packs,
                   "Pressure view must have extent == " << num_vert_packs);
  EKAT_REQUIRE_MSG(rel_hum.extent(0) == num_vert_packs,
                   "Relative humidity view must have extent == " << num_vert_packs);
  EKAT_REQUIRE_MSG(pdel.extent(0) == num_vert_packs,
                   "Hydrostatic pressure thickness must have extent == " << num_vert_packs);
  int num_iface_packs = (num_levels_+1)/HAERO_PACK_SIZE;
  if (num_iface_packs * HAERO_PACK_SIZE < (num_levels_+1)) {
    num_iface_packs++;
  }
  EKAT_REQUIRE_MSG(ht.extent(0) == num_iface_packs,
                   "Height view must have extent == " << num_iface_packs);
}

Atmosphere::~Atmosphere() {
}

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
extern "C" {

void* a_temperature_c(void* a)
{
  auto* atm = static_cast<Atmosphere*>(a);
  auto& T = atm->temperature;
  return (void*)T.data();
}

void* a_pressure_c(void* a)
{
  auto* atm = static_cast<Atmosphere*>(a);
  auto& p = atm->temperature;
  return (void*)p.data();
}

void* a_relative_humidity_c(void* a)
{
  auto* atm = static_cast<Atmosphere*>(a);
  auto& RH = atm->relative_humidity;
  return (void*)RH.data();
}

void* a_height_c(void* a)
{
  auto* atm = static_cast<Atmosphere*>(a);
  auto& h = atm->height;
  return (void*)h.data();
}

void* a_hydrostatic_dp_c(void* a) {
  auto* atm = static_cast<Atmosphere*>(a);
  auto& hdp = atm->hydrostatic_dp;
  return (void*)hdp.data();
}

Real a_pblh_c(void* a)
{
  auto* atm = static_cast<Atmosphere*>(a);
  return atm->planetary_boundary_height;
}

} // extern "C"

} // haero

