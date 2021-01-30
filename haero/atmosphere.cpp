#include "haero/atmosphere.hpp"

namespace haero {

Atmosphere::Atmosphere(int num_levels,
                       const Kokkos::View<PackType*>& temp,
                       const Kokkos::View<PackType*>& press,
                       const Kokkos::View<PackType*>& rel_hum,
                       const Kokkos::View<PackType*>& ht):
  num_levels_(num_levels),
  temperature_(temp),
  pressure_(press),
  relative_humidity_(rel_hum),
  height_(ht) {
  // Make sure the views we're given are properly sized.
  EKAT_REQUIRE_MSG(temp.extent(0) == num_levels/HAERO_PACK_SIZE,
                   "Temperature view has wrong number of vertical level packs!");
  EKAT_REQUIRE_MSG(press.extent(0) == num_levels/HAERO_PACK_SIZE,
                   "Pressure view has wrong number of vertical level packs!");
  EKAT_REQUIRE_MSG(rel_hum.extent(0) == num_levels/HAERO_PACK_SIZE,
                   "Relative humidity view has wrong number of vertical level packs!");
  EKAT_REQUIRE_MSG(ht.extent(0) == (num_levels+1)/HAERO_PACK_SIZE,
                   "Height view has wrong number of vertical interface packs!");
}

Atmosphere::~Atmosphere() {
}

// Interoperable C functions for providing data to Fortran.
// See haero.F90 for details on how these functions are used.
extern "C" {

void* a_temperature_c(void* a)
{
  auto* atm = static_cast<Atmosphere*>(a);
  auto& T = atm->temperature();
  return (void*)T.data();
}

void* a_pressure_c(void* a)
{
  auto* atm = static_cast<Atmosphere*>(a);
  auto& p = atm->temperature();
  return (void*)p.data();
}

void* a_relative_humidity_c(void* a)
{
  auto* atm = static_cast<Atmosphere*>(a);
  auto& RH = atm->relative_humidity();
  return (void*)RH.data();
}

void* a_height_c(void* a)
{
  auto* atm = static_cast<Atmosphere*>(a);
  auto& h = atm->height();
  return (void*)h.data();
}

} // extern "C"

} // haero

