#include "haero/atmosphere.hpp"

namespace haero {

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

