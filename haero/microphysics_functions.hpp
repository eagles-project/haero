#ifndef HAERO_MICROPHYSICS_FUNCTIONS_HPP
#define HAERO_MICROPHYSICS_FUNCTIONS_HPP

#include "haero/haero.hpp"
#include "constants.hpp"

namespace haero {

/** This file contains functions that are used by multiple processes and/or
diagnostics. Its purpose is like that of the constants.hpp file, except that
these are functions, not constant values.

Every function in this file should be analogous to "pure" functions in fortran,
and capable of running on the device.

*/


/** @brief This function defines the "A" constant that describes the effect
 of surface tension and curvature on a water droplet, as
 described by Kelvin's equation.  Frequently,
 this quantity is denoted "A" in the literature.

@warning It may or may not be used in SI units.

Pruppacher/Klett 2nd ed. (6.28) defines
 @f$  A = \frac{2 M_{H_2O} \sigma_{w/a}}{\mathcal{R} T \rho_{H_2O}} @f$,
 where the molar mass of water, M_w, surface tension \sigma, the
 universal gas constant R, and the density of liquid water \rho_w,
 are all defined in c.g.s. units.   Temperature T is given in Kelvin.
 A has units of length, cm, in c.g.s. units.

In legacy MAM4, it is computed as
 - modal_aero_wateruptake.F90: a = 2.e4_r8*mw*surften/(ugascon*tair*rhow)
      with mw = 18 kg/mol, surften = 76 erg/cm2, ugascon = 8.3e7 erg/K/mol,
      tair = 273 K, and rhow = 1 g/cm3
      so that a = 0.001207 = 1.121e-3.  We conclude that changing 2 to 2e4
      is the conversion of centimeters to microns.
 - ndrop.F90: aten = 2._r8*mwh2o*surften/(r_universal*t0*rhoh2o)
      with mwh20 = 18.016 kg/kmol, surften = 0.076 N/m,
      r_universal = 8.31456e3 J/K/kmol, t0 = 273 K, rhoh2o 1000 kg/m3,
      so that aten = 1.206e-9 m.
      ndrop.F90 uses kmol, not mol, so is not quite the same as SI.

  Prupaccher & Klett 2nd ed. eqn. (6.28) approximates A in cgs units
    as 3.3e-5/T with T in K.  For T = 273.15, this gives 1.208e-7 cm.

In Haero, consistent with our SI unit design convention, we use constants.hpp:
  molar mass of water =

  @param T temperature [K]
  @return A [m]
*/
template <typename ScalarType> KOKKOS_INLINE_FUNCTION
ScalarType kelvin_coeff_A(const ScalarType& T) {
  const Real molar_mass_water = Constants::molec_weight_h2o;
  const Real rgas = Constants::r_gas;
  const Real rho_h2o = Constants::density_h2o;
  return 2*molar_mass_water * surface_tension_water_air(T) /
    (rgas * T * rho_h2o);
}


/** @brief Surface tension of water/air interface

  Legacy MAM4 uses T = 273 everywhere surface tension is
  required.

  Defined by Prupaccher/Klett 2nd. ed. eqn (5.12).

  @param T temperature [K]
  @return surface tension [N/m]
*/
template <typename ScalarType> KOKKOS_INLINE_FUNCTION
ScalarType surface_tension_water_air(const ScalarType& T) {
  const Real coeffs[7] = {75.93,    // a0
                          0.115,    // a1
                          6.818e-2, // a2
                          6.511e-3, // a3
                          2.933e-4, // a4
                          6.283e-6, // a5
                          5.285e-8} // a6
  const Real K_to_C = Constants::freezing_pt_h2o;
  const Real erg_per_cm2_to_N_per_m = 1e-3;
  ScalarType result = 0;
  for (int i=0; i<7; ++i) {
    result += coeffs[i] * pow(T - K_to_C, i);
  }
  return erg_per_cm2_to_N_per_m * result;
}

} // namespace haero
#endif
