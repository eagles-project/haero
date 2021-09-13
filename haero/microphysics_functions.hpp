#ifndef HAERO_MICROPHYSICS_FUNCTIONS_HPP
#define HAERO_MICROPHYSICS_FUNCTIONS_HPP

#include "haero/haero.hpp"
#include "constants.hpp"
#include "math.hpp"

namespace haero {

/** This file contains functions that are used by multiple processes and/or
diagnostics. Its purpose is like that of the constants.hpp file, except that
these are functions, not constant values.

Every function in this file should be analogous to "pure" functions in fortran,
and capable of running on the device.

*/


/** @brief Surface tension of water/air interface

  Legacy MAM4 uses T = 273 everywhere surface tension is
  required.

  Primary temperature range is defined by the quadratic best fit coefficients
  for Drop Volume data in Table IV from Gittens, 1969, Variation of
  surface tension of water with temperature, J. Colloid and Interface Sciences
  30(3), 406--412.

  Extrapolation to supercooled regime given by
  Defined by Prupaccher/Klett 2nd. ed. eqn (5.12),
  @f$ \sigma_{w/a} = \sum_{i=0}^6 a_0 T^i @f$
  where T is deg C and the coefficients are a curve fit

  @param T temperature [K]
  @return surface tension [N/m]
*/
template <typename ScalarType> KOKKOS_INLINE_FUNCTION
ScalarType surface_tension_water_air(const ScalarType& T) {
  const Real coeffs_extrap[7] = {75.93,     // a0
                          0.115,     // a1
                          6.818e-2,  // a2
                          6.511e-3,  // a3
                          2.933e-4,  // a4
                          6.283e-6,  // a5
                          5.285e-8}; // a6
  const Real coeffs_interp[3] = {75.93, -0.1365, -0.3827e-3};
  const Real K_to_C = Constants::freezing_pt_h2o;
  const Real erg_per_cm2_to_N_per_m = 1e-3;

  const bool do_extrap = T < Constants::freezing_pt_h2o;

  ScalarType result;
  if (do_extrap) {
    // the curve fit is defined in degrees C, so we have to watch out for
    // 0^0.  Easy fix: move the 0th term of the sum outside the loop.
    result = coeffs_extrap[0];
    for (int i=1; i<7; ++i) {
      result += coeffs_extrap[i] * pow(T - K_to_C, i);
    }
  }
  else {
    result = coeffs_interp[0] + coeffs_interp[1]*(T - K_to_C) +
      coeffs_interp[2]*square(T-K_to_C);
  }
  return erg_per_cm2_to_N_per_m * result;
}


/** @brief This function defines the "A" constant that describes the effect
 of surface tension and curvature on a water droplet, as
 described by Kelvin's equation.  A has units of length. Frequently,
 this quantity is denoted "A" in the literature.

Pruppacher/Klett 2nd ed. (6.28) defines
 @f$  A = \frac{2 M_{H_2O} \sigma_{w/a}}{\mathcal{R} T \rho_{H_2O}} @f$,
 where the molar mass of water, M_w, surface tension \sigma, the
 universal gas constant R, and the density of liquid water \rho_w,
 are all defined in c.g.s. units.   Temperature T is given in Kelvin.
 A is given in cm in c.g.s. units.

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
      ndrop.F90 uses kmol, not mol, so is not quite the same as SI, but its
      output is in meters as is ours, below.

  Prupaccher & Klett 2nd ed. eqn. (6.28) approximate A in cgs units
    as 3.3e-5/T with T in K.  For T = 273.15, this gives 1.208e-7 cm.
    This ignores the variance of surface tension with temperature.

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




} // namespace haero
#endif
