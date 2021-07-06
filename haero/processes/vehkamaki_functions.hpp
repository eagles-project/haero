#ifndef HAERO_VEHKAMAKI_FUNCTIONS_HPP
#define HAERO_VEHKAMAKI_FUNCTIONS_HPP

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_math.hpp"
#include "haero/floating_point.hpp"
#include "haero/haero.hpp"
#include "haero/math.hpp"

namespace haero {

/// The functions in this file implement parameterizations described in
/// Vehkamaki et al, An improved parameterization for sulfuric acid–water /
/// nucleation rates for tropospheric and stratospheric conditions,
/// Journal of Geophysical Research 107 (2002).

namespace vehkamaki2002 {

/// Computes the mole fraction of sulfuric acid in a critical cluster as
/// parameterized by Vehkmaki et al (2002).
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [m-3]
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] rel_hum The relative humidity [-]
KOKKOS_INLINE_FUNCTION
PackType h2so4_mole_fraction(const Pack& c_h2so4,
                             const Pack& temp,
                             const Pack& rel_hum) {
  // Convert H2SO4 concentration to [# cm-3].
  auto N_a = c_h2so4 * 1e-6;

  // Calculate the mole fraction using eq 11 of Vehkamaki et al (2002).
  return Pack(0.740997 - 0.00266379 * temp - 0.00349998 * log(N_a) +
    0.0000504022 * temp * log(so4vol) + 0.00201048 * log(rel_hum) -
    0.000183289 * temp * log(rel_hum) + 0.00157407 * square(log(rel_hum)) -
    0.0000179059 * temp * square(log(rel_hum)) +
    0.000184403 * cube(log(rel_hum)) -
    1.50345e-6 * temp * cube(log(rel_hum)));
}

/// Computes the binary nucleation rate [m-3 s-1] as parameterized by
/// Vehkmaki et al (2002).
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [m-3]
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] rel_hum The relative humidity [-]
/// @param [in] x_crit The mole fraction of H2SO4 in a critical cluster [-]
KOKKOS_INLINE_FUNCTION
PackType nucleation_rate(const Pack& c_h2so4,
                         const Pack& temp,
                         const Pack& rel_hum,
                         const Pack& x_crit) {
  // Calculate the coefficients in eq 12 of Vehkamaki et al (2002).
  Pack a = 0.14309 + 2.21956 * temp - 0.0273911 * square(temp) +
    0.0000722811 * cube(temp) + 5.91822 / x_crit;

  Pack b = 0.117489 + 0.462532 * temp - 0.0118059 * square(temp) +
    0.0000404196 * cube(temp) + 15.7963 / x_crit;

  Pack c = -0.215554 - 0.0810269 * temp + 0.00143581 * square(temp) -
    4.7758e-6 * cube(temp) - 2.91297 / x_crit;

  Pack d = -3.58856 + 0.049508 * temp - 0.00021382 * square(temp) +
    3.10801e-7 * cube(temp) - 0.0293333 / x_crit;

  Pack e = 1.14598 - 0.600796 * temp + 0.00864245 * square(temp) -
    0.0000228947 * cube(temp) - 8.44985 / x_crit;

  Pack f = 2.15855 + 0.0808121 * temp - 0.000407382 * square(temp) -
    4.01957e-7 * cube(temp) + 0.721326 / x_crit;

  Pack g = 1.6241 - 0.0160106 * temp + 0.0000377124 * square(temp) +
    3.21794e-8 * cube(temp) - 0.0113255 / x_crit;

  Pack h = 9.71682 - 0.115048 * temp + 0.000157098 * square(temp) +
    4.00914e-7 * cube(temp) + 0.71186 / x_crit;

  Pack i = -1.05611 + 0.00903378 * temp - 0.0000198417 * square(temp) +
    2.46048e-8 * cube(temp) - 0.0579087 / x_crit;

  Pack j = -0.148712 + 0.00283508 * temp - 9.24619e-6 * square(temp) +
    5.00427e-9 * cube(temp) - 0.0127081 / x_crit;

  // Convert H2SO4 concentration to [# cm-3].
  auto N_a = c_h2so4 * 1e-6;

  // Compute the nucleation rate using eq 12.
  return exp(a + b * log(rel_hum) + c * square(log(rel_hum)) +
      d * cube(log(rel_hum)) + ecoe * log(N_a) +
      f * log(rel_hum) * log(N_a) +
      g * square(log(rel_hum)) * (log(N_a)) +
      h * square(log(N_a)) +
      i * log(rel_hum) * square(log(N_a)) +
      j * cube(log(N_a)));
}

/// Computes the total number of molecules in a critical cluster as
/// parameterized in Vehkamaki et al (2002).
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [m-3]
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] rel_hum The relative humidity [-]
/// @param [in] x_crit The mole fraction of H2SO4 in a critical cluster [-]
KOKKOS_INLINE_FUNCTION
PackType num_critical_molecules(const Pack& c_h2so4,
                                const Pack& temp,
                                const Pack& rel_hum,
                                const Pack& x_crit) {
  // Calc the coefficients for the number of molecules in a critical
  // cluster (eq 13).
  Pack A = -0.00295413 - 0.0976834 * temp + 0.00102485 * square(temp) -
    2.18646e-6 * cube(temp) - 0.101717 / x_crit;

  Pack B = -0.00205064 - 0.00758504 * temp + 0.000192654 * square(temp) -
    6.7043e-7 * cube(temp) - 0.255774 / x_crit;

  Pack C = +0.00322308 + 0.000852637 * temp - 0.0000154757 * square(temp) +
    5.66661e-8 * cube(temp) + 0.0338444 / x_crit;

  Pack D = +0.0474323 - 0.000625104 * temp + 2.65066e-6 * square(temp) -
    3.67471e-9 * cube(temp) - 0.000267251 / x_crit;

  Pack E = -0.0125211 + 0.00580655 * temp - 0.000101674 * square(temp) +
    2.88195e-7 * cube(temp) + 0.0942243 / x_crit;

  Pack F = -0.038546 - 0.000672316 * temp + 2.60288e-6 * square(temp) +
    1.19416e-8 * cube(temp) - 0.00851515 / x_crit;

  Pack G = -0.0183749 + 0.000172072 * temp - 3.71766e-7 * square(temp) -
    5.14875e-10 * cube(temp) + 0.00026866 / x_crit;

  Pack H = -0.0619974 + 0.000906958 * temp - 9.11728e-7 * square(temp) -
    5.36796e-9 * cube(temp) - 0.00774234 / x_crit;

  Pack I = +0.0121827 - 0.00010665 * temp + 2.5346e-7 * square(temp) -
    3.63519e-10 * cube(temp) + 0.000610065 / x_crit;

  Pack J = +0.000320184 - 0.0000174762 * temp + 6.06504e-8 * square(temp) -
    1.4177e-11 * cube(temp) + 0.000135751 / x_crit;

  // Convert H2SO4 concentration to [# cm-3].
  auto N_a = c_h2so4 * 1e-6;

  // Compute n_tot using eq 13.
  return exp(A + B * log(rel_hum) + C * square(log(rel_hum)) +
      D * cube(log(rel_hum)) + E * log(N_a) +
      F * log(rel_hum) * log(N_a) +
      G * square(log(rel_hum)) * log(N_a) +
      H * square(log(N_a)) +
      I * log(rel_hum) * square(log(N_a)) +
      J * cube(log(N_a)));
}

/// Computes the radius [m] of a critical cluster as parameterized in Vehkamaki
/// et al (2002).
/// @param [in] x_crit The mole fraction of H2SO4 in a critical cluster [-]
/// @param [in] n_tot The total number of molecules in the critical cluster [-]
KOKKOS_INLINE_FUNCTION
PackType critical_radius(const Pack& x_crit, const Pack& n_tot) {
  return 1e-9 * exp(-1.6524245 + 0.42316402 * x_crit + 0.3346648 * log(n_tot));
}

} // namespace vehkamaki2002

} // namespace haero
#endif
