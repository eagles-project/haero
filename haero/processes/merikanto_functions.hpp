#ifndef HAERO_MERIKANTO_FUNCTIONS_HPP
#define HAERO_MERIKANTO_FUNCTIONS_HPP

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_math.hpp"
#include "haero/floating_point.hpp"
#include "haero/haero.hpp"
#include "haero/math.hpp"

namespace haero {

/// The functions in this file implement parameterizations described in
/// Merikanto et al, New parameterization of sulfuric acid-ammonia-water ternary
/// nucleation rates at tropospheric conditions, Journal of Geophysical Research
/// 112 (2007).

namespace merikanto2007 {

/// Computes the logarithm of the ternary nucleation rate [cm-3 s-1] as
/// parameterized by Merikanto et al (2007), eq 8. NOTE THE DIFFERENT UNITS in
/// log J as computed by this function!
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] rel_hum The relative humidity [-]
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [m-3]
/// @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
KOKKOS_INLINE_FUNCTION
PackType log_nucleation_rate(const Pack& temp,
                             const Pack& rel_hum,
                             const Pack& c_h2so4,
                             const Pack& xi_nh3) {
  // Compute the 20 polynomial fit functions f_i.

  // Polynomial coefficients (Table 1).
  static const Real a[21][4] = {
    {0, 0, 0, 0}, // <-- padding so we can use one-based indexing
    {-358.234, 4.86304, -0.0217555, 0.0000321287},
    {-980.923, 10.0542, -0.0330664, 0.0000342740},
    {1200.47, -17.3711, 0.0817068, -0.000125345},
    {-14.8330, 0.293263, -0.00164975, 2.84407e-6},
    {-4.39129e6, 56383.9, -239.836, 0.337651},
    {4.90553, -0.0546302, 0.000202584, -2.50241e-7},
    {-231376, 2919.29, -12.2865, 0.0172493},
    {75061.2, -931.880, 3.86327, -0.00534947},
    {-3180.56, 39.0827, -0.160485, 0.000220314},
    {-100.216, 0.977887, -0.00305118, 2.96732e-6},
    {5599.91, -70.7090, 0.297880, -0.000418665},
    {2.36093e6, -29752.1, 125.050, -0.175300},
    {16597.8, -175.237, 0.603322, -0.000673179},
    {-89.3896, 1.15334, -0.00495455, 7.09631e-6},
    {-629.788, 7.77281, -0.0319741, 0.0000438376},
    {-732007, 9100.06, -37.7711, 0.0523546},
    {40751.1, -501.670, 2.06347, -0.00283687},
    {-1911.03, 23.6904, -0.0980787, 0.000135646},
    {2.79231, -0.0342255, 0.000140192, -1.92012e-7},
    {3.17121, -0.0378223, 0.000150056, -1.98284e-7}
  };
  Pack f[21];
  for (int i = 0; i <= 20; ++i) {
    f[i] = Pack(a[i][0] + a[i][1]*temp + a[i][2]*square(temp) + a[i][3]*cube(temp));
  }

  // Compute the logarithm of the nucleation rate [cm-3 s-1] using eq 8.
  auto c = 1e-6 * c_h2so4; // H2SO4 concentration in [cm-3]
  auto xi = xi_nh3;
  return -12.86185 + f[1]*rel_hum + f[2]*log(rel_hum) +
         f[3]*log(c) + f[4]*square(log(c)) + f[5]/square(log(c)) +
         f[6]*xi + f[7]*log(xi) + f[8]*square(log(xi)) + f[9]*cube(log(xi)) +
         f[10]*rel_hum*log(xi) + f[11]*log(c)*log(xi) + f[12]*log(xi)/log(c) +
         f[13]*log(rel_hum)/log(c) + f[14]*log(rel_hum)*log(xi) +
         f[15]/(cube(xi)*log(c)) + f[16]*square(log(xi))/log(c) +
         f[17]*cube(log(xi))/log(c) + f[18]*log(c)*square(log(xi)) +
         f[19]*square(log(c))*cube(log(xi)) + f[20]*log(rel_hum)*cube(log(xi));
}

/// Computes the "onset temperature" [K] (eq 10) below which Merikanto's
/// parameterization for the nucleation rate (eq 8) cannot be used (in which
/// case the authors suggest setting the nucleation rate to zero).
/// @param [in] rel_hum The relative humidity [-]
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [m-3]
/// @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
KOKKOS_INLINE_FUNCTION
PackType onset_temp(const Pack& rel_hum,
                    const Pack& c_h2so4,
                    const Pack& xi_nh3) {
  return Pack(143.6002929064716 + 1.0178856665693992 * rel_hum +
              10.196398812974294 * log(c_h2so4) -
              0.1849879416839113 * square(log(c_h2so4)) -
              17.161783213150173 * log(xi_nh3) +
              (109.92469248546053 * log(xi_nh3)) / log(c_h2so4) +
              0.7734119613144357 * log(c_h2so4) * log(xi_nh3) -
              0.15576469879527022 * square(log(xi_nh3)));
}

/// Computes the radius of a critical cluster [m] as parameterized in Merikanto
/// et al (2007), eq 11
/// @param [in] log_J The logarithm of the nucleation rate ["log (cm-3 s-1)"]
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [# m-3]
/// @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
KOKKOS_INLINE_FUNCTION
PackType critical_radius(const Pack& log_J,
                         const Pack& temp,
                         const Pack& c_h2so4,
                         const Pack& xi_nh3) {
  auto c = 1e-6 * c_h2so4; // H2SO4 concentration in [cm-3]
  auto xi = xi_nh3;
  return 1e-9 *
    Pack(0.328886 - 0.00337417*temp + 0.0000183474 * square(temp)
       + 0.00254198*log(c) - 0.0000949811*temp*log(c)
       + 0.000744627*square(log(c)) + 0.0243034*log(xi)
       + 0.0000158932*temp*log(xi) - 0.00203460*log(c)*log(xi)
       - 0.000559304*square(log(xi)) - 4.88951e-7*temp*square(log(xi))
       + 0.000138470*cube(log(xi)) + 4.14108e-6*log_J - 0.0000268131*temp*log_J
       + 0.00128791*log(xi)*log_J - 3.80352e-6*temp*log(xi)*log_J
       - 0.0000187902*square(log_J));
}

/// Computes the total number of molecules in a critical cluster as
/// parameterized in Merikanto et al (2007), eq 12
/// @param [in] log_J The logarithm of the nucleation rate ["log (cm-3 s-1)"]
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [# m-3]
/// @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
KOKKOS_INLINE_FUNCTION
PackType num_critical_molecules(const Pack& log_J,
                                const Pack& temp,
                                const Pack& c_h2so4,
                                const Pack& xi_nh3) {
  auto c = 1e-6 * c_h2so4; // H2SO4 concentration in [cm-3]
  auto xi = xi_nh3;
  return Pack(57.4009 - 0.299634*temp + 0.000739548*square(temp)
            - 5.09060*log(c) + 0.0110166*temp*log(c) + 0.0675003*square(log(c))
            - 0.810283*log(xi) + 0.0159051*temp*log(xi)
            - 0.204417*log(c)*log(xi) + 0.0891816*square(log(xi))
            - 0.000496903*temp*square(log(xi)) + 0.00570439*cube(log(xi))
            + 3.40987*log_J - 0.0149170*temp*log_J + 0.0845909*log(xi)*log_J
            - 0.000148006*temp*log(xi)*log_J + 0.00503805*square(log_J));
}

/// Computes the total number of H2SO4 molecules in a critical cluster as
/// parameterized in Merikanto et al (2007).
/// @param [in] log_J The logarithm of the nucleation rate ["log (cm-3 s-1)"]
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [# m-3]
/// @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
KOKKOS_INLINE_FUNCTION
PackType num_h2so4_molecules(const Pack& log_J,
                             const Pack& temp,
                             const Pack& c_h2so4,
                             const Pack& xi_nh3) {
  auto c = 1e-6 * c_h2so4; // H2SO4 concentration in [cm-3]
  auto xi = xi_nh3;
  return Pack(-4.71542 + 0.134364*temp - 0.000471847*square(temp)
            - 2.56401*log(c) + 0.0113533*temp*log(c) + 0.00108019*square(log(c))
            + 0.517137*log(xi) - 0.00278825*temp*log(xi)
            + 0.806697*square(log(xi)) - 0.00318491*temp*square(log(xi))
            - 0.0995118*cube(log(xi)) + 0.000400728*temp*cube(log(xi))
            + 1.32765*log_J - 0.00616765*temp*log_J - 0.110614*log(xi)*log_J
            + 0.000436758*temp*log(xi)*log_J + 0.000916366*square(log_J));
}

/// Computes the total number of NH3 molecules in a critical cluster as
/// parameterized in Merikanto et al (2007).
/// @param [in] log_J The logarithm of the nucleation rate ["log (cm-3 s-1)"]
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [# m-3]
/// @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
KOKKOS_INLINE_FUNCTION
PackType num_nh3_molecules(const Pack& log_J,
                           const Pack& temp,
                           const Pack& c_h2so4,
                           const Pack& xi_nh3) {
  auto c = 1e-6 * c_h2so4; // H2SO4 concentration in [cm-3]
  auto xi = xi_nh3;
  return Pack(71.2007 - 0.840960*temp + 0.00248030*square(temp)
            + 2.77986*log(c) - 0.0147502*temp*log(c) + 0.0122645*square(log(c))
            - 2.00993*log(xi) + 0.00868912*temp*log(xi)
            - 0.00914118*log(c)*log(xi) + 0.137412*square(log(xi))
            - 0.000625323*temp*square(log(xi)) + 0.0000937733*cube(log(xi))
            + 0.520297*log_J - 0.00241987*temp*log_J
            + 0.0791639*log(xi)*log_J - 0.000302159*temp*log(xi)*log_J
            + 0.00469770*square(log_J));
}

} // namespace merikanto2007

}  // namespace haero
#endif
