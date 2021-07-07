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
/// 112 (2007). Also included are corrections described in Merikanto et al,
/// Correction to "New parameterization...", Journal of Geophysical Research
/// 114 (2009).

/// These parameterizations are valid for the following ranges:
/// temperature:                235 - 295 K
/// relative humidity:          0.05 - 0.95
/// H2SO4 number concentration: 5e4 - 1e9 cm-3
/// NH3 molar mixing ratio:     0.1 - 1000 ppt

namespace merikanto2007 {

/// Computes the logarithm of the ternary nucleation rate [cm-3 s-1] as
/// parameterized by Merikanto et al (2007), eq 8.
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] rel_hum The relative humidity [-]
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
/// @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
KOKKOS_INLINE_FUNCTION
PackType log_nucleation_rate(const PackType& temp, const PackType& rel_hum,
                             const PackType& c_h2so4, const PackType& xi_nh3) {
  // Compute the 20 polynomial fit functions f_i(T).

  // Polynomial coefficients (Table 1). According to the correction (Merikanto
  // et al, 2009),
  // double precision is required.
  static const double a[21][4] = {
      {0, 0, 0, 0},  // <-- padding so we can use one-based indexing
      {-358.233770505299, 4.8630382337426985, -0.02175548069741675,
       0.00003212869941055865},
      {-980.923146020468, 10.054155220444462, -0.03306644502023841,
       0.000034274041225891804},
      {1200.472096232311, -17.37107890065621, 0.08170681335921742,
       -0.00012534476159729881},
      {-14.833042158178936, 0.2932631303555295, -0.0016497524241142845,
       2.844074805239367e-6},
      {-4.39129415725234e6, 56383.93843154586, -239.835990963361,
       0.33765136625580167},
      {4.905527742256349, -0.05463019231872484, 0.00020258394697064567,
       -2.502406532869512 - 7},
      {-231375.56676032578, 2919.2852552424706, -12.286497122264588,
       0.017249301826661612},
      {75061.15281456841, -931.8802278173565, 3.863266220840964,
       -0.005349472062284983},
      {-3180.5610833308, 39.08268568672095, -0.16048521066690752,
       0.00022031380023793877},
      {-100.21645273730675, 0.977886555834732, -0.0030511783284506377,
       2.967320346100855e-6},
      {5599.912337254629, -70.70896612937771, 0.2978801613269466,
       -0.00041866525019504},
      {2.360931724951942e6, -29752.130254319443, 125.04965118142027,
       -0.1752996881934318},
      {16597.75554295064, -175.2365504237746, 0.6033215603167458,
       -0.0006731787599587544},
      {-89.38961120336789, 1.153344219304926, -0.00495454549700267233,
       7.096309866238719e-6},
      {-629.7882041830943, 7.772806552631709, -0.031974053936299256,
       0.00004383764128775082},
      {-732006.8180571689, 9100.06398573816, -37.771091915932004,
       0.05235455395566905},
      {40751.075322248245, -501.66977622013934, 2.063469732254135,
       -0.002836873785758324},
      {-1911.0303773001353, 23.6903969622286, -0.09807872005428583,
       0.00013564560238552576},
      {2.792313345723013, -0.03422552111802899, 0.00014019195277521142,
       -1.9201227328396297e-7},
      {3.1712136610383244, -0.037822330602328806, 0.0001500555743561457,
       -1.9828365865570703e-7}};
  PackType f[21];
  for (int i = 0; i <= 20; ++i) {
    f[i] = PackType(a[i][0] + a[i][1] * temp + a[i][2] * square(temp) +
                    a[i][3] * cube(temp));
  }

  // Compute the logarithm of the nucleation rate [cm-3 s-1] using eq 8.
  auto c = c_h2so4;
  auto xi = xi_nh3;
  return -12.86185 + f[1] * rel_hum + f[2] * log(rel_hum) + f[3] * log(c) +
         f[4] * square(log(c)) + f[5] / square(log(c)) + f[6] * xi +
         f[7] * log(xi) + f[8] * square(log(xi)) + f[9] * cube(log(xi)) +
         f[10] * rel_hum * log(xi) + f[11] * log(c) * log(xi) +
         f[12] * log(xi) / log(c) + f[13] * log(rel_hum) / log(c) +
         f[14] * log(rel_hum) * log(xi) +
         f[15] * rel_hum / (cube(xi) * log(c)) +
         f[16] * square(log(xi)) / log(c) + f[17] * cube(log(xi)) / log(c) +
         f[18] * log(c) * square(log(xi)) +
         f[19] * square(log(c)) * cube(log(xi)) +
         f[20] * log(rel_hum) * cube(log(xi));
}

/// Computes the "onset temperature" [K] (eq 10) above which Merikanto's
/// parameterization for the nucleation rate (eq 8) cannot be used (in which
/// case the authors suggest setting the nucleation rate to zero).
/// @param [in] rel_hum The relative humidity [-]
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
/// @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
KOKKOS_INLINE_FUNCTION
PackType onset_temperature(const PackType& rel_hum, const PackType& c_h2so4,
                           const PackType& xi_nh3) {
  return PackType(143.600 + 1.01789 * rel_hum + 10.1964 * log(c_h2so4) -
                  0.184988 * square(log(c_h2so4)) - 17.1618 * log(xi_nh3) +
                  109.9247 * log(xi_nh3) / log(c_h2so4) +
                  0.773412 * log(c_h2so4) * log(xi_nh3) -
                  0.155764 * square(log(xi_nh3)));
}

/// Computes the radius of a critical cluster [nm] as parameterized in Merikanto
/// et al (2007), eq 11.
/// @param [in] log_J The logarithm of the nucleation rate ["log (cm-3 s-1)"]
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
/// @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
KOKKOS_INLINE_FUNCTION
PackType critical_radius(const PackType& log_J, const PackType& temp,
                         const PackType& c_h2so4, const PackType& xi_nh3) {
  auto c = c_h2so4;
  auto xi = xi_nh3;
  return PackType(
      0.328886 - 0.00337417 * temp + 0.0000183474 * square(temp) +
      0.00254198 * log(c) - 0.0000949811 * temp * log(c) +
      0.000744627 * square(log(c)) + 0.0243034 * log(xi) +
      0.0000158932 * temp * log(xi) - 0.00203460 * log(c) * log(xi) -
      0.000559304 * square(log(xi)) - 4.88951e-7 * temp * square(log(xi)) +
      0.000138470 * cube(log(xi)) + 4.14108e-6 * log_J -
      0.0000268131 * temp * log_J + 0.00128791 * log(xi) * log_J -
      3.80352e-6 * temp * log(xi) * log_J - 0.0000187902 * square(log_J));
}

/// Computes the total number of molecules in a critical cluster as
/// parameterized in Merikanto et al (2007), eq 12.
/// @param [in] log_J The logarithm of the nucleation rate ["log (cm-3 s-1)"]
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
/// @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
KOKKOS_INLINE_FUNCTION
PackType num_critical_molecules(const PackType& log_J, const PackType& temp,
                                const PackType& c_h2so4,
                                const PackType& xi_nh3) {
  auto c = c_h2so4;
  auto xi = xi_nh3;
  return PackType(
      57.4009 - 0.299634 * temp + 0.000739548 * square(temp) -
      5.09060 * log(c) + 0.0110166 * temp * log(c) +
      0.0675003 * square(log(c)) - 0.810283 * log(xi) +
      0.0159051 * temp * log(xi) - 0.204417 * log(c) * log(xi) +
      0.0891816 * square(log(xi)) - 0.000496903 * temp * square(log(xi)) +
      0.00570439 * cube(log(xi)) + 3.40987 * log_J - 0.0149170 * temp * log_J +
      0.0845909 * log(xi) * log_J - 0.000148006 * temp * log(xi) * log_J +
      0.00503805 * square(log_J));
}

/// Computes the total number of H2SO4 molecules in a critical cluster as
/// parameterized in Merikanto et al (2007), eq 13.
/// @param [in] log_J The logarithm of the nucleation rate ["log (cm-3 s-1)"]
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
/// @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
KOKKOS_INLINE_FUNCTION
PackType num_h2so4_molecules(const PackType& log_J, const PackType& temp,
                             const PackType& c_h2so4, const PackType& xi_nh3) {
  auto c = c_h2so4;
  auto xi = xi_nh3;
  return PackType(
      -4.71542 + 0.134364 * temp - 0.000471847 * square(temp) -
      2.56401 * log(c) + 0.0113533 * temp * log(c) +
      0.00108019 * square(log(c)) + 0.517137 * log(xi) -
      0.00278825 * temp * log(xi) + 0.806697 * square(log(xi)) -
      0.00318491 * temp * square(log(xi)) - 0.0995118 * cube(log(xi)) +
      0.000400728 * temp * cube(log(xi)) + 1.32765 * log_J -
      0.00616765 * temp * log_J - 0.110614 * log(xi) * log_J +
      0.000436758 * temp * log(xi) * log_J + 0.000916366 * square(log_J));
}

/// Computes the total number of NH3 molecules in a critical cluster as
/// parameterized in Merikanto et al (2007), eq 14.
/// @param [in] log_J The logarithm of the nucleation rate ["log (cm-3 s-1)"]
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
/// @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
KOKKOS_INLINE_FUNCTION
PackType num_nh3_molecules(const PackType& log_J, const PackType& temp,
                           const PackType& c_h2so4, const PackType& xi_nh3) {
  auto c = c_h2so4;
  auto xi = xi_nh3;
  return PackType(
      71.2007 - 0.840960 * temp + 0.00248030 * square(temp) + 2.77986 * log(c) -
      0.0147502 * temp * log(c) + 0.0122645 * square(log(c)) -
      2.00993 * log(xi) + 0.00868912 * temp * log(xi) -
      0.00914118 * log(c) * log(xi) + 0.137412 * square(log(xi)) -
      0.000625323 * temp * square(log(xi)) + 0.0000937733 * cube(log(xi)) +
      0.520297 * log_J - 0.00241987 * temp * log_J +
      0.0791639 * log(xi) * log_J - 0.000302159 * temp * log(xi) * log_J +
      0.00469770 * square(log_J));
}

/// Computes the logarithm of the threshold number concentration of H2SO4 [cm-3]
/// that produces a nucleation rate of 1 cm-3 s-1 at the given temperature and
/// molar mixing ratio of NH3.
/// @param [in] temp The atmospherіc temperature [K]
/// @param [in] xi_nh3 The molar mixing ratio of NH3 [ppt]
KOKKOS_INLINE_FUNCTION
PackType log_h2so4_nucleation_threshold(const PackType& temp,
                                        const PackType& xi_nh3) {
  auto xi = xi_nh3;
  return PackType(
      -40.5988 + 5.00845 / xi + 0.00995956 * xi + 0.231207 * temp -
      0.0191883 * temp / xi - 0.0000312301 * temp * xi + 15.4213 * log(xi) -
      0.06366755 * temp * log(xi) - 3.48925 * square(log(xi)) +
      0.0143679 * temp * square(log(xi)) + 0.234708 * cube(log(xi)) -
      0.000995330 * temp * cube(log(xi)));
}

}  // namespace merikanto2007

}  // namespace haero
#endif
