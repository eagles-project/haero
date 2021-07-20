#ifndef HAERO_KERMINEN2002_HPP
#define HAERO_KERMINEN2002_HPP

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_math.hpp"
#include "haero/haero.hpp"
#include "haero/math.hpp"
#include "haero/physical_constants.hpp"

namespace haero {

/// The functions in this file implement parameterizations described in
/// Kerminen and Kulmala, Analytical formulae connecting the "real" and the
/// "apparent" nucleation rate and the nuclei number concentration for
/// atmospheric nucleation events, Aerosol Science 33 (2002) pp. 609-622

namespace kerminen2002 {

/// Computes a factor that, when multiplied by a "real" nucleation rate, yields
/// the "apparent" rate of nucleation for the given molecule size. This
/// formulation assumes a single nucleating/condensing gas species.
/// @param [in] aero_num_conc The number concentration for the aerosol species
///                           to which nuclei are added after growth [#/cc]
/// @param [in] mean_aero_diam The mean diameter of the pre-existing aerosol
///                            population to which growing nuclei are added [nm]
/// @param [in] gas_accom_coeff The accommodation coefficient for the gas [-]
/// @param [in] gas_temperature The temperature of the nucleating gas [K]
/// @param [in] nuc_growth_rate The nuclei growth rate [nm / h]
/// @param [in] nuc_mass_density The mass density of nuclei [kg/m3]
/// @param [in] initial_nuc_diam The initial diameter of a fresh nucleus [nm]
/// @param [in] final_nuc_diam The final size of a condensed nucleus [nm]
KOKKOS_INLINE_FUNCTION
PackType apparent_nucleation_factor(
    const PackType& aero_num_conc, Real mean_aero_diam, Real gas_accom_coeff,
    const PackType& gas_temperature, const PackType& nuc_growth_rate,
    Real nuc_mass_density, Real initial_nuc_diam, Real final_nuc_diam) {
  EKAT_KERNEL_REQUIRE(initial_nuc_diam < final_nuc_diam);

  // Compute beta_m, the transitional correction for the condensational mass
  // flux (Fuchs and Sutugin (1971) via Kerminen et al 2002, eq 4).
  PackType lambda =
      gas_kinetics::mean_free_path(PackType(initial_nuc_diam), aero_num_conc);
  PackType Kn =
      gas_kinetics::knudsen_number(lambda, PackType(0.5 * final_nuc_diam));
  PackType beta_m =
      (1 + Kn) / (1 + 0.377 * Kn + (1.33 * Kn * (1 + Kn) / gas_accom_coeff));

  // Compute the condensation sink CS1 [m-2] using Kerminen et al 2002, eq 3.
  // The conversion factor 1e-3 is the product of two conversions:
  // 1. 1e-9 for final_nuc_diam -> [m]
  // 2. 1e6 for aero_num_conc -> [#/m3].
  PackType CS1 = 0.5 * 1e-3 * final_nuc_diam * beta_m * aero_num_conc;

  // Compute gamma, the proportionality factor for the exponential factor eta.
  Real gamma0 = 0.23;        // [nm2 m2 / h]
  PackType gamma = gamma0 *  // Kerminen et al 2002, eq 22
                   pow(initial_nuc_diam, 0.2) * pow(final_nuc_diam / 3, 0.075) *
                   pow(mean_aero_diam / 150, 0.048) *
                   pow(nuc_mass_density / 1000, -0.33) *
                   pow(gas_temperature / 293, -0.75);
  PackType eta = gamma * CS1 / nuc_growth_rate;  // Kerminen et al 2002, eq 11
  return exp(eta / final_nuc_diam - eta / initial_nuc_diam);
}

}  // namespace kerminen2002

}  // namespace haero
#endif
