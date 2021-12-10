#ifndef HAERO_KERMINEN2002_HPP
#define HAERO_KERMINEN2002_HPP

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_math.hpp"
#include "haero/haero.hpp"
#include "haero/math.hpp"

namespace haero {

/// The functions in this file implement parameterizations described in
/// Kerminen and Kulmala, Analytical formulae connecting the “real” and the
/// “apparent” nucleation rate and the nuclei number concentration for
/// atmospheric nucleation events, Aerosol Science 33 (2002).

namespace kerminen2002 {

/// Computes a factor that transforms the "real" (base) nucleation rate J to
/// an "apparent" nucleation rate that accounts for the growth of nucleated
/// particles in critical clusters (CC) needed to place them into an
/// appropriate nucleation mode.
/// @param [in] c_so4 The number concentration of SO4 aerosol [#/cc]
/// @param [in] c_nh4 The number concentration of NH4 aerosol [#/cc]
/// @param [in] nh4_to_so4_molar_ratio The molar ratio of NH4 to SO4 [-]
/// @param [in] temp The atmospheric temperature [K]
/// @param [in] rel_hum The atmospheric relative humidity [-]
/// @param [in] d_dry_crit The dry diameter of particles in a CC [nm]
/// @param [in] d_wet_crit The wet diameter of particles in a CC [nm]
/// @param [in] d_dry_grown The dry diameter of grown particles [nm]
/// @param [in] rho_grown The mass density of grown particles [kg/m3]
/// @param [in] rho_air The mass density of dry air [kg/m3]
/// @param [in] mw_h2so4 The molecular weight of H2SO4 gas [kg/mol]
KOKKOS_INLINE_FUNCTION
PackType apparent_nucleation_factor(const PackType& c_so4,
                                    const PackType& c_nh4,
                                    const PackType& nh4_to_so4_molar_ratio,
                                    const PackType& temp,
                                    const PackType& rel_hum,
                                    const PackType& d_dry_crit,
                                    const PackType& d_wet_crit,
                                    const PackType& d_dry_grown,
                                    const PackType& rho_grown,
                                    const PackType& rho_air,
                                    Real mw_h2so4) {
  // Compute the wet/dry volume ratio using the simple Kohler approximation
  // for ammonium sulfate and bisulfate.
  const auto bounded_rel_hum = max(0.10, min(0.95, rel_hum));
  const auto wet_dry_vol_ratio = 1.0 - 0.56 / log(bounded_rel_hum);

    // Compute the fraction of the wet volume due to SO4 aerosol.
  PackType V_frac_wet_so4 = 1.0 /
      (wet_dry_vol_ratio * (1.0 + nh4_to_so4_molar_ratio * 17.0 / 98.0));

  // Compute the condensation growth rate gr [nm/h] of new particles from
  // KK2002 eq 21 for H2SO4 uptake and correct for NH3/H2O uptake.
  PackType speed = 14.7 * sqrt(temp);  // molecular speed [m/s]
  PackType gr = 3.0e-9 * speed * mw_h2so4 * c_so4 / (rho_grown * V_frac_wet_so4);

  //--------------------------------------------
  // Compute gamma from KK2002 eq 22 [nm2/m2/h]
  //--------------------------------------------

  // Wet diameter [nm] of grown particles with dry diameter d_dry_grown.
  PackType d_wet_grown = 1e9 * d_dry_grown * pow(wet_dry_vol_ratio, 1.0 / 3.0);

    // Compute gamma, neglecting the (d_mean/150)^0.048 factor.
  PackType gamma = 0.23 * pow(d_wet_crit, 0.2) *
                   pow(d_wet_grown / 3.0, 0.075) *
                   pow(1e-3 * rho_grown, -0.33) * pow(temp / 293.0, -0.75);

  // Compute the condensation sink CS' from KK2002 eqs 3-4. For the purposes
  // of this calculation, we use alpha == 1 and we use the mean free path of
  // air as computed from the air density in the calculation of the Knudsen
  // number for the nucleation mode.
  // NOTE: this differs from the MAM4 calculation, which uses an H2SO4
  // NOTE: uptake rate that assumes a process ordering, which we're no
  // NOTE: longer allowed to do.
  const Real alpha = 1;  // accommodation coefficient

  // The Knudsen number for the nucleated particles is Kn = 2 * lambda / d,
  // where lambda is the mean free path of air, and d is the grown particle
  // diameter. The mean free path is 1/(n * sigma), where n = rho_air/mw_air
  // is the number density of air, and sigma = pi*d^2 is the cross section
  // of a grown particle. Putting everything togther, we have
  //      2 * mw_air
  // Kn = --------------- 3
  //      pi * rho_air * d
  // TODO: should we attempt to estimate the wet number density?
  static const Real mw_air = Constants::molec_weight_dry_air;
  static const Real pi = Constants::pi;
  const PackType Kn = 2 * mw_air / (pi * rho_air * cube(d_wet_grown));

  // Compute the transitional correction for the condensational mass flux
  // (Fuchs and Sutugin, 1971, or KK2002 eq 4).
  const auto beta =
      (1.0 + Kn) / (1.0 + 0.377 * Kn + 1.33 * Kn * (1 + Kn) / alpha);

  // Compute the condensation sink for the nucleated particles from KK2002 eq 3.
  const auto cond_sink = 0.5 * d_wet_grown * beta * (c_so4 + c_nh4);

  // Compute eta [nm] using KK2002 eq 11.
  PackType eta = gamma * cond_sink / gr;

  return exp(eta / d_wet_grown - eta / d_wet_crit);
}

}  // namespace kerminen2002

}  // namespace haero
#endif
