#ifndef HAERO_PHYSICAL_CONSTANTS_HPP
#define HAERO_PHYSICAL_CONSTANTS_HPP

#include "haero/haero_config.hpp"

namespace haero {

/** @defgroup double_precision_required double_precision_required
    @{
*/
/// Avogadro's constant [# / kmol]
static constexpr double avogadro_molec_per_kmole = 6.022214E26;
/// Avogadro's constant [$ / mol]
static constexpr double avogadro_molec_per_mole = 1.0E-3 * avogadro_molec_per_kmole;
/// Boltzmann's constant [J/(K #)]
static constexpr double boltzmann_joule_per_k_per_molec = 1.38065e-23;
/// pi
static constexpr double pi = 3.14159265358979323846264;
/// @}

/** @defgroup single_precision_allowed working_precision
    @{
*/
/// Universal gas constant [J / (K kmol)]
static constexpr Real r_gas_joule_per_k_per_kmole = avogadro_molec_per_kmole*boltzmann_joule_per_k_per_molec;
/// Universal gas constant [J / (K mol)]
static constexpr Real r_gas_joule_per_k_per_mole = 1.0E-3 * r_gas_joule_per_k_per_kmole;
/// acceleration of gravity [m/s^2]

/// acceleration of gravity [m/s^2]
static constexpr Real gravity_m_per_s2 = 9.80616;

/// Molecular weight of water [kg/mole *or* g/mole]
static constexpr Real molec_weight_h2o_g_per_mole = 18.016;

/// Molecular weight of dry air [kg/mole *or* g/mole]
static constexpr Real molec_weight_dry_air_g_per_mole = 28.966;

/// Molecular weight of sulfate ion @f$\text{SO}_4^{2-}@f$ [kg/kmole *or* g/mole]
/// @todo Wikipedia defines this with 4 significant figures, not 2
static constexpr Real molec_weight_so4_g_per_mole = 96.0;

/// Molecular weight of ammonium ion @f$\text{NH}_4^+@f$ [kg/mole *or* g/mole]
static constexpr Real molec_weight_nh4_g_per_mole = 18.0;

/// Density of fresh water [kg/m^3]
static constexpr Real rho_h20_fresh_kg_per_m3 = 1.0e3;

/// Standard pressure [Pa]
static constexpr Real p_std_pa = 101325.0;

/// Freezing point of freshwater [K]
static constexpr Real t_freeze_h2o_k = 273.15;

/// Melting point of freshwater [K]
static constexpr Real t_melt_h2o_k = t_freeze_h2o_k;

/// Molecular diffusion volume of dry air [dimensionless]
static constexpr Real molec_diffuse_dry_air = 20.1;

/// Water to air weight ratio [dimensionless]
static constexpr Real epsilon_h2o_air = molec_weight_h2o_g_per_mole / molec_weight_dry_air_g_per_mole;

/// Latent heat of evaporation [J/kg]
static constexpr Real latent_heat_evap_joule_per_kg = 2.501e6;

/// Latent heat of fustion [J/kg]
static constexpr Real latent_heat_fusion_joule_per_kg = 3.337e5;

/// Water vapor gas constant [J/K/kg]
static constexpr Real r_gas_h2o_vapor_joule_per_k_per_kg = r_gas_joule_per_k_per_kmole / molec_weight_h2o_g_per_mole;

/// Dry air gas constant [J/K/kg]
static constexpr Real r_gas_dry_air_joule_per_k_per_kg = r_gas_joule_per_k_per_kmole / molec_weight_dry_air_g_per_mole;


/// Specific heat (constant pressure) of dry air [J/kg/K]
static constexpr Real cp_dry_air_joule_per_k_per_kg = 1.00464e3;

/** Surface tension at water-air interface at 273K [mN/m]

   @todo
   (Wikipedia)[https://en.wikipedia.org/wiki/Surface_tension]
    gives additional significant figures: 75.64
*/
static constexpr Real surften_sigma_h2o_air_273k_mnewton_per_m = 76.0;

/// Dry adiabatic lapse rate [K/m]
static constexpr Real dry_adiabatic_lapse_rate_K_per_m = 0.0098;

/// @}

} // namespace haero
#endif
