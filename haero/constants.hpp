#ifndef HAERO_PHYSICAL_CONSTANTS_HPP
#define HAERO_PHYSICAL_CONSTANTS_HPP

// clang-format off
// All physical constants in Haero are in SI units unless it's explicitly
// stated otherwise. Units given in brackets are SI units and are included
// here only for clarity.
//
// Implementation note: we store constants in a struct so their static
// declarations work in a CUDA environment.
#include "haero/haero_config.hpp"

namespace haero {

struct Constants {

/** @defgroup Physical_constants_double_precision_required Physical constants (double_precision)

    Force-double precision for constants in this group.

    @{
*/
/// pi
static constexpr double pi = 3.14159265358979323846264;
/// 1/6th of the pi (required for lognormally distributed particle size/volume calculations)
static constexpr double pi_sixth = pi/6;
/// Avogadro's constant [#/mol]
static constexpr double avogadro = 6.022214076e23;
/// Boltzmann's constant [J/K]
static constexpr double boltzmann = 1.380649e-23;
/// @}

/** @defgroup Physical_constants_single_precision_allowed Physical constants (working_precision)

    Configurable precision (single or double) constants.

    @{
*/
/// Universal gas constant [J/(K mol)]
static constexpr Real r_gas = avogadro*boltzmann;

/// acceleration of gravity [m/s^2]
static constexpr Real gravity = 9.80616;

/// Molecular weight of water [kg/mol]
static constexpr Real molec_weight_h2o = 0.18016;

/// Molecular weight of dry air [kg/mol]
static constexpr Real molec_weight_dry_air = 0.028966;

/// Molecular weight of sulfate ion @f$\text{SO}_4^{2-}@f$ [kg/mol]
static constexpr Real molec_weight_so4 = 0.09606;

/// Molecular weight of ammonium ion @f$\text{NH}_4^+@f$ [kg/mol]
static constexpr Real molec_weight_nh4 = 0.018039;

/// Mass density of water [kg/m^3]
static constexpr Real density_h2o = 1.0e3;

/// Pressure at standard conditions (STP) [Pa]
static constexpr Real pressure_stp = 101325.0;

/// Freezing point of water [K]
static constexpr Real freezing_pt_h2o = 273.15;

/// Melting point of water [K]
static constexpr Real melting_pt_h2o = freezing_pt_h2o;

/// Molecular diffusion volume of dry air [-]
static constexpr Real molec_diffusion_dry_air = 20.1;

/// Water-to-dry-air weight ratio [-]
static constexpr Real weight_ratio_h2o_air = molec_weight_h2o / molec_weight_dry_air;

/// Latent heat of evaporation [J/kg]
static constexpr Real latent_heat_evap = 2.501e6;

/// Latent heat of fustion [J/kg]
static constexpr Real latent_heat_fusion = 3.337e5;

/// Water vapor gas constant [J/K/kg]
static constexpr Real r_gas_h2o_vapor = r_gas / molec_weight_h2o;

/// Dry air gas constant [J/K/kg]
static constexpr Real r_gas_dry_air = r_gas / molec_weight_dry_air;

/// Specific heat (at constant pressure) of dry air [J/kg/K]
static constexpr Real cp_dry_air = 1.00464e3;

/// Surface tension at water-air interface at 273K [N/m]
static constexpr Real surface_tension_h2o_air_273k = 0.07564;

/// Dry adiabatic lapse rate [K/m]
static constexpr Real dry_adiabatic_lapse_rate = 0.0098;

/// @}

}; // struct Constants

constexpr double Constants::avogadro;
constexpr Real Constants::r_gas_dry_air;

} // namespace haero
// clang-format on
#endif
