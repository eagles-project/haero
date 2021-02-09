#ifndef HAERO_DRIVER_ATMOSPHERE_HPP
#define HAERO_DRIVER_ATMOSPHERE_HPP

#include "haero/haero_config.hpp"
#include "haero/physical_constants.hpp"
#include "haero/floating_point.hpp"
#include "ekat/ekat_assert.hpp"
#include "kokkos/Kokkos_Core.hpp"
#include <cmath>

namespace haero {
namespace driver {



/** @defgroup ReferenceAtmosphere ReferenceAtmosphere

  Functions and classes related to generating model atmosphere vertical profiles.

  @{
*/

/** @brief This type and its associated functions
  define ambient atmospheric conditions for supported models.

   The model used to generate the atmospheric conditions.
   * uniform: all vertical levels have identical atmospheric conditions
   * hydrostatic: the vertical profile of the atmospheric conditions is
                  determined using the hydrostatic equilibrium approximation.
*/
struct AtmosphericConditions {

  enum { uniform, hydrostatic } model;
  /// Parameters specific to each model.
  union {
    /// Uniform model.
    struct {
      /// Mean molecular weight of air [kg/mol].
      Real mu;
      /// Scaled atmospheric height [m].
      Real H;
      /// Pressure [Pa].
      Real p0;
      /// Temperature [K].
      Real T0;
      /// Clear-sky relative humidity [-]
      Real phi0;
      /// Cloud fraction [-]
      Real N0;
    } uniform;
    /// Hydrostatic model (we assume a constant virtual temperature lapse rate).
    struct {
      /// Reference pressure at z = 0 [Pa].
      Real p0;
      /// Reference virtual temperature [K] at z = 0, p = p_0.
      Real T0;
      /// Virtual temperature lapse rate @f$ \Gamma = -\partial T_v / \partial z @f$ [K/m]
      Real lapse_rate;
      /// Water vapor mixing ratio [kg H<sub>2</sub>O / kg air] at surface
      Real qv0;
      /// Water vapor decay rate with height [1/m]
      Real qv1;
    } hydrostatic;
  } params;
};

/** Construct and return a hydrostatic instance of AtmosphericConditions

  This method checks that the input arguments are within reasonably expected
  bounds for the standard units listed below.

  @param [in] p0 reference pressure [Pa]
  @param [in] T0 reference virtual temperature [K]
  @param [in] Gamma virtual temperature lapse rate [K/m]
  @param [in] qv0 water vapor mixing ratio at z = 0 [kg H<sub>2</sub>O / kg air]
  @param [in] qv1 water vapor decay rate [1/m]
*/
AtmosphericConditions hydrostatic_conditions(const Real p0, const Real T0, const Real Gamma,
                                             const Real qv0, const Real qv1);

/**  virtual temperature appx. factor [K]

  (equations (2.1) and (2.3) from Klemp & Wilhelmson, 1978, J. Atm. Sci. 35)
*/
static constexpr Real alpha_v = 0.61;

/// dry air kappa [1]
static constexpr Real kappa = r_gas_dry_air_joule_per_k_per_kg/cp_dry_air_joule_per_k_per_kg;

/** @brief Computes temperature from virtual temperature and water vapor mixing ratio,
@f$ T(T_v, q_v) @f$.

  @param [in] Tv virtual temperature [K]
  @param [in] qv water vapor mixing ratio [kg H<sub>2</sub>O / kg air]
  @return temperature [K]
*/
KOKKOS_INLINE_FUNCTION
Real temperature_from_virtual_temperature(const Real Tv, const Real qv) {
  return Tv / (1 + alpha_v * qv);
}

/** @brief Virtual temperature profile, @f$T_v(z)@f$.

  @param [in] z height [m]
  @param [in] T0 reference virtual temperature [K]
  @param [in] Gamma virtual temperature lapse rate [K/m]
  @return virtual temperature
*/
KOKKOS_INLINE_FUNCTION
Real virtual_temperature(const Real z, const Real T0, const Real Gamma) {
  return T0 - Gamma * z;
}


/** @brief Virtual temperature profile, @f$T_v(z)@f$.

  @param [in] z height [m]
  @param [in] conds AtmosphericConditions
  @return virtual temperature
*/
KOKKOS_INLINE_FUNCTION
Real virtual_temperature(const Real z, const AtmosphericConditions& conds) {
  EKAT_KERNEL_ASSERT(conds.model == AtmosphericConditions::hydrostatic);
  return virtual_temperature(z, conds.params.hydrostatic.T0, conds.params.hydrostatic.lapse_rate);
}

/** @brief Water vapor mixing ratio profile, @f$ q_v(z) @f$

  @param [in] z height [m]
  @param [in] qv0 mixing ratio at surface [kg H<sub>2</sub>O / kg air]
  @param [in] qv1 decay rate of water vapor [1/m]
  @return @f$ q_v @f$
*/
KOKKOS_INLINE_FUNCTION
Real water_vapor_mixing_ratio(const Real z, const Real qv0, const Real qv1) {
  return qv0 * std::exp(-qv1 * z);
}

/** @brief Water vapor mixing ratio profile, @f$ q_v(z) @f$

  @f$ q_v(z) = q_v^{(0)} e^{-q_v^{(1)}z} @f$

  @param [in] z height [m]
  @param [in] conds AtmosphericConditions
  @return @f$q_v(z)@f$ [kg H<sub>2</sub>O / kg air]
*/
KOKKOS_INLINE_FUNCTION
Real water_vapor_mixing_ratio(const Real z, const AtmosphericConditions& conds) {
  EKAT_KERNEL_ASSERT(conds.model == AtmosphericConditions::hydrostatic);
  return water_vapor_mixing_ratio(z, conds.params.hydrostatic.qv0,
                                     conds.params.hydrostatic.qv1);
}

/** @brief Computes the hydrostatic pressure at a given height, based on a virtual temperature
  profile with constant lapse rate @f$ \Gamma = -\frac{\partial T_v}{\partial z} \ge 0 @f$.

  @param [in] z height [m]
  @param [in] p0 reference presssure [Pa]
  @param [in] T0 reference virtual temperature [K]
  @param [in] Gamma virtual temperature lapse rate [K/m]
  @return p [Pa]
*/
KOKKOS_INLINE_FUNCTION
Real hydrostatic_pressure_at_height(const Real z, const Real p0, const Real T0, const Real Gamma) {
  Real result = 0;
  if (FloatingPoint<Real>::zero(Gamma)) {
    result = p0 * std::exp(-gravity_m_per_s2 * z / (r_gas_dry_air_joule_per_k_per_kg * T0));
  }
  else {
    const Real pwr = gravity_m_per_s2 / (r_gas_dry_air_joule_per_k_per_kg * Gamma);
    result = p0 * std::pow(T0, -pwr)*std::pow(T0 - Gamma*z, pwr);
  }
  return result;
}

/** @brief Computes the hydrostatic pressure at a given height, based on a virtual temperature
  profile with constant lapse rate @f$ \Gamma = -\frac{\partial T_v}{\partial z} \ge 0 @f$.

  @param [in] z height [m]
  @param [in] conds AtmosphericConditions
  @return pressure [Pa]
*/
KOKKOS_INLINE_FUNCTION
Real hydrostatic_pressure_at_height(const Real z, const AtmosphericConditions& conds) {
  EKAT_KERNEL_ASSERT(conds.model == AtmosphericConditions::hydrostatic);
  return hydrostatic_pressure_at_height(z, conds.params.hydrostatic.p0,
                                           conds.params.hydrostatic.T0,
                                           conds.params.hydrostatic.lapse_rate);
}

/** @brief Computes the height based on hydrostatic pressure.

    @param [in] p pressure [Pa]
    @param [in] p0 reference pressure [Pa]
    @param [in] T0 reference virtual temperature
    @param [in] Gamma virtual temperature lapse rate, @f$\Gamma = -\frac{\partial T_v}{\partial z} \ge 0 @f$
    @return height [m]
*/
KOKKOS_INLINE_FUNCTION
Real height_at_pressure(const Real p, const Real p0, const Real T0, const Real Gamma) {
  Real result = 0;
  if (FloatingPoint<Real>::zero(Gamma)) {
    result = -r_gas_dry_air_joule_per_k_per_kg * T0 * std::log(p/p0)/gravity_m_per_s2;
  }
  else {
    const Real pwr = r_gas_dry_air_joule_per_k_per_kg*Gamma/gravity_m_per_s2;
    result = (T0/Gamma)*(1 - std::pow(p/p0, pwr));
  }
  return result;
}

/** @brief Computes the height based on hydrostatic pressure.

    @param [in] p pressure [Pa]
    @param [in] conds AtmosphericConditions
    @return height [m]
*/
KOKKOS_INLINE_FUNCTION
Real height_at_pressure(const Real p, const AtmosphericConditions& conds) {
  EKAT_KERNEL_ASSERT(conds.model == AtmosphericConditions::hydrostatic);
  return height_at_pressure(p, conds.params.hydrostatic.p0,
      conds.params.hydrostatic.T0, conds.params.hydrostatic.lapse_rate);
}

/** @brief Computes the potential temperature or virtual potential temperature

  @param [in] T temperature or virtual temperature [K]
  @param [in] p pressure [Pa]
  @param [in] p0 reference pressure [Pa]
  @return @f$ \theta @f$ or @f$ \theta_v @f$, depending on first argument @f$ T @f$ or @f$ T_v @f$
*/
KOKKOS_INLINE_FUNCTION
Real potential_temperature(const Real T, const Real p, const Real p0) {
  return T * std::pow(p0/p, kappa);
}


/** @brief Computes the potential temperature or virtual potential temperature

  @param [in] T temperature or virtual temperature [K]
  @param [in] p pressure [Pa]
  @param [in] conds AtmosphericConditions
  @return @f$ \theta @f$ or @f$ \theta_v @f$, depending on first argument @f$ T @f$ or @f$ T_v @f$
*/
KOKKOS_INLINE_FUNCTION
Real potential_temperature(const Real T, const Real p, const AtmosphericConditions& conds) {
  EKAT_KERNEL_ASSERT(conds.model == AtmosphericConditions::hydrostatic);
  return potential_temperature(T, p, conds.params.hydrostatic.p0);
}

/** @brief computes the Exner pressure function @f$ \Pi = \left(\frac{p}{p_0}\right)^{\kappa}@f$,
  where @f$\kappa = R_d/c_p@f$ is a dry-air constant.

  @param [in] p [Pa]
  @param [in] conds
*/
KOKKOS_INLINE_FUNCTION
Real exner_function(const Real p, const AtmosphericConditions& conds) {
  EKAT_KERNEL_ASSERT(conds.model == AtmosphericConditions::hydrostatic);
  return std::pow(p/conds.params.hydrostatic.p0, kappa);
}

/// @} group AtmosphericConditions
} // namespace driver
} // namespace haero
#endif