#ifndef HAERO_DRIVER_ATMOSPHERE_HPP
#define HAERO_DRIVER_ATMOSPHERE_HPP

#include "haero/haero_config.hpp"
#include "haero/physical_constants.hpp"
#include "haero/floating_point.hpp"
#include "ekat/ekat_assert.hpp"
#include "kokkos/Kokkos_Core.hpp"
#include <cmath>
#include <string>

namespace haero {
namespace driver {

/** @defgroup ReferenceAtmosphere ReferenceAtmosphere

  Functions and classes related to generating model atmosphere vertical profiles.

  @{
*/

/** @brief This type and its associated functions
  define ambient atmospheric conditions for the driver's dynamics model.
*/
struct AtmosphericConditions {
  static constexpr Real pref = 100000; /// Reference pressure at z=0, Tv=Tv0 [Pa]
  /**  virtual temperature appx. factor [K]

    (equations (2.1) and (2.3) from Klemp & Wilhelmson, 1978, J. Atm. Sci. 35)
  */
  static constexpr Real alpha_v = 0.61;
  /// dry air kappa [1]
  static constexpr Real kappa = r_gas_dry_air_joule_per_k_per_kg/cp_dry_air_joule_per_k_per_kg;
  /// reference virtual potential temperature [K]
  Real Tv0; 
  /// initial virtual temperature lapse rate [K/m]
  Real Gammav; 
  /// maximum magnitude of vertical velocity [m/s]
  Real w0; 
  /// top of model [m]
  int ztop; 
  /// period of velocity oscillation [s]
  int tperiod; 
  /// initial water vapor mass mixing ratio at z = 0 [kg H<sub>2</sub>O / kg air]
  Real qv0; 
  /// initial decay rate of water vapor mass mixing ratio with height [per m]
  Real qv1; 

  /** Construct and return a hydrostatic instance of AtmosphericConditions

    This method checks that the input arguments are within reasonably expected
    bounds for the standard units listed below.

    @param [in] Tv0_ reference virtual temperature [K]
    @param [in] Gammav_ virtual temperature lapse rate [K/m]
    @param [in] w0_ maximum vertical velocity [m/s]
    @param [in] ztop_ model top [m]
    @param [in] tperiod_ period of velocity oscillation [s]
    @param [in] qv0_ water vapor mixing ratio at z = 0 [kg H<sub>2</sub>O / kg air]
    @param [in] qv1_ water vapor decay rate [1/m]
  */
  AtmosphericConditions(const Real Tv0_ = 300, const Real Gammav_ = 0.01, const Real w0_ = 1,
   const int ztop_ = 20E3,  const int tperiod_=900, const Real qv0_=1.5E-3, const Real qv1_ = 1E-3);
  
  KOKKOS_INLINE_FUNCTION
  AtmosphericConditions(const AtmosphericConditions& other) : Tv0(other.Tv0), Gammav(other.Gammav),
    w0(other.w0), ztop(other.ztop), tperiod(other.tperiod), qv0(other.qv0), qv1(other.qv1) {}
   
  /// Write instance info to string 
  std::string info_string(const int tab_level=0) const;
};


/** @brief Computes temperature from virtual temperature and water vapor mixing ratio,
@f$ T(T_v, q_v) @f$.

  @param [in] Tv virtual temperature [K]
  @param [in] qv water vapor mixing ratio [kg H<sub>2</sub>O / kg air]
  @return temperature [K]
*/
KOKKOS_INLINE_FUNCTION
Real temperature_from_virtual_temperature(const Real Tv, const Real qv) {
  return Tv / (1 + AtmosphericConditions::alpha_v * qv);
}

/** @brief Virtual temperature profile, @f$T_v(z)@f$.

  @param [in] z height [m]
  @param [in] T0 reference virtual temperature [K]
  @param [in] Gamma virtual temperature lapse rate [K/m]
  @return virtual temperature
*/
KOKKOS_INLINE_FUNCTION
Real virtual_temperature(const Real z, const Real Tv0, const Real Gammav) {
  return Tv0 - Gammav * z;
}


/** @brief Virtual temperature profile, @f$T_v(z)@f$.

  @param [in] z height [m]
  @param [in] conds AtmosphericConditions
  @return virtual temperature
*/
KOKKOS_INLINE_FUNCTION
Real virtual_temperature(const Real z, const AtmosphericConditions& conds) {
  return virtual_temperature(z, conds.Tv0, conds.Gammav);
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
  return water_vapor_mixing_ratio(z, conds.qv0, conds.qv1);
}

/** @brief Computes the hydrostatic pressure at a given height, based on a virtual temperature
  profile with constant lapse rate @f$ \Gamma = -\frac{\partial T_v}{\partial z} \ge 0 @f$.

  @param [in] z height [m]
  @param [in] p0 reference presssure [Pa]
  @param [in] T0 reference temperature [K]
  @param [in] Gamma temperature lapse rate [K/m]
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
  return hydrostatic_pressure_at_height(z, AtmosphericConditions::pref,
                                           conds.Tv0,
                                           conds.Gammav);
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
  return height_at_pressure(p, AtmosphericConditions::pref,
      conds.Tv0, conds.Gammav);
}

/** @brief Computes the potential temperature or virtual potential temperature

  @param [in] T temperature or virtual temperature [K]
  @param [in] p pressure [Pa]
  @param [in] p0 reference pressure [Pa]
  @return @f$ \theta @f$ or @f$ \theta_v @f$, depending on first argument @f$ T @f$ or @f$ T_v @f$
*/
KOKKOS_INLINE_FUNCTION
Real potential_temperature(const Real T, const Real p, const Real p0) {
  return T * std::pow(p0/p, AtmosphericConditions::kappa);
}


/** @brief Computes the potential temperature or virtual potential temperature

  @param [in] T temperature or virtual temperature [K]
  @param [in] p pressure [Pa]
  @param [in] conds AtmosphericConditions
  @return @f$ \theta @f$ or @f$ \theta_v @f$, depending on first argument @f$ T @f$ or @f$ T_v @f$
*/
KOKKOS_INLINE_FUNCTION
Real potential_temperature(const Real T, const Real p, const AtmosphericConditions& conds) {
  return potential_temperature(T, p, AtmosphericConditions::pref);
}

/** @brief computes the Exner pressure function @f$ \Pi = \left(\frac{p}{p_0}\right)^{\kappa}@f$,
  where @f$\kappa = R_d/c_p@f$ is a dry-air constant.

  @param [in] p [Pa]
*/
KOKKOS_INLINE_FUNCTION
Real exner_function(const Real p) {
  return std::pow(p/AtmosphericConditions::pref, AtmosphericConditions::kappa);
}

/// @} group AtmosphericConditions
} // namespace driver
} // namespace haero
#endif
