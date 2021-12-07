#ifndef HAERO_CONVERSIONS_HPP
#define HAERO_CONVERSIONS_HPP

#include "haero/constants.hpp"
#include "haero/haero.hpp"

/// This file contains functions for converting between various representations
/// of physical quantities in aerosol parameterizations.

namespace haero {

namespace conversions {

/// Given a number concentration for a species or mixture [m-3], computes and
/// returns a mass mixing ratio [kg species/kg dry air] based on its molecular
/// weight and on the density of dry air in the vicinity.
/// @param [in] number_conc The number concentration of the species/mixture
/// [m-3]
/// @param [in] molecular_wt The molecular weight of the species/mixture
/// [kg/kmol]
/// @param [in] dry_air_density The mass density of dry air [kg/m3]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar
mmr_from_number_conc(const Scalar& number_conc, Real molecular_wt,
                     const Scalar& dry_air_density) {
  const auto Na = Constants::avogadro;
  return number_conc * molecular_wt / (dry_air_density * Na);
}

/// Given a mass mixing ratio (mmr) for a species or mixture [kg species/kg
/// dry air], computes and returns a number density [m-3] based on its molecular
/// weight and on the density of dry air in the vicinity.
/// @param [in] mmr The number concentration of the species/mixture [m-3]
/// @param [in] molecular_wt The molecular weight of the species/mixture
/// [kg/kmol]
/// @param [in] dry_air_density The mass density of dry air [kg/m3]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar number_conc_from_mmr(
    const Scalar& mmr, Real molecular_wt, const Scalar& dry_air_density) {
  const auto Na = Constants::avogadro;
  return mmr * (dry_air_density * Na) / molecular_wt;
}

/// Given a molar mixing ratio (vmr) for a species or mixture
/// [kmol species/kmol dry air], computes and returns a mass mixing ratio
/// [kg species/kg dry air] based on its molecular weight.
/// @param [in] vmr The molar mixing ratio of the species/mixture [kmol/kmol
/// air]
/// @param [in] molecular_wt The molecular weight of the species/mixture
/// [kg/kmol]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar mmr_from_vmr(const Scalar& vmr,
                                           Real molecular_wt) {
  const auto mw_dry_air = Constants::molec_weight_dry_air;
  return vmr * molecular_wt / mw_dry_air;
}

/// Given a mass mixing ratio (mmr) for a species or mixture [kg species/kg
/// dry air], computes and returns a molar mixing ratio [kmol species/k dry air]
/// based on its molecular weight.
/// @param [in] mmr The molar mixing ratio of the species/mixture [kmol/kmol
/// dry air]
/// @param [in] molecular_wt The molecular weight of the species/mixture
/// [kg/kmol]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar vmr_from_mmr(const Scalar& mmr,
                                           Real molecular_wt) {
  const auto mw_dry_air = Constants::molec_weight_dry_air;
  return mmr * mw_dry_air / molecular_wt;
}

/// Computes the virtual temperature [K] from the temperature [K] and a water
/// vapor mass mixing ratio [kg vapor/kg dry air]. See Equation (3.1.15) from
/// Gill, 1982, Atmosphere-Ocean Dynamics, Academic Press, San Diego CA.
/// @param [in] Tv virtual temperature [K]
/// @param [in] q1 specific humidity [-]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar
virtual_temperature_from_temperature(const Scalar& T, const Scalar& q1) {
  return T * (1.0 + 0.6078 * q1);
}

/// Computes the temperature [K] from the virtual temperature [K] and a water
/// vapor mass mixing ratio [kg vapor/kg dry air].
/// @param [in] Tv virtual temperature [K]
/// @param [in] q1 specific humidity [-]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar
temperature_from_virtual_temperature(const Scalar& Tv, const Scalar& q1) {
  return Tv / (1.0 + 0.6078 * q1);
}

/// Computes the dry mass density from the total mass density using the specific
/// humidity.
/// @param [in] rho total mass density [kg/m3]
/// @param [in] q1 specific humidity [-]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar dry_from_total_mass_density(const Scalar& rho,
                                                          const Scalar& q1) {
  return rho * (1 - q1);
}

/// Computes the mass density of water vapor from the total mass density using
/// the specific humidity.
/// @param [in] rho total mass density [kg/m3]
/// @param [in] q1 specific humidity [-]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar vapor_from_total_mass_density(const Scalar& rho,
                                                            const Scalar& q1) {
  return rho * q1;
}

/// Computes the water vapor mixing ratio from the specific humidity. This
/// calculation diverges at q1 = 1.
/// @param [in] q1 specific humidity [-]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar
vapor_mixing_ratio_from_specific_humidity(const Scalar& q1) {
  return q1 / (1 - q1);
}

/// Computes the specific humidity from the water vapor mixing ratio. This
/// calculation diverges at q1 = 1.
/// @param [in] qv water vapor mixing ratio [kg vapor / kg dry air]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar
specific_humidity_from_vapor_mixing_ratio(const Scalar& qv) {
  return qv / (qv + 1);
}

/// Computes the saturation vapor pressure of water as a function
/// of temperature.
///
/// The formula is the improved Magnus formula from
///
///  O. A. Alduchov and R. E. Eskridge, 1996, Improved Magnus form approximation
///  of saturation vapor pressure, Journal of Applied Meteorology 35:601--609.
///
///  See eq. (21) from that paper, which improves the original formula's
///  accuracy over the temperature range [-40, 50] C expected of
///  atmospheric conditions throughout an entire vertical column.
///  The paper uses temperature in Celsius and returns
///  pressure in hPa. We have changed to SI units in this implementation.
///
///  See also Lamb & Verlinde section 3.3.
///
///  @param [in] T temperature [K]
///  @return es(T) saturation vapor pressure of water vapor [Pa]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar
vapor_saturation_pressure_magnus(const Scalar& T) {
  static constexpr Real e0 = 610.94;       // Pa
  static constexpr Real exp_num = 17.625;  // nondimensional
  static constexpr Real exp_den = 234.04;  // deg C
  const auto Tf = Constants::freezing_pt_h2o;
  const auto celsius_temp = T - Tf;
  return e0 * exp(exp_num * celsius_temp / (exp_den + celsius_temp));
}

/// Computes the relative humidity from the water vapor mixing ratio and the
/// pressure and temperature, given the relationship between temperature and
/// the water vapor saturation pressure.
/// @param [in] qv water vapor mixing ratio [kg vapor/kg dry air]
/// @param [in] p total pressure [Pa]
/// @param [in] T temperature [K]
/// @param [in] vsp A function that computes the vapor saturation pressure from
///                 the temperature. If not supplied,
///                 @ref vapor_saturation_pressure_magnus is used.
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar relative_humidity_from_vapor_mixing_ratio(
    const Scalar& qv, const Scalar& p, const Scalar& T,
    Scalar (*vsp)(const Scalar&) = vapor_saturation_pressure_magnus<Scalar>) {
  auto es = vsp(T);
  return qv / (es / p);
}

/// Computes the water vapor mixing ratio from the relative humidity and the
/// pressure and temperature, given the relationship between temperature and
/// the water vapor saturation pressure.
/// @param [in] rel_hum relative humidity [-]
/// @param [in] p total pressure [Pa]
/// @param [in] T temperature [K]
/// @param [in] vsp A function that computes the vapor saturation pressure from
///                 the temperature. If not supplied,
///                 @ref vapor_saturation_pressure_magnus is used.
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar vapor_mixing_ratio_from_relative_humidity(
    const Scalar& rel_hum, const Scalar& p, const Scalar& T,
    Scalar (*vsp)(const Scalar&) = vapor_saturation_pressure_magnus<Scalar>) {
  auto es = vsp(T);
  return rel_hum * es / p;
}

}  // namespace conversions

}  // namespace haero

#endif
