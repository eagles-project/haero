#ifndef HAERO_CONVERSIONS_HPP
#define HAERO_CONVERSIONS_HPP

#include "haero/haero.hpp"
#include "haero/physical_constants.hpp"

/// This file contains functions for converting between various representations
/// of physical quantities in aerosol parameterizations.

namespace haero {

namespace conversions {

/// Given a number concentration for a species or mixture [m-3], computes and
/// returns a mass mixing ratio [kg species/kg air] based on its molecular
/// weight and on the density of dry air in the vicinity.
/// @param [in] number_conc The number concentration of the species/mixture
/// [m-3]
/// @param [in] molecular_wt The molecular weight of the species/mixture
/// [kg/kmol]
/// @param [in] air_density The mass density of dry air [kg/m3]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar mmr_from_number_conc(const Scalar& number_conc, Real molecular_wt,
                            const Scalar& air_density) {
  return number_conc * molecular_wt / (air_density * constants::avogadro);
}

/// Given a mass mixing ratio (mmr) for a species or mixture [kg species/kg
/// air], computes and returns a number density [m-3] based on its molecular
/// weight and on the density of dry air in the vicinity.
/// @param [in] mmr The number concentration of the species/mixture [m-3]
/// @param [in] molecular_wt The molecular weight of the species/mixture
/// [kg/kmol]
/// @param [in] air_density The mass density of dry air [kg/m3]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar number_conc_from_mmr(const Scalar& mmr, Real molecular_wt,
                            const Scalar& air_density) {
  return mmr * (air_density * constants::avogadro) / molecular_wt;
}

/// Given a molar mixing ratio (vmr) for a species or mixture
/// [kmol species/kmol air], computes and returns a mass mixing ratio
/// [kg species/kg air] based on its molecular weight.
/// @param [in] vmr The molar mixing ratio of the species/mixture [kmol/kmol
/// air]
/// @param [in] molecular_wt The molecular weight of the species/mixture
template <typename Scalar>
/// [kg/kmol]
KOKKOS_INLINE_FUNCTION
Scalar mmr_from_vmr(const Scalar& vmr, Real molecular_wt) {
  return vmr * molecular_wt / constants::molec_weight_dry_air;
}

/// Given a mass mixing ratio (mmr) for a species or mixture [kg species/kg
/// air], computes and returns a molar mixing ratio [kmol species/k air] based
/// on its molecular weight.
/// @param [in] mmr The molar mixing ratio of the species/mixture [kmol/kmol
/// air]
/// @param [in] molecular_wt The molecular weight of the species/mixture
/// [kg/kmol]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar vmr_from_mmr(const Scalar& mmr, Real molecular_wt) {
  return mmr * constants::molec_weight_dry_air / molecular_wt;
}

/// Computes the virtual temperature [K] from the temperature [K] and a water
/// vapor mass mixing ratio [kg vapor/kg dry air].
/// @param [in] Tv virtual temperature [K]
/// @param [in] q1 specific humidity [-]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar virtual_temperature_from_temperature(const Scalar& T, const Scalar& q1) {
  return T * (1.0 + 0.6078 * q1);
}

/// Computes the temperature [K] from the virtual temperature [K] and a water
/// vapor mass mixing ratio [kg vapor/kg dry air].
/// @param [in] Tv virtual temperature [K]
/// @param [in] q1 specific humidity [-]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar temperature_from_virtual_temperature(const Scalar& Tv,
                                            const Scalar& q1) {
  return Tv / (1.0 + 0.6078 * q1);
}

/// Computes the dry mass density from the total mass density using the specific
/// humidity.
/// @param [in] rho total mass density [kg/m3]
/// @param [in] q1 specific humidity [-]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar dry_from_total_mass_density(const Scalar& rho, const Scalar& q1) {
  return rho * (1 - q1);
}

/// Computes the mass density of water vapor from the total mass density using
/// the specific humidity.
/// @param [in] rho total mass density [kg/m3]
/// @param [in] q1 specific humidity [-]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar vapor_from_total_mass_density(const Scalar& rho, const Scalar& q1) {
  return rho * q1;
}

/// Computes the water vapor mixing ratio from the specific humidity. This
/// calculation diverges at q1 = 1.
/// @param [in] q1 specific humidity [-]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar vapor_mixing_ratio_from_specific_humidity(const Scalar& q1) {
  return q1 / (1 - q1);
}

/// Computes the specific humidity from the water vapor mixing ratio. This
/// calculation diverges at q1 = 1.
/// @param [in] qv water vapor mixing ratio [kg vapor / kg dry air]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar specific_humidity_from_vapor_mixing_ratio(const Scalar& qv) {
  return qv / (qv + 1);
}

/// Computes the water vapor saturation pressure from a temperature using the
/// Tetens formula from Soong-Ogura 1973 equation (A1) or Klemp-Wilhelmson 1978
/// eqn. (2.11).
/// @param [in] T temperature [K]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar vapor_saturation_pressure_tetens(const Scalar& T) {
  static constexpr Real half15ln10 = 17.269388197455342630;
  static constexpr Real tetens_coeff = 380.042;
  return tetens_coeff * exp(Scalar(half15ln10) * (T - 273) / (T - 36));
}

/// Computes the relative humidity from the water vapor mixing ratio and the
/// pressure and temperature, given the relationship between temperature and
/// the water vapor saturation pressure.
/// @param [in] qv water vapor mixing ratio [kg vapor/kg dry air]
/// @param [in] p total pressure [Pa]
/// @param [in] T temperature [K]
/// @param [in] vsp A function that computes the vapor saturation pressure from
///                 the temperature. If not supplied,
///                 @ref vapor_saturation_pressure_tetens is used.
template <typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar relative_humidity_from_vapor_mixing_ratio(
    const Scalar& qv, const Scalar& p, const Scalar& T,
    const std::function<Scalar(const Scalar&)>& vsp =
        vapor_saturation_pressure_tetens<Scalar>) {
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
///                 @ref vapor_saturation_pressure_tetens is used.
template <typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar vapor_mixing_ratio_from_relative_humidity(
    const Scalar& rel_hum, const Scalar& p, const Scalar& T,
    const std::function<Scalar(const Scalar&)>& vsp =
        vapor_saturation_pressure_tetens<Scalar>) {
  auto es = vsp(T);
  return rel_hum * es / p;
}

}  // namespace conversions

}  // namespace haero

#endif
