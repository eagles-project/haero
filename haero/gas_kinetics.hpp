#ifndef HAERO_GAS_KINETICS_HPP
#define HAERO_GAS_KINETICS_HPP

#include "haero/conversions.hpp"
#include "haero/haero.hpp"
#include "haero/physical_constants.hpp"

/// This file contains functions that relate various quantities related to the
/// kinetic theory of gases.

namespace haero {

namespace gas_kinetics {

/// Computes the total mass density of air using the ideal gas law.
/// @param [in] p The total pressure of a parcel of air [Pa]
/// @param [in] T The temperature of the air parcel [K]
/// @param [in] qv The water vapor mixing ratio of the air parcel [-]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar air_mass_density(const Scalar& p, const Scalar& T,
                                               const Scalar& qv) {
  auto q1 = conversions::specific_humidity_from_vapor_mixing_ratio(qv);
  auto Tv = conversions::virtual_temperature_from_temperature(T, q1);
  return p / (constants::r_gas_dry_air * Tv);
}

/// Computes the molecular speed of an ideal gas at the given temperature with
/// the given molecular weight and the given adiabatic constant.
/// @param [in] T The temperature of the gas [K]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar molecular_speed(const Scalar& T,
                                              Real molecular_wt, Real gamma) {
  return sqrt(gamma * constants::r_gas * T / molecular_wt);
}

/// Computes the mean free path [m] of a gas given the diameter of its molecules
/// and its number concentration.
/// @param [in] diameter The diameter of a molecule in the gas [m]
/// @param [in] num_conc The number concentration of the gas [#/m3]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar mean_free_path(const Scalar& diameter,
                                             const Scalar& num_conc) {
  static constexpr double sqrt2 = std::sqrt(2);
  return 1.0 / (sqrt2 * constants::pi * square(diameter) * num_conc);
}

/// Computes the Knudsen number associated with a gas flow with the given mean
/// free path over the given characteristic scale length.
/// @param [in] mean_free_path The mean free path of molecules in the gas [m]
/// @param [in] scale_length The characteristic scale length of the flow [m]
template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar knudsen_number(const Scalar& mean_free_path,
                                             const Scalar& scale_length) {
  return mean_free_path / scale_length;
}

}  // namespace gas_kinetics

}  // namespace haero

#endif
