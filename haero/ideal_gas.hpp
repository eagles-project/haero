#ifndef HAERO_IDEAL_GAS_HPP
#define HAERO_IDEAL_GAS_HPP

#include "haero/conversions.hpp"
#include "haero/haero.hpp"
#include "haero/physical_constants.hpp"

/// This file contains functions that relate various quantities using the ideal
/// gas law.

namespace haero {

namespace ideal_gas {

/// Computes the total mass density of air given the total pressure, the
/// temperature, and the specific humidity.
/// @param [in] p The total pressure of a parcel of air [Pa]
/// @param [in] T The temperature of the air parcel [K]
/// @param [in] qv The water vapor mixing ratio of the air parcel [-]
KOKKOS_INLINE_FUNCTION
PackType mass_density(const PackType& p, const PackType& T,
                      const PackType& qv) {
  auto q1 = conversions::specific_humidity_from_vapor_mixing_ratio(qv);
  auto Tv = conversions::virtual_temperature_from_temperature(T, q1);
  return p / (constants::r_gas_dry_air * Tv);
}

}  // namespace ideal_gas

}  // namespace haero

#endif
