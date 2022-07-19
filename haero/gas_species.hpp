#ifndef HAERO_GAS_SPECIES_HPP
#define HAERO_GAS_SPECIES_HPP

#include <limits>
#include <map>
#include <string>
#include <vector>

#include "haero/haero.hpp"

namespace haero {

/// @struct GasSpecies
/// This type represents a gas that participates in one or more aerosol
/// microphysics parameterizations.
struct GasSpecies final {
  // Default constructor needed for device
  KOKKOS_INLINE_FUNCTION
  GasSpecies() {
  }

  /// Creates a new gas species
  /// @param [in] molecular_wt the molecular weight of this species [kg/mol]
  KOKKOS_INLINE_FUNCTION
  explicit GasSpecies(Real molecular_wt)
      : molecular_weight(molecular_wt) {
  }

  KOKKOS_INLINE_FUNCTION
  GasSpecies(const GasSpecies& g) : molecular_weight(g.molecular_weight) {
  }

  KOKKOS_INLINE_FUNCTION
  GasSpecies& operator=(const GasSpecies& g) {
    if (&g != this) {
      molecular_weight = g.molecular_weight;
    }
    return *this;
  }

  /// Molecular weight [kg/mol]
  Real molecular_weight;

  // Comparison operators.
  bool operator==(const GasSpecies& other) const {
    return (molecular_weight == other.molecular_weight);
  }
  bool operator!=(const GasSpecies& other) const { return !(*this == other); }
};

}  // namespace haero
#endif
