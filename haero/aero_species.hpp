#ifndef HAERO_AERO_SPECIES_HPP
#define HAERO_AERO_SPECIES_HPP

#include <limits>
#include <map>
#include <string>
#include <vector>

#include "haero/haero.hpp"

namespace haero {

/// @struct AeroSpecies
/// This type represents an aerosol species.
struct AeroSpecies final {
  // Default constructor needed to resize Kokkos Views on device before deep
  // copy.
  KOKKOS_INLINE_FUNCTION
  AeroSpecies() {
  }

  /// Creates a new aerosol species.
  /// @param [in] molecular_wt The molecular weight [kg/mol]of the species
  /// @param [in] dry_rad The dry radius [m] of the species' particle size
  /// @param [in] hygro Base hygroscopicity of the species
  KOKKOS_INLINE_FUNCTION
  AeroSpecies(Real molecular_wt, Real dens, Real hygro)
      : molecular_weight(molecular_wt), density(dens), hygroscopicity(hygro) {
  }

  KOKKOS_INLINE_FUNCTION
  AeroSpecies(const AeroSpecies& a)
      : molecular_weight(a.molecular_weight),
        density(a.density),
        hygroscopicity(a.hygroscopicity) {
  }

  KOKKOS_INLINE_FUNCTION
  AeroSpecies& operator=(const AeroSpecies& a) {
    if (&a != this) {
      molecular_weight = a.molecular_weight;
      density = a.density;
      hygroscopicity = a.hygroscopicity;
    }
    return *this;
  }

  // Molecular weight [kg/mol]
  Real molecular_weight;

  /// Material density [kg/m^3]
  Real density;

  /// Hygroscopicity
  Real hygroscopicity;

  // Comparison operators.
  bool operator==(const AeroSpecies& other) const {
    return ((molecular_weight == other.molecular_weight) and
            (density == other.density) and
            (hygroscopicity == other.hygroscopicity));
  }
  bool operator!=(const AeroSpecies& other) const {
    return !(*this == other);
  }
};

}  // namespace haero

#endif
