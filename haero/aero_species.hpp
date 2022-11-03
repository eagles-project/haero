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
  // Molecular weight [kg/mol]
  const Real molecular_weight;

  /// Material density [kg/m^3]
  const Real density;

  /// Hygroscopicity
  const Real hygroscopicity;
};

}  // namespace haero

#endif
