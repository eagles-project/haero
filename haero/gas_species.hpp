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

/// This factory function constructs a set of gas species corresponding to
/// the legacy MAM4 model.
///
/// These specific species names can be found in an e3sm v1 atm.log file.
///
/// for info on data sources, see
/// https://eagles-project.atlassian.net/wiki/spaces/Computation/pages/1125515265/Aerosol+species+and+mode+data
inline std::vector<GasSpecies> create_mam4_gas_species() {
  const std::vector<std::string> gnames = {
      "ozone",
      "hydrogen peroxide",
      "sulfuric acid",
      "sulfur dioxide",
      "dimethylsulfide",
      "secondary organic aerosol precursor",
      "oxygen",
      "carbon dioxide",
      "nitrous oxide",
      "methane",
      "trichlorofluoromethane",
      "dichlorodifluoromethane",
      "NH3"};
  const std::vector<std::string> gsymbs = {
      "O3",  "H2O2", "H2SO4", "SO2",   "DMS",   "SOAG", "O2",
      "CO2", "N2O",  "CH4",   "CFC11", "CFC12", "NH3"};
  const std::vector<Real> gmws = {47.9982, 34.0136, 98.0784, 64.0648, 62.1324,
                                  12.011,  31.988,  44.009,  44.013,  16.04,
                                  137.37,  120.91,  50.0};
  static constexpr Real g_to_kg = 0.001;
  std::vector<GasSpecies> result;
  for (int i = 0; i < gnames.size(); ++i) {
    result.push_back(GasSpecies(gnames[i], gsymbs[i], "(No description)",
                                g_to_kg * gmws[i]));
  }
  return result;
}

}  // namespace haero
#endif
