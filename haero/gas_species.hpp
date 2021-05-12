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
  static const int NAME_LEN = 128;
  static const int DESC_LEN = 512;
  // Default constructor needed for device
  KOKKOS_INLINE_FUNCTION
  GasSpecies() {
    name_view[0] = '\0';
    symbol_view[0] = '\0';
    desc_view[0] = '\0';
  }

  /// Creates a new gas species
  /// @param [in] name A unique descriptive name for this species
  /// @param [in] symbol a unique short name or symbol for this species
  /// @param [in] description A short text description of this species and its
  /// composition
  /// @param [in] molecular_wt the molecular weight of this species [kg/mol]
  GasSpecies(const std::string& name, const std::string& symbol,
             const std::string& description, Real molecular_wt)
      : molecular_weight(molecular_wt) {
    EKAT_ASSERT(name.size() < NAME_LEN);
    EKAT_ASSERT(symbol.size() < NAME_LEN);
    strncpy(name_view, name.c_str(), NAME_LEN);
    strncpy(symbol_view, symbol.c_str(), NAME_LEN);
    strncpy(desc_view, description.c_str(), DESC_LEN);
  }

  KOKKOS_INLINE_FUNCTION
  GasSpecies(const GasSpecies& g) : molecular_weight(g.molecular_weight) {
    for (int i = 0; i < NAME_LEN; ++i) name_view[i] = g.name_view[i];
    for (int i = 0; i < NAME_LEN; ++i) symbol_view[i] = g.symbol_view[i];
    for (int i = 0; i < DESC_LEN; ++i) desc_view[i] = g.desc_view[i];
  }

  KOKKOS_INLINE_FUNCTION
  GasSpecies& operator=(const GasSpecies& g) {
    if (&g != this) {
      molecular_weight = g.molecular_weight;
      for (int i = 0; i < NAME_LEN; ++i) name_view[i] = g.name_view[i];
      for (int i = 0; i < NAME_LEN; ++i) symbol_view[i] = g.symbol_view[i];
      for (int i = 0; i < DESC_LEN; ++i) desc_view[i] = g.desc_view[i];
    }
    return *this;
  }

  /// full species name
  std::string name() const { return std::string(name_view); }

  /// abbreviated name
  std::string symbol() const { return std::string(symbol_view); }

  /// species description
  std::string description() const { return std::string(desc_view); }

  /// Molecular weight [kg/mol]
  Real molecular_weight;

 private:
  char name_view[NAME_LEN];
  char symbol_view[NAME_LEN];
  char desc_view[DESC_LEN];
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
      "dichlorodifluoromethane"};
  const std::vector<std::string> gsymbs = {"O3",  "H2O2", "H2SO4", "SO2",
                                           "DMS", "SOAG", "O2",    "CO2",
                                           "N2O", "CH4",  "CFC11", "CFC12"};
  const std::vector<Real> gmws = {47.9982, 34.0136, 98.0784, 64.0648,
                                  62.1324, 12.011,  31.988,  44.009,
                                  44.013,  16.04,   137.37,  120.91};
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
