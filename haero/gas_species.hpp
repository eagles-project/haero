#ifndef HAERO_GAS_SPECIES_HPP
#define HAERO_GAS_SPECIES_HPP

#include "haero/haero.hpp"

#include <string>
#include <vector>
#include <map>
#include <limits>

namespace haero {

/// @struct GasSpecies
/// This type represents a gas that participates in one or more aerosol microphysics parameterizations.
struct GasSpecies final {

  // Default constructor needed for device
  KOKKOS_INLINE_FUNCTION
  GasSpecies() : name_view(), symbol_view() {}

  /// Creates a new gas species
  /** @param [in] name A unique descriptive name for this species
      @param [in] symbol a unique short name or symbol for this species
      @param [in] molecular_wt the molecular weight of this species [kg/mol]
  */
  GasSpecies(const std::string& name,
    const std::string& symbol,
    Real molecular_wt) :
    molecular_weight(molecular_wt),
    name_view(name), symbol_view(symbol) {}


  /// full species name
  std::string name() const {return name_view.label();}

  /// abbreviated name
  std::string symbol() const {return symbol_view.label();}

  /// Molecular weight [kg/mol]
  Real molecular_weight;

  private:
    Kokkos::View<int> name_view;
    Kokkos::View<int> symbol_view;
};

/// This factory function constructs a set of gas species corresponding to
/// the legacy MAM4 model.
///
/// for info on data sources, see
/// https://eagles-project.atlassian.net/wiki/spaces/Computation/pages/1125515265/Aerosol+species+and+mode+data
inline std::vector<GasSpecies> create_mam4_gas_species() {
  const std::vector<std::string> gnames =
    {"ozone", "hydrogen peroxide", "sulfuric acid", "sulfur dioxide", "dimethylsulfide", "secondary organic aerosol precursor",
    "oxygen", "carbon dioxide", "nitrous oxide", "methane", "trichlorofluoromethane", "dichlorodifluoromethane"};
  const std::vector<std::string> gsymbs = {"O3", "H2O2", "H2SO4", "SO2", "DMS", "SOAG",
    "O2", "CO2", "N2O", "CH4", "CFC11", "CFC12"};
  const std::vector<Real> gmws = {47.9982, 34.0136, 98.0784, 64.0648, 62.1324, 12.011,
    31.988, 44.009, 44.013, 16.04, 137.37, 120.91};
  std::vector<GasSpecies> result;
  for (int i=0; i<gnames.size(); ++i) {
    result.push_back(GasSpecies(gnames[i], gsymbs[i], gmws[i]));
  }
  return result;
}

}
#endif
