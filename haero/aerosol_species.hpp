#ifndef HAERO_AEROSOL_SPECIES_HPP
#define HAERO_AEROSOL_SPECIES_HPP

#include "haero/haero.hpp"

#include <string>
#include <vector>
#include <map>
#include <limits>

namespace haero {

/// @struct AerosolSpecies
/// This type represents an aerosol species.
struct AerosolSpecies final {

  // Default constructor needed to resize Kokkos Views on device before deep copy.
  KOKKOS_INLINE_FUNCTION
  AerosolSpecies(): name_view(), symbol_view() {}

  /// Creates a new aerosol species.
  /// @param [in] name A unique descriptive name for this species.
  /// @param [in] symbol A unique short, symbolic name for this species.
  /// @param [in] molecular_wt The molecular weight [kg/mol]of the species
  /// @param [in] dry_rad The dry radius [m] of the species' particle size
  /// @param [in] hygro Base hygroscopicity of the species
  AerosolSpecies(const std::string& name,
          const std::string& symbol,
          Real molecular_wt,
          Real dry_rad,
          Real dens,
          Real hygro):
    molecular_weight(molecular_wt), dry_radius(dry_rad),
    density(dens), hygroscopicity(hygro), name_view(name), symbol_view(symbol) {}

  /// Full species name.
  std::string name() const { return name_view.label(); }

  /// Abbreviated symbolic name.
  std::string symbol() const { return symbol_view.label(); }

  // Molecular weight [kg/mol]
  Real molecular_weight;

  /// Dry radius [m]
  Real dry_radius;

  /// Material density [kg/m^3]
  Real density;

  /// Hygroscopicity
  Real hygroscopicity;

private:
  Kokkos::View<int> name_view;
  Kokkos::View<int> symbol_view;
};

/// This factory function constructs a set of aerosol species corresponding to
/// the legacy MAM4 model.
///
/// for info on data sources, see
/// https://eagles-project.atlassian.net/wiki/spaces/Computation/pages/1125515265/Aerosol+species+and+mode+data
inline std::vector<AerosolSpecies> create_mam4_aerosol_species() {
  const std::vector<std::string> aer_names = {"sulfate", "primary organic matter",
    "secondary organic aerosol", "black carbon", "dust", "sodium chloride", "marine organic matter"};
  const std::vector<std::string> aer_symbs = {"SO4",  "POM",  "SOA",   "BC",  "DST",   "NaCl",  "MOM"};
  const std::vector<Real> aer_mw =           {115.107, 12.011,  12.011,  12.011,  135.064, 58.4425, 250093};
  const std::vector<Real> aer_dry_rad =      {6.95e-8, 2.12e-8, 2.12e-8, 1.18e-8, 1.51e-6, 2.09e-7, 2.09e-7};
  const std::vector<Real> aer_dens =         {1770,    1000,    1000,    1700,    2600,    1900,    1601};
  const std::vector<Real> aer_hygro =        {0.507,   1e-10,   0.14,    1e-10,   0.14,    1.16,    0.1};

  std::vector<AerosolSpecies> result;
  static constexpr Real g_to_kg = 0.001; /// Convert molecular weights to SI units (g/mol) -> (kg/mol)
  for (int i=0; i<aer_mw.size(); ++i) {
    result.push_back(AerosolSpecies(aer_names[i], aer_symbs[i],
      g_to_kg*aer_mw[i], aer_dry_rad[i], aer_dens[i], aer_hygro[i]));
  }
  return result;
}

/// This factory function constructs a set of gas species corresponding to
/// the legacy MAM4 model. Includes:
/// 1. sulfuric acid (H2SO4)
/// 2. semi-volatile organic gas-phase species (SOAG)
inline std::vector<AerosolSpecies> create_mam4_gas_species() {
  return std::vector<AerosolSpecies>({
    // FIXME: All of the numerical parameters here are wrong and need to be fixed
    AerosolSpecies("sulfuric acid", "H2SO4", 1.0, 1.0, 1.0, 1.0),
    AerosolSpecies("semi-volatile organic gas-phase species", "SOAG", 1.0, 1.0, 1.0, 1.0),
  });
}

/// This factory function constructs a map of mode name -> species symbol
/// associations corresponding to those used in the legacy MAM4 model:
///
/// +----------------+---------------------------------
/// | Mode           ¦ Species
/// +----------------+---------------------------------
/// | aitken         | SO4, NaCl, SOA, MOM
/// +----------------+---------------------------------
/// | accumulation   | SO4, NaCl, DST, BC, SOA, POM, MOM
/// +----------------+---------------------------------
/// | coarse         | SO4, NaCl, DST, BC, SOA, POM, MOM
/// +----------------+---------------------------------
/// | primary_carbon | BC, POM, MOM
/// +----------------+---------------------------------
inline std::map<std::string, std::vector<std::string> > create_mam4_mode_species() {
  return std::map<std::string, std::vector<std::string> >({
    {"aitken", {"SO4", "NaCl", "SOA", "MOM"}},
    {"accumulation", {"SO4", "NaCl", "DST", "BC", "SOA", "POM", "MOM"}},
    {"coarse", {"SO4", "NaCl", "DST", "BC", "SOA", "POM", "MOM"}},
    {"primary_carbon", {"BC", "POM", "MOM"}},
  });
}


}

#endif
