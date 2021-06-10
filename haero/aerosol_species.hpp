#ifndef HAERO_AEROSOL_SPECIES_HPP
#define HAERO_AEROSOL_SPECIES_HPP

#include <limits>
#include <map>
#include <string>
#include <vector>

#include "haero/haero.hpp"

namespace haero {

/// @struct AerosolSpecies
/// This type represents an aerosol species.
struct AerosolSpecies final {
  static const int NAME_LEN = 128;
  static const int DESC_LEN = 512;
  // Default constructor needed to resize Kokkos Views on device before deep
  // copy.
  KOKKOS_INLINE_FUNCTION
  AerosolSpecies() {
    name_view[0] = '\0';
    symbol_view[0] = '\0';
    desc_view[0] = '\0';
  }

  /// Creates a new aerosol species.
  /// @param [in] name A unique descriptive name for this species.
  /// @param [in] symbol A unique short, symbolic name for this species.
  /// @param [in] description A short text description of this species and its
  /// composition
  /// @param [in] molecular_wt The molecular weight [kg/mol]of the species
  /// @param [in] dry_rad The dry radius [m] of the species' particle size
  /// @param [in] hygro Base hygroscopicity of the species
  AerosolSpecies(const std::string& name, const std::string& symbol,
                 const std::string& description, Real molecular_wt, Real dens,
                 Real hygro)
      : molecular_weight(molecular_wt), density(dens), hygroscopicity(hygro) {
    EKAT_ASSERT(name.size() < NAME_LEN);
    EKAT_ASSERT(symbol.size() < NAME_LEN);
    strncpy(name_view, name.c_str(), NAME_LEN);
    strncpy(symbol_view, symbol.c_str(), NAME_LEN);
    strncpy(desc_view, description.c_str(), DESC_LEN);
  }

  KOKKOS_INLINE_FUNCTION
  AerosolSpecies(const AerosolSpecies& a)
      : molecular_weight(a.molecular_weight),
        density(a.density),
        hygroscopicity(a.hygroscopicity) {
    for (int i = 0; i < NAME_LEN; ++i) name_view[i] = a.name_view[i];
    for (int i = 0; i < NAME_LEN; ++i) symbol_view[i] = a.symbol_view[i];
    for (int i = 0; i < DESC_LEN; ++i) desc_view[i] = a.desc_view[i];
  }

  KOKKOS_INLINE_FUNCTION
  AerosolSpecies& operator=(const AerosolSpecies& a) {
    if (&a != this) {
      molecular_weight = a.molecular_weight;
      density = a.density;
      hygroscopicity = a.hygroscopicity;
      for (int i = 0; i < NAME_LEN; ++i) name_view[i] = a.name_view[i];
      for (int i = 0; i < NAME_LEN; ++i) symbol_view[i] = a.symbol_view[i];
      for (int i = 0; i < DESC_LEN; ++i) desc_view[i] = a.desc_view[i];
    }
    return *this;
  }

  /// Full species name.
  std::string name() const { return std::string(name_view); }

  /// Abbreviated symbolic name.
  std::string symbol() const { return std::string(symbol_view); }

  /// Species description.
  std::string description() const { return std::string(desc_view); }

  // Molecular weight [kg/mol]
  Real molecular_weight;

  /// Material density [kg/m^3]
  Real density;

  /// Hygroscopicity
  Real hygroscopicity;

 private:
  char name_view[NAME_LEN];
  char symbol_view[NAME_LEN];
  char desc_view[DESC_LEN];
};

/// This factory function constructs a set of aerosol species corresponding to
/// the legacy MAM4 model.
///
/// for info on data sources, see
/// https://eagles-project.atlassian.net/wiki/spaces/Computation/pages/1125515265/Aerosol+species+and+mode+data
inline std::vector<AerosolSpecies> create_mam4_aerosol_species() {
  const std::vector<std::string> aer_names = {"sulfate",
                                              "primary organic matter",
                                              "secondary organic aerosol",
                                              "black carbon",
                                              "dust",
                                              "sodium chloride",
                                              "marine organic matter"};
  const std::vector<std::string> aer_symbs = {"SO4", "POM",  "SOA", "BC",
                                              "DST", "NaCl", "MOM"};
  const std::vector<Real> aer_mw = {96.060,  12.011,  12.011, 12.011,
                                    135.064, 58.4425, 250093};
  const std::vector<Real> aer_dens = {1770, 1000, 1000, 1700, 2600, 1900, 1601};
  const std::vector<Real> aer_hygro = {0.507, 1e-10, 0.14, 1e-10,
                                       0.14,  1.16,  0.1};

  std::vector<AerosolSpecies> result;
  static constexpr Real g_to_kg =
      0.001;  /// Convert molecular weights to SI units (g/mol) -> (kg/mol)
  for (int i = 0; i < aer_mw.size(); ++i) {
    result.push_back(AerosolSpecies(aer_names[i], aer_symbs[i],
                                    "(No description)", g_to_kg * aer_mw[i],
                                    aer_dens[i], aer_hygro[i]));
  }
  return result;
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
inline std::map<std::string, std::vector<std::string> >
create_mam4_mode_species() {
  return std::map<std::string, std::vector<std::string> >({
      {"aitken", {"SO4", "NaCl", "SOA", "MOM"}},
      {"accumulation", {"SO4", "NaCl", "DST", "BC", "SOA", "POM", "MOM"}},
      {"coarse", {"SO4", "NaCl", "DST", "BC", "SOA", "POM", "MOM"}},
      {"primary_carbon", {"BC", "POM", "MOM"}},
  });
}

}  // namespace haero

#endif
