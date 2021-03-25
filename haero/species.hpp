#ifndef HAERO_SPECIES_HPP
#define HAERO_SPECIES_HPP

#include "haero/haero.hpp"

#include <string>
#include <vector>
#include <map>

namespace haero {

/// @struct Species
/// This type represents an aerosol or gas species.
struct Species final {

  // Default constructor needed to resize Kokkos Views on device before deep copy.
  KOKKOS_INLINE_FUNCTION
  Species() { 
    p_name[0]='\0'; 
    p_symbol[0]='\0';
  }

  /// Creates a new (aerosol or gas) species.
  /// @param [in] name A unique descriptive name for this species.
  /// @param [in] symbol A unique short, symbolic name for this species.
  /// @param [in] molecular_wt The molecular weight of the species
  /// @param [in] crystal_pt The crystalization point of the species
  /// @param [in] deliq_pt The deliquenscence point of the species
  Species(const std::string& name,
          const std::string& symbol,
          Real molecular_wt,
          Real crystal_pt,
          Real deliq_pt):
    molecular_weight(molecular_wt), crystalization_point(crystal_pt),
    deliquescence_point(deliq_pt) {
    EKAT_ASSERT(name.size() < 100);
    EKAT_ASSERT(symbol.size() < 100);
    strncpy(p_name, name.c_str(), 100);
    strncpy(p_symbol, symbol.c_str(), 100);
  }

  /// Full species name.
  std::string name() const { return std::string(p_name); }

  /// Abbreviated symbolic name.
  std::string symbol() const { return std::string(p_symbol); }

  // Molecular weight [kg/mol]
  Real molecular_weight;

  // Crystalization point (relative humidity threshold) [-]
  Real crystalization_point;

  // Deliquescence point (relative humidity threshold) [-]
  Real deliquescence_point;
private:
  char p_name[100];
  char p_symbol[100];
};

/// This factory function constructs a set of aerosol species corresponding to
/// the legacy MAM4 model. Includes:
/// 1. sulfate (SO4)
/// 2. primary organic aerosol (POA)
/// 3. secondary organic aerosol (SOA)
/// 4. black carbon (BC)
/// 5. sea salt (SS)
/// 6. dust (DST)
/// 7. marine organic aerosol (MOA)
inline std::vector<Species> create_mam4_aerosol_species() {
  return std::vector<Species>({
    // FIXME: All of the numerical parameters here are wrong and need to be
    // FIXME: fixed (except perhaps the molecular weight of SO4).
    Species("sulfate", "SO4", 96., 1.0, 1.0),
    Species("primary organic aerosol", "POA", 1.0, 1.0, 1.0),
    Species("secondary organic aerosol", "SOA", 1.0, 1.0, 1.0),
    Species("black carbon", "BC", 1.0, 1.0, 1.0),
    Species("sea salt", "SS", 1.0, 1.0, 1.0),
    Species("dust", "DST", 1.0, 1.0, 1.0),
    Species("marine organic aerosol", "MOA", 1.0, 1.0, 1.0),
  });
}

/// This factory function constructs a set of gas species corresponding to
/// the legacy MAM4 model. Includes:
/// 1. sulfuric acid (H2SO4)
/// 2. semi-volatile organic gas-phase species (SOAG)
inline std::vector<Species> create_mam4_gas_species() {
  return std::vector<Species>({
    // FIXME: All of the numerical parameters here are wrong and need to be fixed
    Species("sulfuric acid", "H2SO4", 1.0, 1.0, 1.0),
    Species("semi-volatile organic gas-phase species", "SOAG", 1.0, 1.0, 1.0),
  });
}

/// This factory function constructs a map of mode name -> species symbol
/// associations corresponding to those used in the legacy MAM4 model:
///
/// +----------------+---------------------------------
/// | Mode           ¦ Species
/// +----------------+---------------------------------
/// | aitken         | SO4, SS, SOA, MOA
/// +----------------+---------------------------------
/// | accumulation   | SO4, SS, DST, BC, SOA, POA, MOA
/// +----------------+---------------------------------
/// | coarse         | SO4, SS, DST, BC, SOA, POA, MOA
/// +----------------+---------------------------------
/// | primary_carbon | BC, POA, MOA
/// +----------------+---------------------------------
inline std::map<std::string, std::vector<std::string> > create_mam4_mode_species() {
  return std::map<std::string, std::vector<std::string> >({
    {"aitken", {"SO4", "SS", "SOA", "MOA"}},
    {"accumulation", {"SO4", "SS", "DST", "BC", "SOA", "POA", "MOA"}},
    {"coarse", {"SO4", "SS", "DST", "BC", "SOA", "POA", "MOA"}},
    {"primary_carbon", {"BC", "POA", "MOA"}},
  });
}


}

#endif
