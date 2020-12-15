#ifndef HAERO_SPECIES_HPP
#define HAERO_SPECIES_HPP

#include <string>
#include <vector>
#include <map>

namespace haero {

/// @struct Species
/// This type represents an aerosol or gas species.
struct Species final {

  /// Creates a new (aerosol or gas) species.
  /// @param [in] name A unique descriptive name for this species.
  /// @param [in] symbol A unique short, symbolic name for this species.
  Species(const std::string& name,
          const std::string& symbol):
    name(name), symbol(symbol) {}

  /// Full species name.
  std::string name;

  /// Abbreviated symbolic name.
  std::string symbol;
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
    Species("sulfate", "SO4"),
    Species("primary organic aerosol", "POA"),
    Species("secondary organic aerosol", "SOA"),
    Species("black carbon", "BC"),
    Species("sea salt", "SS"),
    Species("dust", "DST"),
    Species("marine organic aerosol", "MOA"),
  });
}

/// This factory function constructs a set of gas species corresponding to
/// the legacy MAM4 model. Includes:
/// 1. sulfuric acid (H2SO4)
/// 2. semi-volatile organic gas-phase species (SOAG)
inline std::vector<Species> create_mam4_gas_species() {
  return std::vector<Species>({
    Species("sulfuric acid", "H2SO4"),
    Species("semi-volatile organic gas-phase species", "SOAG"),
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
