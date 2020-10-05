#ifndef HAERO_SPECIES_HPP
#define HAERO_SPECIES_HPP

#include <string>

namespace haero {

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

}

#endif
