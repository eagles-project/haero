#ifndef HAERO_SPECIES_HPP
#define HAERO_SPECIES_HPP

#include <string>

namespace haero {

/// This type represents an aerosol or gas species.
class Species final {
  public:
  /// Full species name.
  std::string name;
  /// Abbreviated symbolic name.
  std::string symbol;
};

}

#endif
