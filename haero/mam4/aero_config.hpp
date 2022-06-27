#ifndef HAERO_MAM4_AERO_CONFIG_HPP
#define HAERO_MAM4_AERO_CONFIG_HPP

#include <haero/mam4/aero_mode.hpp>

#include <haero/view_pack_helpers.hpp>

#include <algorithm>
#include <map>
#include <numeric>
#include <string>

namespace haero {
namespace mam4 {

/// MAM4 column-wise prognostic aerosol fields (also used for tendencies).
struct Prognostics final {
  Prognostics() = default;

  /// modal aerosol number mixing ratios
  ColumnView n_mode[4];

  /// aerosol mass mixing ratios within each mode
  ColumnView q_aero[4][7];

  /// gas mass mixing ratios
  ColumnView q_gas[13];
};

// Tendencies are identical in structure to prognostics.
using Tendencies = Prognostics;

/// MAM4 column-wise diagnostic aerosol fields.
struct Diagnostics final {
  Diagnostics() = default;
};

/// @struct MAM4::AeroConfig: for use with all MAM4 process implementations
struct AeroConfig final {
 public:

  // Types.
  using Prognostics = ::haero::mam4::Prognostics;
  using Diagnostics = ::haero::mam4::Diagnostics;
  using Tendencies  = ::haero::mam4::Tendencies;

  // Default constructor.
  AeroConfig() {}

  // Copy constructor.
  AeroConfig(const AeroConfig&) = default;

  // Destructor.
  ~AeroConfig() = default;

  // Assignment operator.
  AeroConfig& operator=(const AeroConfig&) = default;

  // Comparison operators.
  inline bool operator==(const AeroConfig& other) const {
    return true; // all MAM4 configs are equivalent
  }
  inline bool operator!=(const AeroConfig& other) const {
    return false; // all MAM4 configs are equivalent
  }

  /// Returns the number of aerosol modes.
  static int num_modes() { return 4; }

};

} // namespace mam4
} // namespace haero

#endif
