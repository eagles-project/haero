#ifndef HAERO_MAM4_AERO_CONFIG_HPP
#define HAERO_MAM4_AERO_CONFIG_HPP

#include "aero_mode.hpp"

#include <haero/aero_species.hpp>
#include <haero/gas_species.hpp>
#include <haero/view_pack_helpers.hpp>

#include <algorithm>
#include <map>
#include <numeric>
#include <string>

namespace haero {

namespace mam4 {

using ColumnView = view_2d_pack_type;

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

/// MAM4 column-wise diagnostic aerosol fields.
struct Diagnostics final {
  Diagnostics() = default;
};

/// @struct MAM4::AeroConfig: for use with all MAM4 process implementations
struct AeroConfig final {
 public:

  // Types.
  using Prognostics = Prognostics;
  using Diagnostics = Diagnostics;

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
  }
  inline bool operator!=(const AeroConfig& other) const {
    return (!(*this == other));
  }

  /// Returns the number of aerosol modes.
  static int num_modes() { return 4; }

};

}  // namespace haero

#endif
