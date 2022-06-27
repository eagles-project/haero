#ifndef HAERO_MAM4_AERO_CONFIG_HPP
#define HAERO_MAM4_AERO_CONFIG_HPP

#include <haero/mam4/aero_modes.hpp>
#include <haero/view_pack_helpers.hpp>

#include <algorithm>
#include <map>
#include <numeric>
#include <string>

namespace haero {
namespace mam4 {

/// MAM4 column-wise prognostic aerosol fields (also used for tendencies).
class Prognostics final {
 public:
  explicit Prognostics(int num_levels): nlev_(num_levels) {
    for (int mode = 0; mode < 4; ++mode) {
      n_mode[mode] = ColumnView("n_mode", num_levels);
      for (int spec = 0; spec < 7; ++spec) {
        q_aero[mode][spec] = ColumnView("q_aero", num_levels);
      }
    }
    for (int gas = 0; gas < 13; ++gas) {
      q_gas[gas] = ColumnView("q_gas", num_levels);
    }
  }

  Prognostics(const Prognostics&) = default;
  ~Prognostics() = default;

  /// modal aerosol number mixing ratios (see aero_mode.hpp for indexing)
  ColumnView n_mode[4];

  /// aerosol mass mixing ratios within each mode (see aero_mode.hpp for
  /// indexing)
  ColumnView q_aero[4][7];

  /// gas mass mixing ratios (see aero_mode.hpp for indexing)
  ColumnView q_gas[13];

  int num_levels() const { return nlev_; }

  Prognostics() = delete;
  Prognostics& operator=(const Prognostics&) = delete;

 private:
  int nlev_;
};

// Tendencies are identical in structure to prognostics.
using Tendencies = Prognostics;

/// MAM4 column-wise diagnostic aerosol fields.
class Diagnostics final {
 public:
  explicit Diagnostics(int num_levels): nlev_(num_levels) {}
  Diagnostics(const Diagnostics&) = default;
  ~Diagnostics() = default;

  int num_levels() const { return nlev_; }

  Diagnostics() = delete;
  Diagnostics& operator=(const Diagnostics&) = delete;

 private:
  int nlev_;
};

/// @struct MAM4::AeroConfig: for use with all MAM4 process implementations
class AeroConfig final {
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
