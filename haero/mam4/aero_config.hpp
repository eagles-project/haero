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
  /// Creates a container for prognostic variables on the specified number of
  /// vertical levels.
  explicit Prognostics(int num_levels): nlev_(num_levels) {
    const int nk = PackInfo::num_packs(num_levels);
    for (int mode = 0; mode < 4; ++mode) {
      n_mode[mode] = ColumnView("n_mode", nk);
      for (int spec = 0; spec < 7; ++spec) {
        q_aero[mode][spec] = ColumnView("q_aero", nk);
      }
    }
    for (int gas = 0; gas < 13; ++gas) {
      q_gas[gas] = ColumnView("q_gas", nk);
    }
  }

  Prognostics() = default; // Careful! Only for creating placeholders in views
  Prognostics(const Prognostics&) = default;
  ~Prognostics() = default;
  Prognostics& operator=(const Prognostics&) = default;

  /// modal aerosol number mixing ratios (see aero_mode.hpp for indexing)
  ColumnView n_mode[4];

  /// aerosol mass mixing ratios within each mode (see aero_mode.hpp for
  /// indexing)
  ColumnView q_aero[4][7];

  /// gas mass mixing ratios (see aero_mode.hpp for indexing)
  ColumnView q_gas[13];

  /// For gas-aerosol exchange process (probably temporary)
  ColumnView uptkrate_h2so4;

  KOKKOS_INLINE_FUNCTION
  int num_levels() const { return nlev_; }

  /// Returns true iff all prognostic quantities are nonnegative, using the
  /// given thread team to parallelize the check.
  KOKKOS_INLINE_FUNCTION
  bool quantities_nonnegative(const TeamType& team) const {
    const int nk = PackInfo::num_packs(num_levels());
    int violations = 0;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, nk),
      KOKKOS_LAMBDA(int k, int& violation) {
        for (int mode = 0; mode < 4; ++mode) { // check mode mmrs
          if ((n_mode[mode](k) < 0).any()) {
            ++violation;
          } else {
            for (int spec = 0; spec < 7; ++spec) { // check aerosol mmrs
              if ((q_aero[mode][spec](k) < 0).any()) {
                ++violation;
                break;
              }
            }
          }
          if (violation > 0) break;
        }
        if (violation == 0) {
          for (int gas = 0; gas < 13; ++gas) { // check gas mmrs
            if ((q_gas[gas](k) < 0).any()) ++violation;
          }
        }
      }, violations);
    return (violations == 0);
  }

 private:
  int nlev_;
};

// Tendencies are identical in structure to prognostics.
using Tendencies = Prognostics;

/// MAM4 column-wise diagnostic aerosol fields.
class Diagnostics final {
 public:
  explicit Diagnostics(int num_levels): nlev_(num_levels) {}
  Diagnostics() = default; // Careful! Only for creating placeholders in views
  Diagnostics(const Diagnostics&) = default;
  ~Diagnostics() = default;
  Diagnostics& operator=(const Diagnostics&) = default;

  int num_levels() const { return nlev_; }

  ColumnView dry_geometric_mean_diameter;
  ColumnView wet_geometric_mean_diameter;
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

  bool calculate_gas_uptake_coefficient = false;
  int  number_gauss_points_for_integration = 2;
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