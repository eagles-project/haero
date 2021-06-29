#ifndef HAERO_REGION_OF_VALIDITY_HPP
#define HAERO_REGION_OF_VALIDITY_HPP

#include <map>
#include <vector>

#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "haero/atmosphere.hpp"
#include "haero/haero.hpp"
#include "haero/prognostics.hpp"
#include "haero/view_pack_helpers.hpp"
#include "kokkos/Kokkos_Core.hpp"

namespace haero {

/// @class RegionOfValidity
/// This type expresses the "region of validity" for a given aerosol process
/// in terms of bounds on atmospheric and aerosol prognostic states.
class RegionOfValidity {
 public:
  using Bounds = std::pair<Real, Real>;

  /// Gas- and aerosol-specific bounds are identified by unique tokens. Each
  /// token is generated when a set of bounds is added, and can be retrieved
  /// given the symbolic name of the gas or aerosol. See the Diagnostics
  /// container for a more complete explanation of how tokens allow us to
  /// fake the std::map<std::string, ...> type on GPUs.
  typedef int Token;

  /// This token indicates that a requested set of bounds was not found.
  static const Token BOUNDS_NOT_FOUND = -1;

  /// Constructor
  RegionOfValidity();

  /// Destructor.
  KOKKOS_FUNCTION
  ~RegionOfValidity();

  /// Minimum and maximum bounds on atmospheric temperature [K]
  Bounds temp_bounds;
  /// Minimum and maximum bounds on relative humidity [-]
  Bounds rel_hum_bounds;

  /// Minimum and maximum bounds on the mass mixing ratio for the gas species
  /// corresponding to the given token.
  KOKKOS_INLINE_FUNCTION
  const Bounds &gas_bounds(const Token token) const {
    EKAT_KERNEL_REQUIRE_MSG(token < gas_bounds_.extent(0),
                            "Gas species token not found!");
    return gas_bounds_(token);
  }

  /// Minimum and maximum bounds on the mass mixing ratio for the gas species
  /// corresponding to the given token (non-const).
  KOKKOS_INLINE_FUNCTION
  Bounds &gas_bounds(const Token token) {
    EKAT_KERNEL_REQUIRE_MSG(token < gas_bounds_.extent(0),
                            "Gas species token not found!");
    return gas_bounds_(token);
  }

  /// Returns true if all state data in the given atmosphere container falls
  /// within this region of validity, false otherwise.
  KOKKOS_INLINE_FUNCTION
  bool contains(const Atmosphere &atmosphere) const { return true; }

  /// Returns true if all state data in the given prognostics container falls
  /// within this region of validity, false otherwise.
  KOKKOS_INLINE_FUNCTION
  bool contains(const Prognostics &prognostics) const { return true; }

 protected:
  /// Minimum and maximum bounds on specific gas species, indexed by (case-
  /// insensitive) symbols.
  using BoundsArray = kokkos_device_type::view_1d<Bounds>;
  BoundsArray gas_bounds_;
};

/// @class HostRegionOfValidity
/// This type allows the construction of a RegionOfValidity on a host using
/// simple types that don't work on GPUs.
class HostRegionOfValidity final : public RegionOfValidity {
 public:
  /// Creates an empty region of validity for which bounds can be specified.
  HostRegionOfValidity();

  /// Destructor.
  ~HostRegionOfValidity();

  /// Returns a RegionOfValidity object in which all bounds can be accessed
  /// with Tokens. This allows a region of validity to be copied to a device.
  const RegionOfValidity &getRegionOfValidity() const;

  /// Returns a unique token that identifies bounds on a gas species, given its
  /// symbol. Returns BOUNDS_NOT_FOUND if such bounds do not exist.
  /// @param [in] symbol The symbolic name of the gas species of interest
  Token find_gas_bounds(const std::string &symbol) const;

  /// Adds a set of bounds for the gas species with the given name, returning
  /// a unique token for the new bounds.
  /// @param [in] symbol The symbolic name of the gas species of interest
  Token add_gas_bounds(const std::string &symbol);

 private:
  using BoundsArray = kokkos_device_type::view_1d<Bounds>;

  // Set named string into map and return corresponding token.
  static Token set_string_to_token(
      std::map<std::string, Token> &registered_strings, const std::string &name,
      const Token token);

  // Given named string search map and return corresponding token.
  static Token get_string_to_token(
      const std::map<std::string, Token> &registered_strings,
      const std::string &name);

  // Functions that call the two functions above with the correct map.
  Token set_string_to_token_gases(const std::string &name, const Token token);
  Token get_string_to_token_gases(const std::string &name) const;

  // Maps of Diagnostic variable names to the assigned tokens which are just
  // indexes into the array of views.
  std::map<std::string, Token> registered_gases_;
};

}  // namespace haero

#endif
