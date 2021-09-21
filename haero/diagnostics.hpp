#ifndef HAERO_DIAGNOSTICS_HPP
#define HAERO_DIAGNOSTICS_HPP

#include <map>
#include <vector>

#include "haero/haero.hpp"
#include "haero/view_pack_helpers.hpp"
#include "kokkos/Kokkos_Core.hpp"

namespace haero {

/// @class Diagnostics
/// This type stores a set of named diagnostic variables an aerosol system.
/// The set of diagnostic variables for such a system is determined by the
/// parameterizations selected for that system.
class Diagnostics {
 public:
  /// Diagnostic variables are identified by a unique token. This token is
  /// generated when a field is created, and can be retrieved given the name of
  /// an existing field. Tokens are useful for retrieving fields on a device
  /// (as opposed to a host), since std::string is not usable on GPUs and other
  /// accelerators. Tokens, by contracts, can be obtained on a host and then
  /// passed to a device.
  typedef int Token;

  /// This token indicates that a requested variable was not found.
  static const Token VAR_NOT_FOUND = -1;

  /// Creates an empty Diagnostics to which data can be added.
  /// @param [in] num_aerosol_modes The number of aerosol modes in the system
  /// @param [in] num_aerosol_species A vector of length num_aerosol_modes whose
  ///                                 ith entry is the number of aerosol species
  ///                                 in the ith mode
  /// @param [in] num_gases The number of gas species in the atmosphere
  /// @param [in] num_levels The number of vertical levels per column stored by
  ///                        the state
  Diagnostics(int num_aerosol_modes,
              const std::vector<int> &num_aerosol_species, int num_gases,
              int num_levels);

  /// Default constructor and copy constructor is needed to define Views
  /// of Diagnostics classes and then copy data into the views.
  Diagnostics() = default;
  Diagnostics(const Diagnostics &d) = default;

  /// Destructor.
  KOKKOS_INLINE_FUNCTION
  ~Diagnostics() {}

  // --------------------------------------------------------------------------
  //                                Metadata
  // --------------------------------------------------------------------------

  /// Returns the number of aerosol modes in the system.
  KOKKOS_FUNCTION
  int num_aerosol_modes() const;

  /// Returns the number of aerosol species in the mode with the given index.
  /// @param [in] mode_index The index of the desired mode.
  KOKKOS_IMPL_DEVICE_FUNCTION
  int num_aerosol_species(int mode_index) const;

  /// Returns the total number of distinct aerosol species populations in the
  /// model, counting appearances of one species in different modes separately.
  KOKKOS_FUNCTION
  int num_aerosol_populations() const;

  /// Returns the number of gas species in the system.
  KOKKOS_FUNCTION
  int num_gases() const;

  /// Returns the number of vertical levels per column in the system.
  KOKKOS_FUNCTION
  int num_levels() const { return num_levels_; }

  // --------------------------------------------------------------------------
  //                                  Data
  // --------------------------------------------------------------------------

  /// Returns a const view storing the diagnostic variable with the given name.
  /// If the variable does not yet exist, this throws an exception.
  /// Returns a const view storing the diagnostic variable with a name
  /// corresponding to the given token. If such a variable does not exist, this
  /// throws an exception.
  /// @param [in] token A unique token identifying a diagnostic variable.
  KOKKOS_INLINE_FUNCTION
  ColumnView var(const Token token) const {
    EKAT_KERNEL_REQUIRE_MSG(token < vars_.extent(0),
                            "Diagnostic variable token not found!");
    const ColumnView vars = Kokkos::subview(vars_, token, Kokkos::ALL);
    return vars;
  }

  /// Returns a const view storing the modal aerosol diagnostic variable with a
  /// name corresponding to the given token. If such a variable does not exist,
  /// this throws an exception.
  /// @param [in] name The name of the diagnostic variable.
  KOKKOS_INLINE_FUNCTION
  SpeciesColumnView aerosol_var(const Token token) const {
    EKAT_KERNEL_REQUIRE_MSG(token < aero_vars_.extent(0),
                            "Aerosol diagnostic variable token not found!");
    const SpeciesColumnView vars =
        Kokkos::subview(aero_vars_, token, Kokkos::ALL, Kokkos::ALL);
    return vars;
  }

  /// Returns a const view storing the gas diagnostic variable with a name
  /// corresponding to the given token. If such a variable does not exist, this
  /// throws an exception.
  /// @param [in] token A unique token identifying a diagnostic variable.
  KOKKOS_INLINE_FUNCTION
  SpeciesColumnView gas_var(const Token token) const {
    EKAT_KERNEL_REQUIRE_MSG(token < gas_vars_.extent(0),
                            "Gas diagnostic variable token not found!");
    const SpeciesColumnView vars =
        Kokkos::subview(gas_vars_, token, Kokkos::ALL, Kokkos::ALL);
    return vars;
  }

  /// Returns a const view storing the mode-specific diagnostic variable with a
  /// name corresponding to the given token. If such a variable does not exist,
  /// this throws an exception.
  /// @param [in] token A unique token identifying a diagnostic variable.
  KOKKOS_FUNCTION
  ModeColumnView modal_var(const Token token) const;

 protected:
  // Views that store arrays of views
  using ColumnViewArray = kokkos_device_type::view_2d<PackType>;
  using SpeciesColumnViewArray = kokkos_device_type::view_3d<PackType>;
  using ModeColumnViewArray = kokkos_device_type::view_3d<PackType>;

  // Number of aerosol species in each mode.
  view_1d_int_type num_aero_species_;

  // Number of distinct aerosol populations.
  int num_aero_populations_;

  // Number of gas species.
  int num_gases_;

  /// Number of vertical levels per column in the system.
  int num_levels_;

  // Named diagnostic variables.  These are arrays of views in which the
  // assigned token can be used to index to the proper sub-view.
  ColumnViewArray vars_;
  SpeciesColumnViewArray aero_vars_;
  SpeciesColumnViewArray gas_vars_;
  ModeColumnViewArray modal_vars_;
};

/// @class HostDiagnostics
/// This type stores a set of named diagnostic variables an aerosol system.
/// The set of diagnostic variables for such a system is determined by the
/// parameterizations selected for that system.
class HostDiagnostics final : public Diagnostics {
 public:
  /// Creates an empty HostDiagnostics to which data can be added.
  /// @param [in] num_aerosol_modes The number of aerosol modes in the system
  /// @param [in] num_aerosol_species A vector of length num_aerosol_modes whose
  ///                                 ith entry is the number of aerosol species
  ///                                 in the ith mode
  /// @param [in] num_gases The number of gas species in the atmosphere
  /// @param [in] num_levels The number of vertical levels per column stored by
  ///                        the state
  HostDiagnostics(int num_aerosol_modes,
                  const std::vector<int> &num_aerosol_species, int num_gases,
                  int num_levels);

  /// Destructor.
  ~HostDiagnostics();

  // --------------------------------------------------------------------------
  //                                  Data
  // --------------------------------------------------------------------------

  /// Returns the Diagnostics class in which all of the diagnostic views can
  /// be accessed using Tokens. Since Token is a POD, it is usable on device.
  /// Therefore the Diagnostics class can be copied to device without error.
  const Diagnostics &GetDiagnostics();

  /// Returns the number of aerosol species in the mode with the given index.
  /// @param [in] mode_index The index of the desired mode.
  int num_aerosol_species(int mode_index) const;

  /// Returns a unique token that identifies the given (non-modal) variable
  /// within this object. Returns VAR_NOT_FOUND if this variable does not exist.
  /// @param [in] name The name of the diagnostic variable of interest.
  Token find_var(const std::string &name) const;

  /// Creates a diagnostic variable with the given name within this object,
  /// returning a unique token for the new variable.
  /// @param [in] name The name of the diagnostic variable to be created.
  Token create_var(const std::string &name);

  /// Returns a unique token that identifies the given modal aerosol variable
  /// within this object. Returns VAR_NOT_FOUND if this variable does not exist.
  /// @param [in] name The name of the diagnostic variable of interest.
  Token find_aerosol_var(const std::string &name) const;

  /// Creates a diagnostic modal aerosol variable with the given name within
  /// this object, returning a unique token for the new variable.
  /// @param [in] name The name of the modal diagnostic variable to be created.
  Token create_aerosol_var(const std::string &name);

  /// Returns a unique token that identifies the given gas variable within this
  /// object. Returns VAR_NOT_FOUND if this variable does not exist.
  /// @param [in] name The name of the diagnostic variable of interest.
  Token find_gas_var(const std::string &name) const;

  /// Creates a diagnostic gas species variable with the given name within this
  /// object, returning a unique token for the new variable.
  /// @param [in] name The name of the modal diagnostic variable to be created.
  Token create_gas_var(const std::string &name);

  /// Returns a unique token that identifies the given modal variable within
  /// this object. Returns VAR_NOT_FOUND if this variable does not exist.
  /// @param [in] name The name of the diagnostic variable of interest.
  Token find_modal_var(const std::string &name) const;

  /// Creates a diagnostic modal variable with the given name within this
  /// object, returning a unique token for the new variable.
  /// @param [in] name The name of the modal diagnostic variable to be created.
  Token create_modal_var(const std::string &name);

 private:
  // Set named string into map and return corresponding token.
  static Token set_string_to_token(
      std::map<std::string, Token> &registered_strings, const std::string &name,
      const Token token);

  // Given named string search map and return corresponding token.
  static Token get_string_to_token(
      const std::map<std::string, Token> &registered_strings,
      const std::string &name);

  void clear_maps();
  // Functions that call the two functions above with the correct map.
  Token set_string_to_token_vars(const std::string &name, const Token token);
  Token set_string_to_token_aero(const std::string &name, const Token token);
  Token set_string_to_token_gas(const std::string &name, const Token token);
  Token set_string_to_token_modal(const std::string &name, const Token token);
  Token get_string_to_token_vars(const std::string &name) const;
  Token get_string_to_token_aero(const std::string &name) const;
  Token get_string_to_token_gas(const std::string &name) const;
  Token get_string_to_token_modal(const std::string &name) const;

  // Maps of Diagnostic variable names to the assigned tokens which are just
  // indexes into the array of views.
  std::map<std::string, Token> registered_strings_vars;
  std::map<std::string, Token> registered_strings_aero;
  std::map<std::string, Token> registered_strings_gas;
  std::map<std::string, Token> registered_strings_modal;
};

}  // namespace haero

#endif
