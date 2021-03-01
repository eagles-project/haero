#ifndef HAERO_DIAGNOSTICS_HPP
#define HAERO_DIAGNOSTICS_HPP

#include "haero/haero_config.hpp"
#include "haero/view_pack_helpers.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "kokkos/Kokkos_Core.hpp"
#include <map>
#include <vector>

namespace haero {

/// @class Diagnostics
/// This type stores a set of named diagnostic variables an aerosol system.
/// The set of diagnostic variables for such a system is determined by the
/// parameterizations selected for that system.
class Diagnostics final {
  public:

  /// Diagnostic variables are retrieved with at TOKEN.  A token is returned
  /// when the field is created or can be retrieved given the name of a field.
  /// The token is useful when retrieving a field on device since std::string
  /// is not usable on device.  A TOKEN can be obtained on host and then passed
  /// to device.
  typedef int TOKEN;
  static const TOKEN NOT_FOUND = -1;

  /// This type is used to define Kokkos views on device which might
  /// be either the host or the device if there is one.
  using kokkos_device_type = ekat::KokkosTypes<ekat::DefaultDevice>;

  /// This type represents an array mapping a vertical level index to a pack.
  /// The vertical level(s) are identified by the index.
  using ColumnView = kokkos_device_type::view_1d<PackType>;

  /// This type represents a multidimensional array mapping a species and
  /// vertical level index to a pack.
  /// * The species is identified by the index s.
  /// * The vertical level index is identified by the index k.
  /// So view[s][k] yields the desired pack.
  using SpeciesColumnView = kokkos_device_type::view_2d<PackType>;

  /// This type represents a multidimensional array mapping a mode and a
  /// vertical level index to a pack.
  /// * The mode is identified by the index m.
  /// * The vertical level index is identified by the index k.
  /// So view[m][k] yields the desired pack.
  using ModalColumnView = kokkos_device_type::view_2d<PackType>;

  /// Creates an empty Diagnostics to which data can be added.
  /// @param [in] num_aerosol_modes The number of aerosol modes in the system
  /// @param [in] num_aerosol_species A vector of length num_aerosol_modes whose
  ///                                 ith entry is the number of aerosol species
  ///                                 in the ith mode
  /// @param [in] num_gases The number of gas species in the atmosphere
  /// @param [in] num_levels The number of vertical levels per column stored by
  ///                        the state
  Diagnostics(int num_aerosol_modes,
              const std::vector<int>& num_aerosol_species,
              int num_gases,
              int num_levels);

  /// Destructor.
  KOKKOS_FUNCTION
  ~Diagnostics();

  // --------------------------------------------------------------------------
  //                                Metadata
  // --------------------------------------------------------------------------

  /// Returns the number of aerosol modes in the system.
  int num_aerosol_modes() const;

  /// Returns the number of aerosol species in the mode with the given index.
  /// @param [in] mode_index The index of the desired mode.
  int num_aerosol_species(int mode_index) const;

  /// Returns the total number of distinct aerosol species populations in the
  /// model, counting appearances of one species in different modes separately.
  int num_aerosol_populations() const;

  /// Returns the number of gas species in the system.
  int num_gases() const;

  /// Returns the number of vertical levels per column in the system.
  int num_levels() const;

  // --------------------------------------------------------------------------
  //                                  Data
  // --------------------------------------------------------------------------

  /// Returns token if the given (non-modal) variable exists within this object,
  /// NOT_FOUND otherwise.
  /// @param [in] name The name of the diagnostic variable of interest.
  TOKEN find_var(const std::string& name) const;

  /// Creates a diagnostic variable with the given name within this object.
  /// @param [in] name The name of the diagnostic variable to be created.
  TOKEN create_var(const std::string& name);

  /// Returns the view storing the diagnostic variable with the given name.
  /// If the variable does not yet exist, this throws an exception.
  /// @param [in] name The name of the diagnostic variable.
  KOKKOS_FUNCTION
  ColumnView var(const TOKEN token);

  /// Returns a const view storing the diagnostic variable with the given name.
  /// If the variable does not yet exist, this throws an exception.
  /// @param [in] name The name of the diagnostic variable.
  KOKKOS_INLINE_FUNCTION
  const ColumnView var(const TOKEN token) const {
    EKAT_KERNEL_REQUIRE_MSG(token < vars_.extent(0),
      "Diagnostic variable token not found!");
    const ColumnView vars = Kokkos::subview(vars_, token, Kokkos::ALL);
    return vars;
  }

  /// Returns token if the given modal aerosol variable exists within this object,
  /// NOT_FOUND otherwise.
  /// @param [in] name The name of the diagnostic variable of interest.
  TOKEN find_aerosol_var(const std::string& name) const;

  /// Creates a diagnostic modal aerosol variable with the given name within this
  /// object.
  /// @param [in] name The name of the modal diagnostic variable to be created.
  TOKEN create_aerosol_var(const std::string& name);

  /// Returns the view storing the modal aerosol diagnostic variable with the
  /// given name and mode index. If the variable does not yet exist, this throws
  /// an exception.
  /// @param [in] name The name of the diagnostic variable.
  KOKKOS_FUNCTION
  SpeciesColumnView aerosol_var(const TOKEN token);

  /// Returns a const view storing the modal aerosol diagnostic variable with
  /// the given name and mode index. If the variable does not yet exist, this
  /// throws an exception.
  /// @param [in] name The name of the diagnostic variable.
  KOKKOS_INLINE_FUNCTION
  const SpeciesColumnView aerosol_var(const TOKEN token) const {
    EKAT_KERNEL_REQUIRE_MSG(token < aero_vars_.extent(0),
      "Aerosol diagnostic variable token  not found!");
    const SpeciesColumnView vars = Kokkos::subview(aero_vars_, token, Kokkos::ALL, Kokkos::ALL);
    return vars;
  }

  /// Returns token if a variable defined for each gas species exists within
  /// this object with the given name, NOT_FOUND otherwise.
  /// @param [in] name The name of the diagnostic variable of interest.
  TOKEN find_gas_var(const std::string& name) const;

  /// Creates a diagnostic gas species variable with the given name within this
  /// object.
  /// @param [in] name The name of the modal diagnostic variable to be created.
  TOKEN create_gas_var(const std::string& name);

  /// Returns the view storing a diagnostic variable, defined for each gas
  /// species, with the given name. If the variable does not yet exist, this
  /// throws an exception.
  /// @param [in] name The name of the diagnostic variable.
  KOKKOS_FUNCTION
  SpeciesColumnView gas_var(const TOKEN token);

  /// Returns the const view storing a diagnostic variable, defined for each gas
  /// species, with the given name. If the variable does not yet exist, this
  /// throws an exception.
  /// @param [in] name The name of the diagnostic variable.
  KOKKOS_FUNCTION
  const SpeciesColumnView gas_var(const TOKEN token) const;

  /// Returns token if the given modal variable exists within this object,
  /// NOT_FOUND otherwise.
  /// @param [in] name The name of the diagnostic variable of interest.
  TOKEN find_modal_var(const std::string& name) const;

  /// Creates a diagnostic modal variable with the given name within this object.
  /// @param [in] name The name of the modal diagnostic variable to be created.
  TOKEN create_modal_var(const std::string& name);

  /// Returns the view storing the mode-specific diagnostic variable with the
  /// given name. If the variable does not yet exist, this throws an exception.
  /// @param [in] name The name of the diagnostic variable.
  KOKKOS_FUNCTION
  ModalColumnView modal_var(const TOKEN token);

  /// Returns a const view storing the mode-specific diagnostic variable with
  /// the given name. If the variable does not yet exist, this throws an
  /// exception.
  /// @param [in] name The name of the diagnostic variable.
  KOKKOS_FUNCTION
  const ModalColumnView modal_var(const TOKEN token) const;

  private:
  // Views that store arrays of views
  using ColumnViewArray        = kokkos_device_type::view_2d<PackType>;
  using SpeciesColumnViewArray = kokkos_device_type::view_3d<PackType>;
  using ModalColumnViewArray   = kokkos_device_type::view_3d<PackType>;

  // Set named string into map and return corresponding token.
  static TOKEN
  set_string_to_token(std::map<std::string,TOKEN> &registered_strings,
                      const std::string &name,
                      const TOKEN token);

  // Given named string search map and return corresponding token.
  static TOKEN
  get_string_to_token(const std::map<std::string,TOKEN> &registered_strings,
                      const std::string &name);

  // Functions that call the two functions above with the correct map.
  TOKEN set_string_to_token_vars (const std::string &name, const TOKEN token) ;
  TOKEN set_string_to_token_aero (const std::string &name, const TOKEN token) ;
  TOKEN set_string_to_token_gas  (const std::string &name, const TOKEN token) ;
  TOKEN set_string_to_token_modal(const std::string &name, const TOKEN token) ;
  TOKEN get_string_to_token_vars (const std::string &name) const;
  TOKEN get_string_to_token_aero (const std::string &name) const;
  TOKEN get_string_to_token_gas  (const std::string &name) const;
  TOKEN get_string_to_token_modal(const std::string &name) const;

  // Maps of Diagnostic variable names to the assigned tokens which are just
  // indexes into the array of views.
//  std::map<std::string,TOKEN> registered_strings_vars;
//  std::map<std::string,TOKEN> registered_strings_aero;
//  std::map<std::string,TOKEN> registered_strings_gas;
//  std::map<std::string,TOKEN> registered_strings_modal;

  // Number of aerosol species in each mode.
  const view_1d_int_type num_aero_species_;

  // Number of distinct aerosol populations.
  int num_aero_populations_;

  // Number of gas species.
  int num_gases_;

  /// Number of vertical levels per column in the system.
  const int num_levels_;

  // Named diagnostic variables.  These are arrays of views in which the
  // assigned token can be used to index to the proper sub-view.
  ColumnViewArray        vars_;
  SpeciesColumnViewArray aero_vars_;
  SpeciesColumnViewArray gas_vars_;
  ModalColumnViewArray   modal_vars_;
};

}

#endif
