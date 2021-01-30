#ifndef HAERO_DIAGNOSTICS_HPP
#define HAERO_DIAGNOSTICS_HPP

#include "haero/haero_config.hpp"
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

  /// This type represents an array mapping a vertical level index to a pack.
  /// The vertical level(s) are identified by the index.
  using ColumnView = Kokkos::View<PackType*>;

  /// This type represents a multidimensional array mapping a species and
  /// vertical level index to a pack.
  /// * The species is identified by the index s.
  /// * The vertical level index is identified by the index k.
  /// So view[s][k] yields the desired pack.
  using SpeciesColumnView = Kokkos::View<PackType**>;

  /// This type represents a multidimensional array mapping a mode and a
  /// vertical level index to a pack.
  /// * The mode is identified by the index m.
  /// * The vertical level index is identified by the index k.
  /// So view[m][k] yields the desired pack.
  using ModalColumnView = Kokkos::View<PackType**>;

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

  /// Returns true if the given (non-modal) variable exists within this object,
  /// false otherwise.
  /// @param [in] name The name of the diagnostic variable of interest.
  bool has_var(const std::string& name) const;

  /// Creates a diagnostic variable with the given name within this object.
  /// @param [in] name The name of the diagnostic variable to be created.
  void create_var(const std::string& name);

  /// Returns the view storing the diagnostic variable with the given name.
  /// If the variable does not yet exist, this throws an exception.
  /// @param [in] name The name of the diagnostic variable.
  ColumnView& var(const std::string& name);

  /// Returns a const view storing the diagnostic variable with the given name.
  /// If the variable does not yet exist, this throws an exception.
  /// @param [in] name The name of the diagnostic variable.
  const ColumnView& var(const std::string& name) const;

  /// Returns true if the given modal aerosol variable exists within this object,
  /// false otherwise.
  /// @param [in] name The name of the diagnostic variable of interest.
  bool has_aerosol_var(const std::string& name) const;

  /// Creates a diagnostic modal aerosol variable with the given name within this
  /// object.
  /// @param [in] name The name of the modal diagnostic variable to be created.
  void create_aerosol_var(const std::string& name);

  /// Returns the view storing the modal aerosol diagnostic variable with the
  /// given name and mode index. If the variable does not yet exist, this throws
  /// an exception.
  /// @param [in] name The name of the diagnostic variable.
  SpeciesColumnView& aerosol_var(const std::string& name);

  /// Returns a const view storing the modal aerosol diagnostic variable with
  /// the given name and mode index. If the variable does not yet exist, this
  /// throws an exception.
  /// @param [in] name The name of the diagnostic variable.
  const SpeciesColumnView& aerosol_var(const std::string& name) const;

  /// Returns true if a variable defined for each gas species exists within
  /// this object with the given name, false otherwise.
  /// @param [in] name The name of the diagnostic variable of interest.
  bool has_gas_var(const std::string& name) const;

  /// Creates a diagnostic gas species variable with the given name within this
  /// object.
  /// @param [in] name The name of the modal diagnostic variable to be created.
  void create_gas_var(const std::string& name);

  /// Returns the view storing a diagnostic variable, defined for each gas
  /// species, with the given name. If the variable does not yet exist, this
  /// throws an exception.
  /// @param [in] name The name of the diagnostic variable.
  SpeciesColumnView& gas_var(const std::string& name);

  /// Returns the const view storing a diagnostic variable, defined for each gas
  /// species, with the given name. If the variable does not yet exist, this
  /// throws an exception.
  /// @param [in] name The name of the diagnostic variable.
  const SpeciesColumnView& gas_var(const std::string& name) const;

  /// Returns true if the given modal variable exists within this object,
  /// false otherwise.
  /// @param [in] name The name of the diagnostic variable of interest.
  bool has_modal_var(const std::string& name) const;

  /// Creates a diagnostic modal variable with the given name within this object.
  /// @param [in] name The name of the modal diagnostic variable to be created.
  void create_modal_var(const std::string& name);

  /// Returns the view storing the mode-specific diagnostic variable with the
  /// given name. If the variable does not yet exist, this throws an exception.
  /// @param [in] name The name of the diagnostic variable.
  ModalColumnView& modal_var(const std::string& name);

  /// Returns a const view storing the mode-specific diagnostic variable with
  /// the given name. If the variable does not yet exist, this throws an
  /// exception.
  /// @param [in] name The name of the diagnostic variable.
  const ModalColumnView& modal_var(const std::string& name) const;

  private:

  // Number of aerosol species in each mode.
  std::vector<int> num_aero_species_;

  // Number of distinct aerosol populations.
  int num_aero_populations_;

  // Number of gas species.
  int num_gases_;

  /// Number of vertical levels per column in the system.
  const int num_levels_;

  // Named diagnostic variables.
  std::map<std::string, ColumnView> vars_;
  std::map<std::string, SpeciesColumnView> aero_vars_;
  std::map<std::string, SpeciesColumnView> gas_vars_;
  std::map<std::string, ModalColumnView> modal_vars_;
};

}

#endif
