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
/// This type stores state information for an aerosol system. It stores
/// * mixing ratios of for each mode-specific interstitial and cloud-borne
///   aerosol species
/// * interstitial and cloud-borne number concentrations for each aerosol mode
/// * mole fractions for gas species
class Diagnostics final {
  public:

  /// This is the device on which the Diagnostics stores its data.
  using DeviceType = ekat::KokkosTypes<ekat::DefaultDevice>;

  /// This type represents vectorizable packs of Reals of length HAERO_PACK_SIZE.
  using PackType = ekat::Pack<Real, HAERO_PACK_SIZE>;

  /// This type represents a multidimensional array mapping a column and
  /// vertical level to a pack.
  /// * The column is identified by the index i.
  /// * The vertical level identified by the index k.
  /// So view[i][k] yields the desired pack.
  /// Our views are unmanaged in general, to allow a host model to assume
  /// responsibility for managing resources.
  using ColumnView = Kokkos::View<PackType**>;

  /// This type represents a multidimensional array mapping a column, species,
  /// and vertical level to a pack.
  /// * The column is identified by the index i.
  /// * The species is identified by the index s.
  /// * The vertical level identified by the index k.
  /// So view[i][s][k] yields the desired pack.
  /// Our views are unmanaged in general, to allow a host model to assume
  /// responsibility for managing resources.
  using ColumnSpeciesView = Kokkos::View<PackType***>;

  /// This type represents a multidimensional array mapping a mode, column, and
  /// vertical level to a pack.
  /// * The mode is identified by the index m.
  /// * The column is identified by the index i.
  /// * The vertical level identified by the index k.
  /// So view[m][i][k] yields the desired pack.
  /// Our views are unmanaged in general, to allow a host model to assume
  /// responsibility for managing resources.
  using ModalColumnView = Kokkos::View<PackType***>;

  /// Creates an empty Diagnostics to which data can be added.
  /// @param [in] num_columns The number of vertical columns stored by the state
  /// @param [in] num_levels The number of vertical levels per column stored by
  ///                        the state
  /// @param [in] num_aero_species A vector whose ith entry is the number of
  ///                              aerosol species in the ith mode.
  /// @param [in] num_gas_species The number of gas species in the atmosphere.
  Diagnostics(int num_columns, int num_levels,
              const std::vector<int>& num_aero_species,
              int num_gas_species);

  /// Destructor.
  ~Diagnostics();

  // --------------------------------------------------------------------------
  //                                Metadata
  // --------------------------------------------------------------------------

  /// Returns the number of aerosol modes in the state.
  int num_aerosol_modes() const;

  /// Returns the number of aerosol species in the mode with the given index.
  /// @param [in] mode_index The index of the desired mode.
  int num_aerosol_species(int mode_index) const;

  /// Returns the number of gas species in the state.
  int num_gas_species() const;

  /// Returns the number of independent atmospheric columns in the state.
  int num_columns() const;

  /// Returns the number of vertical levels per column in the state.
  int num_levels() const;

  // --------------------------------------------------------------------------
  //                                  Data
  // --------------------------------------------------------------------------

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
  /// @param [in] mode_index The index of the desired mode.
  const ModalColumnView& modal_var(const std::string& name) const;

  private:

  /// Number of columns in the state.
  int num_columns_;

  /// Number of vertical levels per column in the state.
  int num_levels_;

  /// Aerosol species names within each mode.
  std::vector<int> num_aero_species_;

  /// Gas species names.
  int num_gas_species_;

  // Named diagnostic variables.
  std::map<std::string, ColumnView> vars_;
  std::map<std::string, ModalColumnView> modal_vars_;
};

}

#endif
