#ifndef HAERO_AEROSOL_STATE_HPP
#define HAERO_AEROSOL_STATE_HPP

#include "haero/mode.hpp"
#include "haero/species.hpp"
#include "ekat/ekat_pack.hpp"
#include "haero/mode.hpp"
#include "haero/species.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "kokkos/Kokkos_Core.hpp"
#include <map>
#include <vector>

namespace haero {

/// This type stores state information for an aerosol system. It stores
/// * mixing ratios of for each mode-specific interstitial and cloud-borne
///   aerosol species
/// * interstitial and cloud-borne number concentrations for each aerosol mode
/// * mole fractions for gas species
class AerosolState final {
  public:

  /// This is the device on which the AerosolState stores its data.
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
  using ColumnView = ekat::Unmanaged<Kokkos::View<Real**> >;

  /// This type represents a multidimensional array mapping a column, species,
  /// and vertical level to a pack.
  /// * The column is identified by the index i.
  /// * The species is identified by the index s.
  /// * The vertical level identified by the index k.
  /// So view[i][s][k] yields the desired pack.
  /// Our views are unmanaged in general, to allow a host model to assume
  /// responsibility for managing resources.
  using ColumnSpeciesView = ekat::Unmanaged<Kokkos::View<Real***> >;

  /// Creates an empty AerosolState to which data can be added.
  /// @param [in] num_columns the number of vertical columns stored by the state
  /// @param [in] num_levels the number of vertical levels per column stored by
  ///                        the state
  AerosolState(int num_columns, int num_levels);

  /// Destructor.
  ~AerosolState();

  // --------------------------------------------------------------------------
  //                               State Setup
  // --------------------------------------------------------------------------

  /// Adds an aerosol mode to the state, with a set of species belonging to that
  /// mode. This method creates managed views for the modal aerosol data within
  /// the state. Use this when your host model doesn't need to manage resources
  /// for your aerosol species.
  /// @param [in] mode The name of the mode to which this aerosol species
  ///                       is added.
  /// @param [in] aero_species A list of aerosol species belonging to the mode.
  /// @returns the index of the aerosol mode within the state.
  /// @throws
  int add_aerosol_mode(const Mode& mode,
                       const std::vector<Species>& aero_species);

  /// Adds an aerosol mode to the state, with a set of species belonging to that
  /// mode. This method accepts unmanaged views for the modal aerosol data to
  /// be used by the state, but whose resources are managed elsewhere. Use this
  /// when your host model manages resources for your aerosol species.
  /// @param [in] mode The name of the mode to which this aerosol species
  ///                       is added.
  /// @param [in] aero_species A list of aerosol species belonging to the mode.
  /// @param [in] int_aero_data An unmanaged view to use for storing
  ///                           interstitial aerosol data in this mode.
  /// @param [in] cld_aero_data An unmanaged view to use for storing
  ///                           cloud-borne aerosol data in this mode.
  /// @returns the index of the aerosol mode within the state.
  /// @throws
  int add_aerosol_mode(const Mode& mode,
                       const std::vector<Species>& aero_species,
                       ColumnSpeciesView& int_aero_data,
                       ColumnSpeciesView& cld_aero_data);

  /// Adds a set of gas species to the state. This method creates managed views
  /// for the gas species data within the state. Use this when your host model
  /// doesn't need to manage resources for your gas species.
  /// @param [in] gas_species An list of gas species belonging to the mode.
  ///                         The species are appended to any existing gas
  ///                         species within the state.
  /// @returns the index of the first new gas species added to the state.
  /// @throws
  int add_gas_species(const std::vector<Species>& gas_species);

  /// Adds a set of gas species to the state. This method accepts an unmanaged
  /// view for the gas species data to be used by the state, but whose resources
  /// are managed elsewhere. Use this when your host model manages resources for
  /// your gas species.
  /// @param [in] gas_species An list of gas species to be added to the state.
  ///                         The species are appended to any existing gas
  ///                         species within the state.
  /// @returns the index of the first new gas species added to the state.
  /// @throws
  int add_gas_species(const std::vector<Species>& gas_species,
                      ColumnSpeciesView& int_aero_data);

  /// Finalizes the state, doing any necessary allocations. No data may be
  /// accessed within the state before this is done.
  void finalize();

  /// Returns true if the state has been finalized, false if not.
  bool is_finalized() const;

  // --------------------------------------------------------------------------
  //                           Access to Data and Metadata
  // --------------------------------------------------------------------------

  /// Returns the number of aerosol modes in the state.
  int num_aerosol_modes() const;

  /// Returns the number of aerosol species in the mode with the given index.
  /// @param [in] mode_index The index of the desired mode.
  int num_aerosol_species(int mode_index) const;

  /// Returns the number of gas species in the state.
  int num_gas_species() const;

  /// Returns the view storing interstitial aerosol species mixing fraction data
  /// for the mode with the given index.
  /// @param [in] mode_index The index of the desired mode.
  ColumnSpeciesView& interstitial_aerosols(int mode_index);

  /// Returns the view storing interstitial aerosol species mixing fraction data
  /// for the mode with the given index (const).
  /// @param [in] mode_index The index of the desired mode.
  const ColumnSpeciesView& interstitial_aerosols(int mode_index) const;

  /// Returns the view storing cloud-borne aerosol species mixing fraction data
  /// for the mode with the given index.
  /// @param [in] mode_index The index of the desired mode.
  ColumnSpeciesView& cloudborne_erosols(int mode_index);

  /// Returns the view storing cloud-borne aerosol species data for the mode
  /// with the given index (const).
  /// @param [in] mode_index The index of the desired mode.
  const ColumnSpeciesView& cloudborne_aerosols(int mode_index) const;

  /// Returns the view storing the mole fractions of gas species in the state.
  ColumnSpeciesView& gas_mole_fractions();

  /// Returns the view storing the mole fractions of gas species in the state
  /// (const).
  const ColumnSpeciesView& gas_mole_fractions() const;

  /// Returns the view storing the modal number densities for the state.
  ColumnView& modal_densities();

  /// Returns the view storing the modal number densities for the state (const).
  const ColumnView& modal_densities() const;

  private:

  // These are managed versions of the views above, for views created by the
  // state itself.
  using ManagedColumnView = Kokkos::View<Real**>;
  using ManagedColumnSpeciesView = Kokkos::View<Real***>;

  /// Number of columns in the state.
  int num_columns_;

  /// Number of vertical levels per column in the state.
  int num_levels_;

  /// Number of aerosol modes.
  int num_modes_;

  /// Modal interstitial aerosol species mixing ratios.
  /// interstitial_aerosols_[m][s][i][k] -> mixing ratio of aerosol species s
  /// within mode m, located at vertical level k within column i.
  std::vector<ColumnSpeciesView> int_aero_species_;

  /// Modal cloud-borne aerosol species mixing ratios.
  /// cloud_borne_aerosols_[m][s][i][k] -> mixing ratio of aerosol species s
  /// within mode m, located at vertical level k within column i.
  std::vector<ColumnSpeciesView> cld_aero_species_;

  /// Gas mole fractions.
  /// gases_[s][i][k] -> mole fraction of gas species s located at vertical
  /// level k within column i.
  ColumnSpeciesView gas_mole_fractions_;

  /// Modal number densities.
  /// modal_n_[m][i][k] -> number density of mode m located at vertical level k
  /// within column i.
  ColumnView modal_num_densities_;

  // Lists of managed views to be destroyed with the state.
  std::vector<ManagedColumnView> managed_column_views_;
  std::vector<ManagedColumnSpeciesView> managed_column_species_views_;
};

}

#endif
