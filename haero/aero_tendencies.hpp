#ifndef HAERO_AERO_TENDENCIES_HPP
#define HAERO_AERO_TENDENCIES_HPP

#include "haero/aero_state.hpp"

namespace haero {

/// @class AeroTendencies
/// This type stores time derivatives ("tendencies") for the prognostic
/// variables in an aerosol system. These variables are:
/// * mixing ratios of for each mode-specific interstitial and cloud-borne
///   aerosol species
/// * interstitial and cloud-borne number concentrations for each aerosol mode
/// * mole fractions for gas species
class AeroTendencies final {
  public:

  /// This is the device on which the AeroTendencies stores its data.
  using DeviceType = AeroState::DeviceType;

  /// This type represents vectorizable packs of Reals of length HAERO_PACK_SIZE.
  using PackType = AeroState::PackType;

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

  /// Creates a fully-functional set of tendencies that work with the given
  /// aerosol state.
  /// @param [in] state A finalized AeroState object that provides all the
  ///                   necessary information to create a set of corresponding
  ///                   tendencies.
  explicit AeroTendencies(const AeroState& state);

  /// Destructor.
  ~AeroTendencies();

  /// Returns the number of aerosol mode tendencies.
  int num_aerosol_modes() const;

  /// Returns the number of aerosol species tendencies in the mode with the
  /// given index.
  /// @param [in] mode_index The index of the desired mode.
  int num_aerosol_species(int mode_index) const;

  /// Returns the number of gas species tendencies in the state.
  int num_gas_species() const;

  /// Returns the number of independent atmospheric columns.
  int num_columns() const;

  /// Returns the number of vertical levels per column.
  int num_levels() const;

  /// Returns the view storing interstitial aerosol species mixing fraction
  /// tendencies for the mode with the given index.
  /// @param [in] mode_index The index of the desired mode.
  ColumnSpeciesView& interstitial_aerosols(int mode_index);

  /// Returns the view storing interstitial aerosol species mixing fraction
  /// tendencies for the mode with the given index (const).
  /// @param [in] mode_index The index of the desired mode.
  const ColumnSpeciesView& interstitial_aerosols(int mode_index) const;

  /// Returns the view storing cloud-borne aerosol species mixing fraction
  /// tendencies for the mode with the given index.
  /// @param [in] mode_index The index of the desired mode.
  ColumnSpeciesView& cloudborne_aerosols(int mode_index);

  /// Returns the view storing cloud-borne aerosol species tendencies for the
  /// mode with the given index (const).
  /// @param [in] mode_index The index of the desired mode.
  const ColumnSpeciesView& cloudborne_aerosols(int mode_index) const;

  /// Returns the view storing the mole fraction tendencies for gas species in
  /// the state.
  ColumnSpeciesView& gas_mole_fractions();

  /// Returns the view storing the mole fraction tendencies for gas species in
  /// the state (const).
  const ColumnSpeciesView& gas_mole_fractions() const;

  /// Returns the view storing the number density tendencies for the mode with
  /// the given index.
  /// @param [in] mode_index The index of the desired mode.
  ColumnView& modal_num_density(int mode_index);

  /// Returns the view storing the number density tendencies for the mode with
  /// the given index (const).
  /// @param [in] mode_index The index of the desired mode.
  const ColumnView& modal_num_density(int mode_index) const;

  // --------------------------------------------------------------------------
  //                         Mathematical Operations
  // --------------------------------------------------------------------------

  /// Scales all tendencies by the given constant factor, in place.
  /// @param [in] factor The factor by which the tendencies are scaled.
  /// @returns a reference to the tendencies, which have been scaled.
  AeroTendencies& scale(Real factor);

  /// Acculumates the given set of tendencies into this set, summing the values
  /// of the prognostic variables in place.
  /// @param [in] tendencies The tendencies to be summed into this object.
  void accumulate(const AeroTendencies& tendencies);

  private:

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
  std::vector<ColumnView> modal_num_densities_;
};

}

#endif
