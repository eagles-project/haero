#ifndef HAERO_STATE_HPP
#define HAERO_STATE_HPP

#include "haero/mode.hpp"
#include "haero/species.hpp"
#include "haero/view_pack_helpers.hpp"

namespace haero {

/// This type stores state information for an aerosol system. It stores
/// * mixing ratios of for each mode-specific interstitial and cloud-borne
///   aerosol species
/// * interstitial and cloud-borne number concentrations for each aerosol mode
/// * mole fractions for gas species
class State final {
  public:

  /// This type represents a multidimensional array mapping a set of columns
  /// to a pack, with a column identified by the index i and a given vertical
  /// level identified by the index k: view[i][k] -> pack
  using MultiColumnView = view_2d_pack_type;

  /// This type represents an array storing pack data for a single column,
  /// with a given vertical level identified by the index k: view[k] -> pack
  using SingleColumnView = view_1d_pack_type;

  /// Creates a State using information supplied. Usually, you don't want to
  /// call this constructor explicitly--create a State using a Context object
  /// instead.
  /// @param [in] num_modes the number of aerosol modes stored by the state
  /// @param [in] num_modal_aero_species an array mapping each mode to the
  ///                                    number of aerosol species it contains
  ///                                    (both interstitial and cloud-borne)
  /// @param [in] num_gas_species the number of gas species stored by the state
  /// @param [in] num_columns the number of vertical columns stored by the state
  /// @param [in] num_levels the number of vertical levels per column stored by
  ///                        the state
  State(int num_modes,
        const std::vector<int>& num_modal_species,
        int num_gas_species,
        int num_columns,
        int num_levels);

  /// Destructor.
  ~State();

  private:

  /// Number of aerosol modes in the state.
  int num_modes_;

  /// Index offsets of aerosol species w.r.t. to mode index.
  std::vector<int> modal_species_offsets_;

  /// Number of gas species in the state.
  int num_gas_species_;

  /// Number of columns in the state.
  int num_columns_;

  /// Number of vertical levels per column in the state.
  int num_levels_;

  /// Modal interstitial aerosol species mixing ratios.
  /// interstitial_aerosols_[m][s][i][k] -> mixing ratio of aerosol species s
  /// within mode m, located at vertical level k within column i.
  std::vector<std::vector<MultiColumnView> > interstitial_aerosols_;

  /// Modal cloud-borne aerosol species mixing ratios.
  /// cloud_borne_aerosols_[m][s][i][k] -> mixing ratio of aerosol species s
  /// within mode m, located at vertical level k within column i.
  std::vector<std::vector<MultiColumnView> > cloud_borne_aerosols_;

  /// Gas mole fractions.
  /// gases_[s][i][k] -> mole fraction of gas species s located at vertical
  /// level k within column i.
  std::vector<MultiColumnView> gases_;

  /// Modal number densities.
  /// modal_n_[m][i][k] -> number density of mode m located at vertical level k
  /// within column i.
  std::vector<MultiColumnView> modal_n_;
};

}

#endif
