#ifndef HAERO_PROGNOSTICS_HPP
#define HAERO_PROGNOSTICS_HPP

#include <map>
#include <vector>

#include "haero/haero.hpp"
#include "haero/mode.hpp"
#include "haero/aerosol_species.hpp"
#include "haero/view_pack_helpers.hpp"

namespace haero {

/// This is a forward declaration, since we refer to Tendencies below.
class Tendencies;

/// @class Prognostics
/// This type stores the prognostic variables for an atmospheric column in an
/// aerosol system. Specifically,It stores
/// * mass mixing ratios of for each mode-specific interstitial and cloud-borne
///   aerosol species
/// * interstitial and cloud-borne number concentrations for each aerosol mode
/// * mass mixing ratios for gas species
class Prognostics final {
  public:

  /// Creates a Prognostics object that can store aerosol data can be added.
  /// This constructor accepts a number of Kokkos View objects, managed by the
  /// host model, that store aerosol data. The Prognostics object
  /// creates its own views for these views.
  /// @param [in] num_aerosol_modes The number of aerosol modes in the system
  /// @param [in] num_aerosol_species A vector of length num_aerosol_modes whose
  ///                                 ith entry is the number of aerosol species
  ///                                 in the ith mode
  /// @param [in] num_gases The number of gas species present in the column
  /// @param [in] num_levels the number of vertical levels per column stored by
  ///                        the state
  /// @param [in] int_aerosols A rank-2 Kokkos View:
  ///                          * the first index uniquely identifies a (mode,
  ///                            species) combination
  ///                          * the second index uniquely identifies a Pack
  ///                            whose data are interstitial aerosol mass mixing
  ///                            ratios [kg aerosol/kg dry air] for a number of
  ///                            adjacent vertical levels equal to
  ///                            HAERO_PACK_SIZE (possibly padded)
  /// @param [in] cld_aerosols A rank-2 Kokkos View similar to int_aerosols, but
  ///                          storing data for cloud-borne aerosol mass mixing
  ///                          ratios
  /// @param [in] gases A rank-2 Kokkos View;
  ///                   * the first index uniquely identifies a gas species
  ///                   * the second index uniquely identifies a Pack whose data
  ///                     are gas mass mixing ratios [kg gas/kg dry air] for a
  ///                     number of adjacent vertical levels equal to
  ///                     HAERO_PACK_SIZE (possibly padded)
  /// @param [in] modal_num_concs A rank-2 Kokkos View;
  ///                             * the first index uniquely identifies an
  ///                               aerosol particle size mode with an
  ///                               associated probability distribution function
  ///                               (PDF)
  ///                             * the second index uniquely identifies a Pack
  ///                               whose data are modal number concentrations
  ///                               [# particles/kg dry air] for a number of
  ///                               adjacent vertical levels equal to
  ///                               HAERO_PACK_SIZE (possibly padded)
  Prognostics(int num_aerosol_modes,
              const std::vector<int>& num_aerosol_species,
              int num_gases,
              int num_levels,
              SpeciesColumnView int_aerosols,
              SpeciesColumnView cld_aerosols,
              SpeciesColumnView gases,
              ModalColumnView   modal_num_concs);

  /// Destructor.
  KOKKOS_FUNCTION
  ~Prognostics();

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
  //                                 Data
  // --------------------------------------------------------------------------

  /// Returns the view storing interstitial aerosol mass mixing ratios
  /// [kg aerosol / kg dry air].
  KOKKOS_INLINE_FUNCTION
  SpeciesColumnView interstitial_aerosols() {
    return int_aero_species_;
  }

  /// Returns the view storing interstitial aerosol mass mixing ratios
  /// [kg aerosol / kg dry air] (const).
  KOKKOS_INLINE_FUNCTION
  const SpeciesColumnView interstitial_aerosols() const {
    return int_aero_species_;
  }

  /// Returns the view storing cloud-borne aerosol mass mixing ratios
  /// [kg aerosol / kg dry air].
  SpeciesColumnView cloudborne_aerosols();

  /// Returns the view storing cloud-borne aerosol mass mixing ratios
  /// [kg aerosol / kg dry air] (const).
  const SpeciesColumnView cloudborne_aerosols() const;

  /// Returns the view storing the mass mixing ratios of gas species
  /// [kg gas / kg dry air].
  SpeciesColumnView gases();

  /// Returns the view storing the mass mixing ratios of gas species
  /// [kg gas / kg dry air] (const).
  const SpeciesColumnView gases() const;

  /// Returns the view storing the modal number concentrations [# / kg dry air].
  ModalColumnView modal_num_concs();

  /// Returns the view storing the modal number concentrations [# / kg dry air]
  /// (const).
  const ModalColumnView modal_num_concs() const;

  // --------------------------------------------------------------------------
  //                         Mathematical Operations
  // --------------------------------------------------------------------------

  /// Scales the given set of tendencies and adds it into this state, summing
  /// the values of the prognostic variables in place.
  /// @param [in] scale_factor The factor by which the tendecies are scaled.
  /// @param [in] tendencies The tendencies to be summed into the state.
  void scale_and_add(Real scale_factor, const Tendencies& tendencies);

  private:

  // Aerosol species names within each mode.
  const view_1d_int_type num_aero_species_;

  // Number of distinct aerosol populations.
  int num_aero_populations_;

  // Number of gas species.
  const int num_gases_;

  // Number of vertical levels.
  const int num_levels_;

  // Modal interstitial aerosol species mixing ratios.
  // int_aerosols_[s][k] -> mixing ratio of aerosol mode/species s within the
  // kth pack of vertical levels.
  SpeciesColumnView int_aero_species_;

  // Modal cloud-borne aerosol species mixing ratios.
  // cld_aerosols_[s][k] -> mixing ratio of aerosol mode/species s within the
  // kth pack of vertical levels.
  SpeciesColumnView cld_aero_species_;

  /// Gas mixing ratios.
  /// gases_[s][k] -> mixing ratio of gas species s within the kth pack of
  /// vertical levels.
  SpeciesColumnView gases_;

  /// Modal number concentrations.
  /// modal_num_concs_[m][k] -> number concentration of mode m within the kth
  /// pack of vertical levels.
  ModalColumnView modal_num_concs_;
};

}

#endif
