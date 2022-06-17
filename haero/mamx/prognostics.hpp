#ifndef HAERO_PROGNOSTICS_HPP
#define HAERO_PROGNOSTICS_HPP

#include <map>
#include <vector>

#include "haero/modal_aerosol_config.hpp"

namespace haero {

/// This is a forward declaration, since we refer to Tendencies below.
class Tendencies;

/// @class Prognostics
/// This type stores the prognostic variables for an atmospheric column in an
/// aerosol system. Specifically,It stores
/// * mass mixing ratios of for each mode-specific interstitial and cloud-borne
///   aerosol species
/// * interstitial and cloud-borne number mixing ratios for each aerosol mode
/// * mass mixing ratios for gas species
class Prognostics final {
 public:
  /// Creates a Prognostics object associated with the given configuration for a
  /// column with the given number of vertical levels.
  /// @param [in] config A modal aerosol configuration that defines the aerosol
  ///                    gas species and modes present in the system of interest
  /// @param [in] num_levels the number of vertical levels per column stored by
  ///                        the state
  Prognostics(const ModalAerosolConfig& config, int num_levels);

  /// Creates a Prognostics object that maintains the aerosol data in the given
  /// set of Kokkos Views (which are managed by a host model).
  /// @param [in] config A modal aerosol configuration that defines the aerosol
  ///                    gas species and modes present in the system of interest
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
  /// @param [in] int_mode_num_mix_ratios A rank-2 Kokkos View;
  ///                                     * the first index uniquely identifies
  ///                                       an aerosol particle size mode with
  ///                                       an associated probability
  ///                                       distribution function (PDF)
  ///                                     * the second index uniquely identifies
  ///                                       a Pack whose data are interstitial
  ///                                       number mixing ratios [# particles/
  ///                                       kg dry air] for a number of adjacent
  ///                                       vertical levels equal to
  ///                                       HAERO_PACK_SIZE (possibly padded)
  /// @param [in] cld_mode_num_mix_ratios A rank-2 Kokkos View;
  ///                                     * the first index uniquely identifies
  ///                                       an aerosol particle size mode with
  ///                                       an associated probability
  ///                                       distribution function (PDF)
  ///                                     * the second index uniquely identifies
  ///                                       a Pack whose data are cloudborne
  ///                                       aerosols number mixing ratios [#
  ///                                       particles/kg dry air] for a number
  ///                                       of adjacent vertical levels equal to
  ///                                       HAERO_PACK_SIZE (possibly padded)
  /// @param [in] gases A rank-2 Kokkos View;
  ///                   * the first index uniquely identifies a gas species
  ///                   * the second index uniquely identifies a Pack whose data
  ///                     are gas mass mixing ratios [kg gas/kg dry air] for a
  ///                     number of adjacent vertical levels equal to
  ///                     HAERO_PACK_SIZE (possibly padded)
  Prognostics(const ModalAerosolConfig& config, int num_levels,
              SpeciesColumnView int_aerosols, SpeciesColumnView cld_aerosols,
              ModeColumnView int_mode_num_mix_ratios,
              ModeColumnView cld_mode_num_mix_ratios, SpeciesColumnView gases);

  /// Default constructor and copy constructor is needed to define Views
  /// of Prognostics classes and then copy data into the views.
  Prognostics() = default;
  Prognostics(const Prognostics&) = default;

  /// Destructor.
  KOKKOS_INLINE_FUNCTION
  ~Prognostics() {}

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

  /// Returns the number of vertical levels in the column.
  int num_levels() const;

  // --------------------------------------------------------------------------
  //                                 Data
  // --------------------------------------------------------------------------

  // Modal interstitial aerosol species mixing ratios.
  // int_aerosols_[s][k] -> mixing ratio of aerosol mode/species s within the
  // kth pack of vertical levels.
  SpeciesColumnView interstitial_aerosols;

  // Modal cloud-borne aerosol species mixing ratios.
  // cld_aerosols_[s][k] -> mixing ratio of aerosol mode/species s within the
  // kth pack of vertical levels.
  SpeciesColumnView cloud_aerosols;

  /// Interstitial aerosol modal number mixing ratios [#/kg dry air].
  /// interstitial_num_mix_ratios[m][k] -> interstitial aerosol number
  /// mixing ratio of mode m within the kth pack of vertical levels.
  ModeColumnView interstitial_num_mix_ratios;

  /// Cloudborne aerosol modal number mixing ratios [#/kg dry air]
  /// cloud_num_mix_ratios[m][k] -> cloudborne aerosol number mixing ratio of
  /// mode m within the kth pack of vertical levels.
  ModeColumnView cloud_num_mix_ratios;

  /// Gas mixing ratios.
  /// gases_[s][k] -> mixing ratio of gas species s within the kth pack of
  /// vertical levels.
  SpeciesColumnView gases;

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
  view_1d_int_type num_aero_species_;

  // Number of distinct aerosol populations.
  int num_aero_populations_;

  // Number of gas species.
  int num_gases_;

  // Number of vertical levels.
  int num_levels_;
};

}  // namespace haero

#endif
