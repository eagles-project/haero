#ifndef HAERO_TENDENCIES_HPP
#define HAERO_TENDENCIES_HPP

#include "haero/prognostics.hpp"

namespace haero {

/// @class Tendencies
/// This type stores time derivatives ("tendencies") for the prognostic
/// variables in an aerosol system. These variables are:
/// * mixing ratios of for each mode-specific interstitial and cloud-borne
///   aerosol species
/// * mass mixing ratios for gas species
/// * number concentrations for each aerosol mode
/// All views used by a Tendencies object are managed by that object itself.
class Tendencies final {
  public:

  using kokkos_device_type = ekat::KokkosTypes<ekat::DefaultDevice>;
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

  /// Creates a fully-functional set of tendencies that work with the given
  /// aerosol state. All tendencie views are managed
  /// @param [in] prognostics An assembled Prognostics object that provides the
  ///                         necessary information to create a set of
  ///                         corresponding tendencies.
  explicit Tendencies(const Prognostics& prognostics);

  /// Destructor.
  KOKKOS_FUNCTION
  ~Tendencies();

  /// Returns the view storing interstitial aerosol mass mixing ratio
  /// tendencies [kg aerosol / kg dry air / s].
  SpeciesColumnView interstitial_aerosols();

  /// Returns the view storing interstitial aerosol mass mixing ratio
  /// tendencies [kg aerosol / kg dry air / s] (const).
  const SpeciesColumnView interstitial_aerosols() const;

  /// Returns the view storing cloud-borne aerosol mass mixing ratio
  /// tendencies [kg aerosol / kg dry air / s].
  SpeciesColumnView cloudborne_aerosols();

  /// Returns the view storing cloud-borne aerosol mass mixing ratio
  /// tendencies [kg aerosol / kg dry air / s] (const).
  const SpeciesColumnView cloudborne_aerosols() const;

  /// Returns the view storing the mass mixing ratio tendencies for gas species
  /// [kg gas / kg dry air / s].
  SpeciesColumnView gases();

  /// Returns the view storing the mass mixing ratio tendencies for gas species
  /// [kg gas / kg dry air / s] (const).
  const SpeciesColumnView gases() const;

  /// Returns the view storing the modal number concentration tendencies
  /// [# / kg dry air / s].
  ModalColumnView modal_num_concs();

  /// Returns the view storing the modal number concentration tendencies
  /// [# / kg dry air / s] (const).
  const ModalColumnView modal_num_concs() const;

  // --------------------------------------------------------------------------
  //                         Mathematical Operations
  // --------------------------------------------------------------------------

  /// Scales all tendencies by the given constant factor, in place.
  /// @param [in] factor The factor by which the tendencies are scaled.
  /// @returns a reference to the tendencies, which have been scaled.
  Tendencies& scale(Real factor);

  /// Acculumates the given set of tendencies into this set, summing the values
  /// of the prognostic variables in place.
  /// @param [in] tendencies The tendencies to be summed into this object.
  void accumulate(const Tendencies& tendencies);

  private:

  /// Modal interstitial aerosol species mixing ratio tendencies.
  SpeciesColumnView int_aero_species_;

  /// Modal cloud-borne aerosol species mixing ratio tendencies.
  SpeciesColumnView cld_aero_species_;

  /// Gas mole fraction tendencies.
  SpeciesColumnView gases_;

  /// Modal number concentration tendencies.
  ModalColumnView modal_num_concs_;
};

}

#endif
