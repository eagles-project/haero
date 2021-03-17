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

  /// Creates a fully-functional set of tendencies that work with the given
  /// aerosol state.
  /// @param [in] prognostics An assembled Prognostics object that provides the
  ///                         necessary information to create a set of
  ///                         corresponding tendencies.
  explicit Tendencies(const Prognostics& prognostics);

  /// Destructor.
  KOKKOS_FUNCTION
  ~Tendencies();

  /// Returns the view storing interstitial aerosol mass mixing ratio
  /// tendencies [kg aerosol / kg dry air / s].
  KOKKOS_INLINE_FUNCTION
  SpeciesColumnView interstitial_aerosols() {
    return int_aero_species_;
  }

  /// Returns the view storing interstitial aerosol mass mixing ratio
  /// tendencies [kg aerosol / kg dry air / s] (const).
  KOKKOS_FUNCTION
  const SpeciesColumnView interstitial_aerosols() const;

  /// Returns the view storing cloud-borne aerosol mass mixing ratio
  /// tendencies [kg aerosol / kg dry air / s].
  KOKKOS_FUNCTION
  SpeciesColumnView cloudborne_aerosols();

  /// Returns the view storing cloud-borne aerosol mass mixing ratio
  /// tendencies [kg aerosol / kg dry air / s] (const).
  KOKKOS_FUNCTION
  const SpeciesColumnView cloudborne_aerosols() const;

  /// Returns the view storing the mass mixing ratio tendencies for gas species
  /// [kg gas / kg dry air / s].
  KOKKOS_FUNCTION
  SpeciesColumnView gases();

  /// Returns the view storing the mass mixing ratio tendencies for gas species
  /// [kg gas / kg dry air / s] (const).
  KOKKOS_FUNCTION
  const SpeciesColumnView gases() const;

  /// Returns the view storing the modal number concentration tendencies
  /// [# / kg dry air / s].
  KOKKOS_FUNCTION
  ModalColumnView modal_num_concs();

  /// Returns the view storing the modal number concentration tendencies
  /// [# / kg dry air / s] (const).
  KOKKOS_FUNCTION
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
