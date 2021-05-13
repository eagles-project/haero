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

  /// Modal interstitial aerosol species mixing ratio tendencies.
  SpeciesColumnView interstitial_aerosols;

  /// Modal cloud-borne aerosol species mixing ratio tendencies.
  SpeciesColumnView cloud_aerosols;

  /// Gas mole fraction tendencies.
  SpeciesColumnView gases;

  /// Modal number concentration tendencies.
  ModalColumnView modal_num_concs;

  // --------------------------------------------------------------------------
  //                         Mathematical Operations
  // --------------------------------------------------------------------------

  /// Scales all tendencies by the given constant factor, in place.
  /// @param [in] factor The factor by which the tendencies are scaled.
  /// @returns a reference to the tendencies, which have been scaled.
  Tendencies& scale(Real factor);

  inline Tendencies& set_zero() { return scale(0); }

  /// Acculumates the given set of tendencies into this set, summing the values
  /// of the prognostic variables in place.
  /// @param [in] tendencies The tendencies to be summed into this object.
  void accumulate(const Tendencies& tendencies);

 private:
};

}  // namespace haero

#endif
