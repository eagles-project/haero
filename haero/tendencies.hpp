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

  Tendencies() = default;
  /// Destructor.
  KOKKOS_FUNCTION
  ~Tendencies();

  /// Interstitial aerosol species mixing ratio tendencies.
  SpeciesColumnView interstitial_aerosols;

  /// Cloud-borne aerosol species mixing ratio tendencies.
  SpeciesColumnView cloud_aerosols;

  /// Gas mole fraction tendencies.
  SpeciesColumnView gases;

  /// Interstitial number mixing ratio tendencies per mode.
  ModeColumnView interstitial_num_mix_ratios;

  /// Cloudborne number mixing ratio tendencies per mode.
  ModeColumnView cloud_num_mix_ratios;

  // --------------------------------------------------------------------------
  //                         Mathematical Operations
  // --------------------------------------------------------------------------

  /// Scales all tendencies by the given constant factor, in place.
  /// @param [in] factor The factor by which the tendencies are scaled.
  /// @returns a reference to the tendencies, which have been scaled.
  const Tendencies& scale(Real factor) const;

  inline const Tendencies& set_zero() { return scale(0); }

  /// Acculumates the given set of tendencies into this set, summing the values
  /// of the prognostic variables in place.
  /// @param [in] tendencies The tendencies to be summed into this object.
  void accumulate(const Tendencies& tendencies);

 private:
};

}  // namespace haero

#endif
