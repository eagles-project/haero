#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "haero/constants.hpp"
#include "haero/mode.hpp"

using namespace haero;

/// @brief Test the following methods by calculating expected values by hand:
/// - mean_particle_volume_from_diameter
/// - min_vol_to_num_ratio
/// - max_vol_to_num_ratio
TEST_CASE("mode_ctor", "") {
  // Set constants
  static constexpr Real min_diameter = 8.7e-9;
  static constexpr Real nom_diameter = 2.6e-08;
  static constexpr Real max_diameter = 5.2e-8;
  static constexpr Real mean_std_dev = 1.6;
  static constexpr Real crystallization_pt = 0.35;
  static constexpr Real deliquescence_pt = 0.8;

  // Create a representation of the Aitken mode.
  Mode aitken("Aitken", min_diameter, nom_diameter, max_diameter, mean_std_dev,
              crystallization_pt, deliquescence_pt);
  REQUIRE(aitken.name() == "Aitken");
  REQUIRE(aitken.min_diameter == Real(min_diameter));
  REQUIRE(aitken.nom_diameter == Real(nom_diameter));
  REQUIRE(aitken.max_diameter == Real(max_diameter));
  REQUIRE(aitken.mean_std_dev == Real(mean_std_dev));

  // Verify `mean_particle_volume_from_diameter` calculation
  REQUIRE(FloatingPoint<Real>::rel(
      aitken.mean_particle_volume_from_diameter<Real>(max_diameter),
      (Constants::pi_sixth * cube(max_diameter) *
       exp(4.5 * square(log(mean_std_dev)))),
      5 * std::numeric_limits<Real>::epsilon()));

  // Verify `min_vol_to_num_ratio` calculation
  {
    // compute min_vol_to_num ratio
    const Real comp_min_vol_to_num_ratio =
        1 / aitken.mean_particle_volume_from_diameter<Real>(max_diameter);

    // compute relative difference for min_vol_to_num_ratio
    REQUIRE(FloatingPoint<Real>::rel(comp_min_vol_to_num_ratio,
                                     aitken.min_vol_to_num_ratio<Real>(),
                                     5 * std::numeric_limits<Real>::epsilon()));
  }

  // Verify `max_vol_to_num_ratio` calculation
  {
    // compute max_vol_to_num ratio
    const Real comp_max_vol_to_num_ratio =
        1 / aitken.mean_particle_volume_from_diameter<Real>(min_diameter);

    // compute relative difference for max_vol_to_num_ratio
    REQUIRE(FloatingPoint<Real>::rel(comp_max_vol_to_num_ratio,
                                     aitken.max_vol_to_num_ratio<Real>(),
                                     5 * std::numeric_limits<Real>::epsilon()));
  }

  // Verify `nom_vol_to_num_ratio` calculation
  {
    // compute nom_vol_to_num ratio
    const Real comp_nom_vol_to_num_ratio =
        1 / Constants::pi_sixth * cube(aitken.nom_diameter) *
        exp(4.5 * (log(aitken.mean_std_dev)));

    // compute relative difference for nom_vol_to_num_ratio
    REQUIRE(FloatingPoint<Real>::rel(comp_nom_vol_to_num_ratio,
                                     aitken.nom_vol_to_num_ratio<Real>(),
                                     5 * std::numeric_limits<Real>::epsilon()));
  }
}
