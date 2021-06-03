#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "haero/mode.hpp"
#include "haero/physical_constants.hpp"

using namespace haero;

/// @breif Test the following methods by calculating expected values by hand:
/// - mean_particle_volume_from_diameter
/// - min_vol_to_num_ratio
/// - max_vol_to_num_ratio
TEST_CASE("mode_ctor", "") {
  // Create a representation of the Aitken mode.
  Mode aitken("Aitken", /*min_diameter=*/8.7e-9, /*max_diameter=*/5.2e-8,
              /*mean_std_dev=*/1.6, /*crystallization_pt=*/0.35,
              /*deliquescence_pt=*/0.8);
  REQUIRE(aitken.name() == "Aitken");
  REQUIRE(aitken.min_diameter == Real(8.7e-9));
  REQUIRE(aitken.max_diameter == Real(5.2e-8));
  REQUIRE(aitken.mean_std_dev == Real(1.6));

  // Verify `mean_particle_volume_from_diameter` calculation
  {
    static constexpr Real geom_diam = 5.2e-8;
    static constexpr Real mean_std_dev = 1.6;
    REQUIRE(FloatingPoint<Real>::rel(
        aitken.mean_particle_volume_from_diameter<Real>(geom_diam),
        (constants::pi_sixth * cube(geom_diam) *
         exp(4.5 * square(log(mean_std_dev)))),
        5 * std::numeric_limits<Real>::epsilon()));
  }

  // Verify `min_vol_to_num_ratio` calculation
  {
    // compute min_vol_to_num ratio
    const Real comp_min_vol_to_num_ratio =
        1 / aitken.mean_particle_volume_from_diameter<Real>(5.2e-8);

    // compute relative difference for min_vol_to_num_ratio
    REQUIRE(FloatingPoint<Real>::rel(comp_min_vol_to_num_ratio,
                                     aitken.min_vol_to_num_ratio<Real>(),
                                     5 * std::numeric_limits<Real>::epsilon()));
  }

  // Verify `max_vol_to_num_ratio` calculation
  {
    // compute max_vol_to_num ratio
    const Real comp_max_vol_to_num_ratio =
        1 / aitken.mean_particle_volume_from_diameter<Real>(8.7e-9);

    // compute relative difference for max_vol_to_num_ratio
    REQUIRE(FloatingPoint<Real>::rel(comp_max_vol_to_num_ratio,
                                     aitken.max_vol_to_num_ratio<Real>(),
                                     5 * std::numeric_limits<Real>::epsilon()));
  }
}
