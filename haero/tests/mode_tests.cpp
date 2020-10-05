#include "haero/mode.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("mode_ctor", "") {
  // Create a representation of the Aitken mode.
  Mode aitken("Aitken", 0.0087, 0.052, 1.6);
  REQUIRE(aitken.name == "Aitken");
  REQUIRE(aitken.min_diameter == 0.0087);
  REQUIRE(aitken.max_diameter == 0.052);
  REQUIRE(aitken.mean_std_dev == 1.6);
}

