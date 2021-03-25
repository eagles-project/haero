#include "haero/mode.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("mode_ctor", "") {
  // Create a representation of the Aitken mode.
  Mode aitken("Aitken", 8.7e-9, 5.2e-8, 1.6, 0.35, 0.8);
  REQUIRE(aitken.name() == "Aitken");
  REQUIRE(aitken.min_diameter == 8.7e-9);
  REQUIRE(aitken.max_diameter == 5.2e-8);
  REQUIRE(aitken.mean_std_dev == 1.6);
  REQUIRE(aitken.arithmetic_mean_diam() == 0.5*(8.7e-9 + 5.2e-8));
}

