#include "haero/species.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("aero_species_ctor", "") {
  // Create a representation of sulfate.
  Species sulfate("Sulfate", "SO4");
  REQUIRE(sulfate.name == "Sulfate");
  REQUIRE(sulfate.symbol == "SO4");
}

TEST_CASE("gas_species_ctor", "") {
  // Create a representation of sulfuric acid.
  Species sulfuric_acid("Sulfuric acid", "H2SO4");
  REQUIRE(sulfuric_acid.name == "Sulfuric acid");
  REQUIRE(sulfuric_acid.symbol == "H2SO4");
}

