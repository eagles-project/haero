#include "haero/species.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("aero_species_ctor", "") {
  // Create a representation of sulfate.
  // FIXME: Maybe use better numbers.
  Species sulfate("Sulfate", "SO4", 96.0, 1.0, 1.0);
  REQUIRE(sulfate.name() == "Sulfate");
  REQUIRE(sulfate.symbol() == "SO4");
  REQUIRE(sulfate.molecular_weight == 96.0);
  REQUIRE(sulfate.crystalization_point == 1.0);
  REQUIRE(sulfate.deliquescence_point == 1.0);
}

TEST_CASE("gas_species_ctor", "") {
  // Create a representation of sulfuric acid.
  // FIXME: Maybe use better numbers.
  Species sulfuric_acid("Sulfuric acid", "H2SO4", 1.0, 1.0, 1.0);
  REQUIRE(sulfuric_acid.name() == "Sulfuric acid");
  REQUIRE(sulfuric_acid.symbol() == "H2SO4");
  REQUIRE(sulfuric_acid.molecular_weight == 1.0);
  REQUIRE(sulfuric_acid.crystalization_point == 1.0);
  REQUIRE(sulfuric_acid.deliquescence_point == 1.0);
}

