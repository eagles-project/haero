#include "haero/aerosol_species.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("aero_species_ctor", "") {
  // Create a representation of sulfate.
  AerosolSpecies sulfate("Sulfate", "SO4", 115.107, 6.95e-8, 1770, 0.507);
  REQUIRE(sulfate.name() == "Sulfate");
  REQUIRE(sulfate.symbol() == "SO4");
  REQUIRE(sulfate.molecular_weight == 115.107);
  REQUIRE(sulfate.dry_radius == 6.95e-8);
  REQUIRE(sulfate.density == 1770);
  REQUIRE(sulfate.hygroscopicity == 0.507);
}



