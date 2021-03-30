#include "haero/aerosol_species.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

static constexpr Real g_to_kg = 0.001;

TEST_CASE("aero_species_ctor", "") {
  // Create a representation of sulfate.
  AerosolSpecies sulfate("Sulfate", "SO4", g_to_kg*115.107, 6.95e-8, 1770, 0.507);
  REQUIRE(sulfate.name() == "Sulfate");
  REQUIRE(sulfate.symbol() == "SO4");
  REQUIRE(sulfate.molecular_weight == g_to_kg*115.107);
  REQUIRE(sulfate.dry_radius == 6.95e-8);
  REQUIRE(sulfate.density == 1770);
  REQUIRE(sulfate.hygroscopicity == 0.507);
}



