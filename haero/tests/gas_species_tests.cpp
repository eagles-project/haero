#include "haero/gas_species.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("gas_species_ctor", "") {
  // Create a representation of sulfate.
  GasSpecies co2("carbon dioxide", "CO2", 44.009);
  REQUIRE(co2.name() == "carbon dioxide");
  REQUIRE(co2.symbol() == "CO2");
  REQUIRE(co2.molecular_weight == 44.009);
}



