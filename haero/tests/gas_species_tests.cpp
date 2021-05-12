#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "haero/gas_species.hpp"

using namespace haero;

static constexpr Real g_to_kg = 0.001;

TEST_CASE("gas_species_ctor", "") {
  // Create a representation of sulfate.
  GasSpecies co2("carbon dioxide", "CO2", "Our favorite greenhouse gas!",
                 g_to_kg * 44.009);
  REQUIRE(co2.name() == "carbon dioxide");
  REQUIRE(co2.symbol() == "CO2");
  REQUIRE(co2.description() == "Our favorite greenhouse gas!");
  REQUIRE(co2.molecular_weight == Real(g_to_kg * 44.009));
}
