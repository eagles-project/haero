#include "haero/prognostics.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("prognostics_ctor", "") {

  SECTION("empty_and_unassembled") {
    Prognostics progs(10, 72);
    REQUIRE(not progs.is_assembled());
    REQUIRE(progs.num_aerosol_modes() == 0);
    REQUIRE(progs.num_gas_species() == 0);
    REQUIRE(progs.num_columns() == 10);
    REQUIRE(progs.num_levels() == 72);

    // We can't access views before the state is assembled (especially for
    // non-existent modes, etc!)
    REQUIRE_THROWS(progs.num_aerosol_species(0) == 0);
    REQUIRE_THROWS(progs.interstitial_aerosols(0));
    REQUIRE_THROWS(progs.cloudborne_aerosols(0));
    REQUIRE_THROWS(progs.gas_mole_fractions());
    REQUIRE_THROWS(progs.modal_num_densities());
  }

  SECTION("with_modes_and_species") {
  }
}

