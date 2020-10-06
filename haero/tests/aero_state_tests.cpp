#include "haero/aero_state.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("aero_state_ctor", "") {

  SECTION("empty_and_unassembled") {
    AeroState state(10, 72);
    REQUIRE(not state.is_assembled());
    REQUIRE(state.num_aerosol_modes() == 0);
    REQUIRE(state.num_gas_species() == 0);
    REQUIRE(state.num_columns() == 10);
    REQUIRE(state.num_levels() == 72);

    // We can't access views before the state is assembled (especially for
    // non-existent modes, etc!)
    REQUIRE_THROWS(state.num_aerosol_species(0) == 0);
    REQUIRE_THROWS(state.interstitial_aerosols(0));
    REQUIRE_THROWS(state.cloudborne_aerosols(0));
    REQUIRE_THROWS(state.gas_mole_fractions());
    REQUIRE_THROWS(state.modal_num_density(0));
  }

  SECTION("with_modes_and_species") {
  }
}

