#include "haero/aero_tendencies.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("aero_tendencies_ctor", "") {

  // Create a state and add some modes and species to it.
  AeroState state(10, 72);
  // TODO

  // Try to create tendencies for this state before it's assembled.
  REQUIRE_THROWS(AeroTendencies(state));

  // Assemble the state, create tendencies for it, and make sure the vitals
  // match up.
  state.assemble();
  AeroTendencies tends(state);
  REQUIRE(tends.num_aerosol_modes() == state.num_aerosol_modes());
  for (int m = 0; m < tends.num_aerosol_modes(); ++m) {
    REQUIRE(tends.num_aerosol_species(m) == state.num_aerosol_species(m));
  }
  REQUIRE(tends.num_gas_species() == state.num_gas_species());
  REQUIRE(tends.num_columns() == state.num_columns());
  REQUIRE(tends.num_levels() == state.num_levels());
  for (int m = 0; m < tends.num_aerosol_modes(); ++m) {
    const auto& tends_int_aeros = tends.interstitial_aerosols(m);
    const auto& state_int_aeros = state.interstitial_aerosols(m);
    REQUIRE(tends_int_aeros.extent(1) == state_int_aeros.extent(1));
    REQUIRE(tends_int_aeros.extent(2) == state_int_aeros.extent(2));
    REQUIRE(tends_int_aeros.extent(3) == state_int_aeros.extent(3));

    const auto& tends_cld_aeros = tends.cloudborne_aerosols(m);
    const auto& state_cld_aeros = state.cloudborne_aerosols(m);
    REQUIRE(tends_cld_aeros.extent(1) == state_cld_aeros.extent(1));
    REQUIRE(tends_cld_aeros.extent(2) == state_cld_aeros.extent(2));
    REQUIRE(tends_cld_aeros.extent(3) == state_cld_aeros.extent(3));

    const auto& tends_modal_num_density = tends.modal_num_density(m);
    const auto& state_modal_num_density = state.modal_num_density(m);
    REQUIRE(tends_modal_num_density.extent(1) == state_modal_num_density.extent(1));
    REQUIRE(tends_modal_num_density.extent(2) == state_modal_num_density.extent(2));
  }
  REQUIRE(tends.num_levels() == state.num_levels());
}

