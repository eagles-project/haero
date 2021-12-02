#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "haero/prognostics.hpp"
#include "haero/view_pack_helpers.hpp"

using namespace haero;

TEST_CASE("prognostics_ctor", "") {
  SECTION("single_mode_single_species") {
    // Create a set of prognostics and add some modes and species to it.
    using kokkos_device_type = ekat::KokkosTypes<ekat::DefaultDevice>;
    // Create a test configuration and create prognostics for it.
    auto test_config = ModalAerosolConfig::create_simple_test_config();
    int num_levels = 72;
    Prognostics progs(test_config, num_levels);

    REQUIRE(progs.num_aerosol_modes() == test_config.num_modes());
    REQUIRE(progs.num_gases() == test_config.num_gases());
    REQUIRE(progs.num_aerosol_populations() == test_config.num_aerosol_populations);
    REQUIRE(progs.num_levels() == num_levels);
  }
}
