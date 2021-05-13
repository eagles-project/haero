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
    int num_levels = 72;
    int num_vert_packs = num_levels / HAERO_PACK_SIZE;
    if (num_vert_packs * HAERO_PACK_SIZE < num_levels) {
      num_vert_packs++;
    }
    kokkos_device_type::view_2d<PackType> int_aerosols("interstitial aerosols",
                                                       1, num_vert_packs);
    kokkos_device_type::view_2d<PackType> cld_aerosols("cloudborne aerosols", 1,
                                                       num_vert_packs);
    int num_gases = 1;
    kokkos_device_type::view_2d<PackType> gases("gases", num_gases,
                                                num_vert_packs);
    int num_modes = 1;
    kokkos_device_type::view_2d<PackType> int_num_concs(
        "interstitial number concs", num_modes, num_vert_packs);
    kokkos_device_type::view_2d<PackType> cld_num_concs(
        "cloudborne number concs", num_modes, num_vert_packs);

    Prognostics progs(num_modes, {1}, num_gases, num_levels, int_aerosols,
                      cld_aerosols, gases, int_num_concs, cld_num_concs);
    REQUIRE(progs.num_aerosol_modes() == num_modes);
    REQUIRE(progs.num_gases() == num_gases);
    REQUIRE(progs.num_levels() == num_levels);
    REQUIRE(progs.num_aerosol_populations() == 1);
  }
}
