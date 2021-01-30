#include "haero/prognostics.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("prognostics_ctor", "") {

  SECTION("single_mode_single_species") {
    // Create a set of prognostics and add some modes and species to it.
    int num_levels = 72;
    Kokkos::View<PackType**> int_aerosols("interstitial aerosols", 1, num_levels);
    Kokkos::View<PackType**> cld_aerosols("cloudborne aerosols", 1, num_levels);
    int num_gases = 1;
    Kokkos::View<PackType**> gases("gases", num_gases, num_levels);
    int num_modes = 1;
    Kokkos::View<PackType**> modal_concs("modal number concs", num_modes,
                                         num_levels);
    Prognostics progs(num_modes, {1}, num_gases, num_levels,
                      int_aerosols, cld_aerosols, gases, modal_concs);
    REQUIRE(progs.num_aerosol_modes() == num_modes);
    REQUIRE(progs.num_gases() == num_gases);
    REQUIRE(progs.num_levels() == num_levels);
    REQUIRE(progs.num_aerosol_populations() == 1);
  }
}

