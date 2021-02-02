#include "haero/tendencies.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("tendencies_ctor", "") {

  // Create a set of prognostics and add some modes and species to it.
  int num_levels = 72;
  int num_vert_packs = num_levels / HAERO_PACK_SIZE;
  Kokkos::View<PackType**> int_aerosols("interstitial aerosols", 1,
                                        num_vert_packs);
  Kokkos::View<PackType**> cld_aerosols("cloudborne aerosols", 1,
                                        num_vert_packs);
  int num_gases = 1;
  Kokkos::View<PackType**> gases("gases", num_gases, num_vert_packs);
  int num_modes = 1;
  Kokkos::View<PackType**> modal_concs("modal number concs", num_modes,
                                       num_vert_packs);
  Prognostics progs(num_modes, {1}, num_gases, num_levels,
                    int_aerosols, cld_aerosols, gases, modal_concs);

  // Now create tendencies for it, and make sure the vitals match up.
  Tendencies tends(progs);
  const auto& tends_int_aeros = tends.interstitial_aerosols();
  const auto& progs_int_aeros = progs.interstitial_aerosols();
  REQUIRE(tends_int_aeros.extent(0) == progs_int_aeros.extent(0));
  REQUIRE(tends_int_aeros.extent(1) == progs_int_aeros.extent(1));

  const auto& tends_cld_aeros = tends.cloudborne_aerosols();
  const auto& progs_cld_aeros = progs.cloudborne_aerosols();
  REQUIRE(tends_cld_aeros.extent(0) == progs_cld_aeros.extent(0));
  REQUIRE(tends_cld_aeros.extent(1) == progs_cld_aeros.extent(1));

  const auto& tends_modal_num_concs = tends.modal_num_concs();
  const auto& progs_modal_num_concs = progs.modal_num_concs();
  REQUIRE(tends_modal_num_concs.extent(0) == progs_modal_num_concs.extent(0));
  REQUIRE(tends_modal_num_concs.extent(1) == progs_modal_num_concs.extent(1));
}

