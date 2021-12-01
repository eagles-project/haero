#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "haero/tendencies.hpp"
#include "haero/view_pack_helpers.hpp"

using namespace haero;

TEST_CASE("tendencies_ctor", "") {
  using kokkos_device_type = ekat::KokkosTypes<ekat::DefaultDevice>;
  // Create a set of prognostics and add some modes and species to it.
  int num_levels = 72;
  int num_vert_packs = PackInfo::num_packs(num_levels);
  kokkos_device_type::view_2d<PackType> int_aerosols("interstitial aerosols", 1,
                                                     num_vert_packs);
  kokkos_device_type::view_2d<PackType> cld_aerosols("cloudborne aerosols", 1,
                                                     num_vert_packs);
  int num_gases = 1;
  kokkos_device_type::view_2d<PackType> gases("gases", num_gases,
                                              num_vert_packs);
  int num_modes = 1;
  kokkos_device_type::view_2d<PackType> int_num_mix_ratios(
      "interstitial number mix_ratios", num_modes, num_vert_packs);
  kokkos_device_type::view_2d<PackType> cld_num_mix_ratios(
      "cloud borne number mix_ratios", num_modes, num_vert_packs);

  Prognostics progs(num_levels, int_aerosols,
                    cld_aerosols, gases, int_num_mix_ratios,
                    cld_num_mix_ratios);

  // Now create tendencies for it, and make sure the vitals match up.
  Tendencies tends(progs);
  const auto& tends_int_aeros = tends.interstitial_aerosols;
  const auto& progs_int_aeros = progs.interstitial_aerosols;
  REQUIRE(tends_int_aeros.extent(0) == progs_int_aeros.extent(0));
  REQUIRE(tends_int_aeros.extent(1) == progs_int_aeros.extent(1));

  const auto& tends_cld_aeros = tends.cloud_aerosols;
  const auto& progs_cld_aeros = progs.cloud_aerosols;
  REQUIRE(tends_cld_aeros.extent(0) == progs_cld_aeros.extent(0));
  REQUIRE(tends_cld_aeros.extent(1) == progs_cld_aeros.extent(1));

  const auto& tends_interstitial_num_mix_ratios =
      tends.interstitial_num_mix_ratios;
  const auto& progs_interstitial_num_mix_ratios =
      progs.interstitial_num_mix_ratios;
  REQUIRE(tends_interstitial_num_mix_ratios.extent(0) ==
          progs_interstitial_num_mix_ratios.extent(0));
  REQUIRE(tends_interstitial_num_mix_ratios.extent(1) ==
          progs_interstitial_num_mix_ratios.extent(1));
}
