#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/floating_point.hpp"
#include "haero/haero.hpp"
#include "haero/region_of_validity.hpp"

using namespace haero;

TEST_CASE("region_of_validity", "") {
  // Create some prognostics and atmosphere state data.
  int num_levels = 72;
  int num_vert_packs = num_levels / HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels) {
    num_vert_packs++;
  }

  ModalAerosolConfig config = create_simple_test_config();
  int num_modes = config.num_modes();
  int num_gases = config.num_gases();
  int num_pops = config.num_aerosol_populations;

  SpeciesColumnView int_aerosols("interstitial aerosols", num_pops,
                                 num_vert_packs);
  SpeciesColumnView cld_aerosols("cloudborne aerosols", num_pops,
                                 num_vert_packs);
  ModeColumnView int_num_mix_ratios("interstitial number mix ratios", num_modes,
                                    num_vert_packs);
  ModeColumnView cld_num_mix_ratios("cloudborne number mix ratios", num_modes,
                                    num_vert_packs);
  SpeciesColumnView gases("gases", num_gases, num_vert_packs);
  std::vector<int> num_species_per_mode(num_modes);
  for (int m = 0; m < num_modes; ++m) {
    auto species_for_mode = config.aerosol_species_for_mode(m);
    num_species_per_mode[m] = species_for_mode.size();
  }
  Prognostics progs(num_modes, num_species_per_mode, num_gases, num_levels,
                    int_aerosols, cld_aerosols, int_num_mix_ratios,
                    cld_num_mix_ratios, gases);

  ColumnView temp("temperature", num_vert_packs),
      press("pressure", num_vert_packs),
      qv("vapor mixing ratio", num_vert_packs),
      ht("height", num_vert_packs + 1),
      pdel("pressure thickness", num_vert_packs);
  const Real pblh = 1000.0;
  Atmosphere atm(num_levels, temp, press, qv, ht, pdel, pblh);

  // Zero the state data.
  for (int p = 0; p < progs.num_aerosol_populations(); ++p) {
    Kokkos::parallel_for(
        "zero aerosol data", num_vert_packs, KOKKOS_LAMBDA(const int k) {
          int_aerosols(p, k) = 0;
          cld_aerosols(p, k) = 0;
        });
  }
  for (int m = 0; m < num_modes; ++m) {
    Kokkos::parallel_for(
        "zero aerosol mode data", num_vert_packs, KOKKOS_LAMBDA(const int k) {
          int_num_mix_ratios(m, k) = 0;
          cld_num_mix_ratios(m, k) = 0;
        });
  }
  for (int g = 0; g < num_gases; ++g) {
    Kokkos::parallel_for(
        "zero gas data", num_vert_packs,
        KOKKOS_LAMBDA(const int k) { gases(g, k) = 0; });
  }

  Kokkos::parallel_for(
      "zero atmosphere data", num_vert_packs, KOKKOS_LAMBDA(const int k) {
        atm.temperature(k) = 0;
        atm.pressure(k) = 0;
        atm.vapor_mixing_ratio(k) = 0;
      });

  SECTION("ctor") {
    RegionOfValidity rov;

    // Check default bounds.
    REQUIRE(FloatingPoint<Real>::equiv(rov.temp_bounds.first, 0.0));
    REQUIRE(FloatingPoint<Real>::equiv(rov.temp_bounds.second, 500.0));
    REQUIRE(FloatingPoint<Real>::equiv(rov.rel_hum_bounds.first, 0.0));
    REQUIRE(FloatingPoint<Real>::equiv(rov.rel_hum_bounds.second, 1.0));

    // There are defaults for any aerosol population, even those that don't
    // exist!
    for (int p = 0; p < 10; ++p) {
      auto int_aero_bounds = rov.interstitial_aerosol_bounds(p);
      REQUIRE(FloatingPoint<Real>::equiv(int_aero_bounds.first, 0.0));
      REQUIRE(FloatingPoint<Real>::equiv(int_aero_bounds.second, 1.0));

      auto cld_aero_bounds = rov.cloud_aerosol_bounds(p);
      REQUIRE(FloatingPoint<Real>::equiv(cld_aero_bounds.first, 0.0));
      REQUIRE(FloatingPoint<Real>::equiv(cld_aero_bounds.second, 1.0));
    }

    // Same for gases!
    for (int g = 0; g < 10; ++g) {
      auto gas_bounds = rov.gas_bounds(g);
      REQUIRE(FloatingPoint<Real>::equiv(gas_bounds.first, 0.0));
      REQUIRE(FloatingPoint<Real>::equiv(gas_bounds.second, 1.0));
    }
  }

  SECTION("contains") {
    // By default, a region of validity accepts zero data, so the above
    // prognostics and atmosphere should pass muster.
    RegionOfValidity rov;
    rov.init(config);
    int valid = 0;
    Kokkos::parallel_reduce(
        "rov.contains 1", 1,
        KOKKOS_LAMBDA(const int, int& v) { v = rov.contains(atm, progs); },
        valid);
    REQUIRE(valid);

    // Suppose now that we require data to be positive!
    rov.set_interstitial_aerosol_bounds(0, 1e4, 1e11);
    rov.temp_bounds.first = 273.0;
    Kokkos::parallel_reduce(
        "rov.contains 2", 1,
        KOKKOS_LAMBDA(const int, int& v) { v = rov.contains(atm, progs); },
        valid);
    REQUIRE(not valid);

    // Fix the data and make sure it passes the next time.
    Kokkos::parallel_for(
        "rov.contains 3", num_vert_packs, KOKKOS_LAMBDA(const int k) {
          progs.interstitial_aerosols(0, k) = 1e-14;
          atm.temperature(k) = 325.0;
          atm.pressure(k) = 1e4;
        });
    Kokkos::parallel_reduce(
        "rov.contains 4", 1,
        KOKKOS_LAMBDA(const int, int& v) { v = rov.contains(atm, progs); },
        valid);
    REQUIRE(valid);
  }

  SECTION("intersection") {
    // Create two regions of validity and then intersect them.
    RegionOfValidity rov1;
    rov1.temp_bounds = {200, 300};
    rov1.rel_hum_bounds = {0.1, 0.9};
    rov1.set_interstitial_aerosol_bounds(0, 1e-8, 1e-1);
    rov1.set_interstitial_aerosol_bounds(1, 1e-5, 1e-3);

    RegionOfValidity rov2;
    rov2.temp_bounds = {250, 350};
    rov2.rel_hum_bounds = {0.2, 0.95};
    rov2.set_interstitial_aerosol_bounds(1, 1e-4, 1e-2);
    rov2.set_interstitial_aerosol_bounds(2, 1e-7, 1e-3);

    auto rov3 = RegionOfValidity::intersection(rov1, rov2);
    REQUIRE(FloatingPoint<Real>::equiv(rov3.temp_bounds.first, 250.0));
    REQUIRE(FloatingPoint<Real>::equiv(rov3.temp_bounds.second, 300.0));
    REQUIRE(FloatingPoint<Real>::equiv(rov3.rel_hum_bounds.first, 0.2));
    REQUIRE(FloatingPoint<Real>::equiv(rov3.rel_hum_bounds.second, 0.9));

    auto aero_bounds0 = rov3.interstitial_aerosol_bounds(0);
    REQUIRE(FloatingPoint<Real>::equiv(aero_bounds0.first, 1e-8));
    REQUIRE(FloatingPoint<Real>::equiv(aero_bounds0.second, 1e-1));

    auto aero_bounds1 = rov3.interstitial_aerosol_bounds(1);
    REQUIRE(FloatingPoint<Real>::equiv(aero_bounds1.first, 1e-4));
    REQUIRE(FloatingPoint<Real>::equiv(aero_bounds1.second, 1e-3));

    auto aero_bounds2 = rov3.interstitial_aerosol_bounds(2);
    REQUIRE(FloatingPoint<Real>::equiv(aero_bounds2.first, 1e-7));
    REQUIRE(FloatingPoint<Real>::equiv(aero_bounds2.second, 1e-3));
  }
}
