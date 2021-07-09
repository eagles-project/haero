#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/floating_point.hpp"
#include "haero/region_of_validity.hpp"

using namespace haero;

TEST_CASE("region_of_validity", "") {
  // Create some prognostics and atmosphere state data.
  int num_levels = 72;
  int num_vert_packs = num_levels / HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels) {
    num_vert_packs++;
  }
  SpeciesColumnView int_aerosols("interstitial aerosols", 1, num_vert_packs);
  SpeciesColumnView cld_aerosols("cloudborne aerosols", 1, num_vert_packs);
  int num_modes = 1;
  ModalColumnView int_num_concs("interstitial number concs", num_modes,
                                num_vert_packs);
  ModalColumnView cld_num_concs("cloudborne number concs", num_modes,
                                num_vert_packs);
  int num_gases = 1;
  SpeciesColumnView gases("gases", num_gases, num_vert_packs);

  Prognostics progs(num_modes, {1}, num_gases, num_levels, int_aerosols,
                    cld_aerosols, gases, int_num_concs, cld_num_concs);

  ColumnView temp("temperature", num_vert_packs),
      press("pressure", num_vert_packs),
      rel_hum("relative humidity", num_vert_packs),
      ht("height", num_vert_packs + 1),
      pdel("pressure thickness", num_vert_packs);
  const Real pblh = 1000.0;
  Atmosphere atm(num_levels, temp, press, rel_hum, ht, pdel, pblh);

  // Zero the state data.
  for (int p = 0; p < progs.num_aerosol_populations(); ++p) {
    for (int k = 0; k < num_vert_packs; ++k) {
      int_aerosols(p, k) = 0;
      cld_aerosols(p, k) = 0;
    }
  }
  for (int m = 0; m < num_modes; ++m) {
    for (int k = 0; k < num_vert_packs; ++k) {
      int_num_concs(m, k) = 0;
      cld_num_concs(m, k) = 0;
    }
  }
  for (int g = 0; g < num_gases; ++g) {
    for (int k = 0; k < num_vert_packs; ++k) {
      gases(g, k) = 0;
    }
  }

  for (int k = 0; k < num_vert_packs; ++k) {
    atm.temperature(k) = 0;
    atm.pressure(k) = 0;
    atm.relative_humidity(k) = 0;
  }

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
      auto int_aero_mmr_bounds = rov.interstitial_aerosol_mmr_bounds(p);
      REQUIRE(FloatingPoint<Real>::equiv(int_aero_mmr_bounds.first, 0.0));
      REQUIRE(FloatingPoint<Real>::equiv(int_aero_mmr_bounds.second, 1.0));

      auto cld_aero_mmr_bounds = rov.cloud_aerosol_mmr_bounds(p);
      REQUIRE(FloatingPoint<Real>::equiv(cld_aero_mmr_bounds.first, 0.0));
      REQUIRE(FloatingPoint<Real>::equiv(cld_aero_mmr_bounds.second, 1.0));

      auto int_n_bounds = rov.cloud_aerosol_num_conc_bounds(p);
      REQUIRE(FloatingPoint<Real>::equiv(int_n_bounds.first, 0.0));
      REQUIRE(FloatingPoint<Real>::equiv(int_n_bounds.second, 1e20));

      auto cld_n_bounds = rov.cloud_aerosol_num_conc_bounds(p);
      REQUIRE(FloatingPoint<Real>::equiv(cld_n_bounds.first, 0.0));
      REQUIRE(FloatingPoint<Real>::equiv(cld_n_bounds.second, 1e20));
    }

    // Same for gases!
    for (int g = 0; g < 10; ++g) {
      auto gas_mmr_bounds = rov.gas_mmr_bounds(g);
      REQUIRE(FloatingPoint<Real>::equiv(gas_mmr_bounds.first, 0.0));
      REQUIRE(FloatingPoint<Real>::equiv(gas_mmr_bounds.second, 1.0));
    }
  }

  SECTION("contains") {
    // By default, a region of validity accepts zero data, so the above
    // prognostics and atmosphere should pass muster.
    RegionOfValidity rov;
    REQUIRE(rov.contains(progs));
    REQUIRE(rov.contains(atm));

    // Suppose now that we require data to be positive!
    rov.set_interstitial_aerosol_mmr_bounds(0, 1e4, 1e11);
    REQUIRE(not rov.contains(progs));
    rov.temp_bounds.first = 273.0;
    REQUIRE(not rov.contains(atm));

    // Fix the data and make sure it passes the next time.
    for (int k = 0; k < num_vert_packs; ++k) {
      progs.interstitial_aerosols(0, k) = 1e4;
      atm.temperature(k) = 273.0;
    }
    REQUIRE(rov.contains(progs));
    REQUIRE(rov.contains(atm));
  }

  SECTION("intersection") {
    // Create two regions of validity and then intersect them.
    RegionOfValidity rov1;
    rov1.temp_bounds = {200, 300};
    rov1.rel_hum_bounds = {0.1, 0.9};
    rov1.set_interstitial_aerosol_mmr_bounds(0, 1e-8, 1e-1);
    rov1.set_interstitial_aerosol_mmr_bounds(1, 1e-5, 1e-3);

    RegionOfValidity rov2;
    rov2.temp_bounds = {250, 350};
    rov2.rel_hum_bounds = {0.2, 0.95};
    rov2.set_interstitial_aerosol_mmr_bounds(1, 1e-4, 1e-2);
    rov2.set_interstitial_aerosol_mmr_bounds(2, 1e-7, 1e-3);

    auto rov3 = RegionOfValidity::intersection(rov1, rov2);
    REQUIRE(FloatingPoint<Real>::equiv(rov3.temp_bounds.first, 250.0));
    REQUIRE(FloatingPoint<Real>::equiv(rov3.temp_bounds.second, 300.0));
    REQUIRE(FloatingPoint<Real>::equiv(rov3.rel_hum_bounds.first, 0.2));
    REQUIRE(FloatingPoint<Real>::equiv(rov3.rel_hum_bounds.second, 0.9));

    auto ralf_bounds0 = rov1.interstitial_aerosol_mmr_bounds(0);
    auto aero_bounds0 = rov3.interstitial_aerosol_mmr_bounds(0);
    REQUIRE(FloatingPoint<Real>::equiv(aero_bounds0.first, 1e-8));
    REQUIRE(FloatingPoint<Real>::equiv(aero_bounds0.second, 1e-1));

    auto aero_bounds1 = rov3.interstitial_aerosol_mmr_bounds(1);
    REQUIRE(FloatingPoint<Real>::equiv(aero_bounds1.first, 1e-4));
    REQUIRE(FloatingPoint<Real>::equiv(aero_bounds1.second, 1e-3));

    auto aero_bounds2 = rov3.interstitial_aerosol_mmr_bounds(2);
    REQUIRE(FloatingPoint<Real>::equiv(aero_bounds2.first, 1e-7));
    REQUIRE(FloatingPoint<Real>::equiv(aero_bounds2.second, 1e-3));
  }
}
