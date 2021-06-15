#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "mam_calcsize_test_bridge.hpp"
#include "haero/processes/mam_calcsize_process.hpp"

using namespace haero;

// These tests demonstrate our minimal Fortran-backed prognostic process stub.
TEST_CASE("mam_calcsize_run", "") {
  static_assert(HAERO_PACK_SIZE == 1,
                "Fortran not supported for HAERO_PACK_SIZE != 1.");

  // We create a phony model to be used for these tests.
  auto modes        = create_mam4_modes();
  auto aero_species = create_mam4_aerosol_species();
  auto gas_species  = create_mam4_gas_species();
  auto mode_species = create_mam4_mode_species();
  int num_levels = 72;
  ModalAerosolConfig aero_config(modes, aero_species, mode_species,
                                 gas_species);
  auto* model = Model::ForUnitTests(aero_config, num_levels);

  // Set up some prognosics aerosol data views‥
  int num_aero_populations = model->num_aerosol_populations();
  Kokkos::View<PackType**> int_aerosols("interstitial aerosols",
                                        num_aero_populations, num_levels);
  Kokkos::View<PackType**> cld_aerosols("cloudborne aerosols",
                                        num_aero_populations, num_levels);
  int num_gases = gas_species.size();
  Kokkos::View<PackType**> gases("gases", num_gases, num_levels);
  int num_modes = modes.size();
  Kokkos::View<PackType**> int_num_concs("interstitial number concs", num_modes,
                                         num_levels);
  Kokkos::View<PackType**> cld_num_concs("cloud borne number concs", num_modes,
                                         num_levels);

  // Set up atmospheric data and populate it with some views. It's not
  // important for this data to be valid, since it's unused by these stubs.
  Kokkos::View<PackType*> temp("temperature", num_levels);
  Kokkos::View<PackType*> press("pressure", num_levels);
  Kokkos::View<PackType*> rel_hum("relative humidity", num_levels);
  Kokkos::View<PackType*> pdel("hydrostatic_dp", num_levels);
  Kokkos::View<PackType*> ht("height", num_levels + 1);
  Real pblh = 100.0;
  auto* atm = new Atmosphere(num_levels, temp, press, rel_hum, ht, pdel, pblh);



  // Test process tendencies.
  SECTION("tendencies") {
    //auto* stub = new MAMCalcsizeProcess();

    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    auto* progs = model->create_prognostics(int_aerosols, cld_aerosols, gases,
                                            int_num_concs, cld_num_concs);
    auto* diags = model->create_diagnostics();
    auto* tends = new Tendencies(*progs);

    // Set initial conditions.

    // Set aerosol mass mix fractions. We start with 100% clouds
    // for simplicity.
    for (int p = 0; p < num_aero_populations; ++p) {
      for (int k = 0; k < num_levels; ++k) {
        cld_aerosols(p, k) = 1.0 / num_aero_populations;
        int_aerosols(p, k) = 0.0;
      }
    }

    // Set modal number concentrations.
    Real n0 = 1e6;
    for (int m = 0; m < num_modes; ++m) {
      for (int k = 0; k < num_levels; ++k) {
        int_num_concs(m, k) = n0;
        cld_num_concs(m, k) = 0;
      }
    }

    // Set gas mass mixing rations.
    for (int g = 0; g < num_gases; ++g) {
      for (int k = 0; k < num_levels; ++k) {
        gases(g, k) = 1.0 / num_gases;
      }
    }

    // Now compute the tendencies by running the process.
    Real t = 0.0, dt = 0.01;
    run_bridge(model->modal_aerosol_config(), t, dt, *progs, *atm, *diags,
               *tends);
    /*
    // --------------------------------------------------
    // Check the tendencies to make sure they make sense.
    // --------------------------------------------------

    // Cloudborne aerosol mix fractions have negative tendencies, interstitial
    // mix fractions have positive tendencies, and their sums are zero.
    auto dqdt_c = tends->cloud_aerosols;
    auto dqdt_i = tends->interstitial_aerosols;
    for (int p = 0; p < num_aero_populations; ++p) {
      for (int k = 0; k < num_levels; ++k) {
        REQUIRE(dqdt_c(p, k)[0] < 0.0);
        REQUIRE(dqdt_i(p, k)[0] > 0.0);
        Real sum = dqdt_c(p, k)[0] + dqdt_i(p, k)[0];
        REQUIRE(FloatingPoint<Real>::equiv(sum, 0.0));
      }
    }

    // Aerosol modal number concentrations are unchanged.
    auto dndt = tends->interstitial_num_concs;
    for (int m = 0; m < num_modes; ++m) {
      for (int k = 0; k < num_levels; ++k) {
        REQUIRE(FloatingPoint<Real>::equiv(dndt(m, k)[0], 0.0));
      }
    }

    // Gas mix ratios are unchanged.
    const auto& dqdt_g = tends->gases;
    for (int g = 0; g < num_gases; ++g) {
      for (int k = 0; k < num_levels; ++k) {
        REQUIRE(FloatingPoint<Real>::equiv(dqdt_g(g, k)[0], 0.0));
      }
    }

    // Clean up.
    delete progs;
    delete diags;
    delete tends;
    delete stub;*/
  }

  delete atm;
  delete model;
}

TEST_CASE("compute_diameter", "mam_calcsize_fprocess") {

  static constexpr Real vol2num = 2;
  //compute diameter
  //MAMCalcsizeProcess::compute_diameter(vol2num);
  Real diameter = compute_diameter_bridge(vol2num);

  REQUIRE(diameter==3);
}
