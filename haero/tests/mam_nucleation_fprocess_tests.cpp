#include "haero/model.hpp"
#include "haero/floating_point.hpp"
#include "haero/processes/mam_nucleation_fprocess.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

// These tests exercise our transplant of the MAM nucleation process.
TEST_CASE("mam_nucleation_fprocess", "mam_nucleation_fprocess") {

  // We create a phony model to be used for these tests.
  auto aero_config = create_mam4_modal_aerosol_config();
  int num_levels = 72;
  auto* model = Model::ForUnitTests(aero_config, num_levels);
  int num_gases = aero_config.h_gas_species.size();
  int num_modes = aero_config.h_aerosol_modes.size();

  // Set up some prognosics aerosol data views‥
  int num_aero_populations = model->num_aerosol_populations();
  Kokkos::View<PackType**> int_aerosols("interstitial aerosols",
                                        num_aero_populations, num_levels);
  Kokkos::View<PackType**> cld_aerosols("cloudborne aerosols",
                                        num_aero_populations, num_levels);
  Kokkos::View<PackType**> gases("gases", num_gases, num_levels);
  Kokkos::View<PackType**> modal_concs("modal number concs", num_modes,
                                       num_levels);

  // Set up atmospheric data and populate it with some views.
  Kokkos::View<PackType*> temp("temperature", num_levels);
  Kokkos::View<PackType*> press("pressure", num_levels);
  Kokkos::View<PackType*> rel_hum("relative humidity", num_levels);
  Kokkos::View<PackType*> ht("height", num_levels+1);
  Real pblh = 100.0;
  auto* atm = new Atmosphere(num_levels, temp, press, rel_hum, ht, pblh);

  // Test basic construction.
  SECTION("construct") {
    auto* process = new MAMNucleationFProcess();
    REQUIRE(process->type() == haero::NucleationProcess);
    REQUIRE(process->name() == "MAMNucleationFProcess (Fortran prognostic NucleationProcess)");
    delete process;
  }

  // Test process initialization.
  SECTION("init_process") {
    auto* process = new MAMNucleationFProcess();
    process->init(aero_config);
    delete process;
  }

  // Test process tendencies.
  SECTION("nucleate_without_existing_aerosols") {
    auto* process = new MAMNucleationFProcess();
    process->init(aero_config);

    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    auto* progs = model->create_prognostics(int_aerosols, cld_aerosols, gases,
                                            modal_concs);
    auto* diags = model->create_diagnostics();
    auto* tends = new Tendencies(*progs);

    // Set initial conditions.

    // atmospheric state
    Real h0 = 3e3, dz = h0/num_levels;
    for (int k = 0; k < num_levels; ++k) {
      temp(k) = 273.0;
      press(k) = 1e5;
      rel_hum(k) = 0.95;
      ht(k) = h0 - k*dz;
    }

    // aerosols (none)
    for (int p = 0; p < aero_config.num_aerosol_populations; ++p) {
      for (int k = 0; k < num_levels; ++k) {
        int_aerosols(p, k) = 0.0;
      }
    }

    // gases
    int h2so4_index = aero_config.gas_index("H2SO4");
    printf("h2so4 index: %d\n", h2so4_index);
    for (int k = 0; k < num_levels; ++k) {
      gases(h2so4_index, k) = 1e-13;
    }
    for (int g = 0; g < num_gases; ++g) {
      if (g != h2so4_index) {
        for (int k = 0; k < num_levels; ++k) {
          gases(g, k) = 0.0;
        }
      }
    }

    // Now compute the tendencies by running the process.
    Real t = 0.0, dt = 30.0;
    process->run(aero_config, t, dt, *progs, *atm, *diags, *tends);

    // --------------------------------------------------
    // Check the tendencies to make sure they make sense.
    // --------------------------------------------------

    // SO4 nucleates within the aitken mode. All other tendencies are zero.
    const auto aero_tends = tends->interstitial_aerosols();
    int aitken_index = aero_config.aerosol_mode_index("aitken");
    int aitken_so4_index = aero_config.aerosol_species_index(aitken_index, "SO4");
    int so4_pop_index = aero_config.population_index(aitken_index, aitken_so4_index);
    for (int p = 0; p < aero_config.num_aerosol_populations; ++p) {
      if (p == so4_pop_index) {
        for (int k = 0; k < num_levels; ++k) {
          REQUIRE(int_aerosols(p, k)[0] > 0.0);
        }
      } else {
        for (int k = 0; k < num_levels; ++k) {
          REQUIRE(int_aerosols(p, k)[0] == 0.0);
        }
      }
    }

    // The tendency for H2SO4 should be negative, and the rest should be zero.
    const auto gas_tends = tends->gases();
    for (int k = 0; k < num_levels; ++k) {
      REQUIRE(gas_tends(h2so4_index, k)[0] <= 0.0);
    }
    for (int g = 0; g < num_gases; ++g) {
      if (g != h2so4_index) {
        for (int k = 0; k < num_levels; ++k) {
          REQUIRE(gas_tends(g, k)[0] == 0.0);
        }
      }
    }

    // Clean up.
    delete progs;
    delete diags;
    delete tends;
    delete process;
  }

  delete model;
}

