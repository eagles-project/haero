#include "haero/model.hpp"
#include "haero/floating_point.hpp"
#include "haero/processes/mam_nucleation_fprocess.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

// These tests demonstrate our minimal Fortran-backed prognostic process stub.
TEST_CASE("mam_nucleation_fprocess", "") {

  // We create a phony model to be used for these tests.
  auto modal_aero_config = create_mam4_modal_aerosol_config();
  int num_levels = 72;
  auto* model = Model::ForUnitTests(modal_aero_config, num_levels);
  int num_gases = modal_aero_config.h_gas_species.size();
  int num_modes = modal_aero_config.h_aerosol_modes.size();

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
    process->init(modal_aero_config);
    delete process;
  }

  // Test process tendencies.
  SECTION("run_process") {
    auto* process = new MAMNucleationFProcess();
    process->init(modal_aero_config);

    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    auto* progs = model->create_prognostics(int_aerosols, cld_aerosols, gases,
                                            modal_concs);
    auto* diags = model->create_diagnostics();
    auto* tends = new Tendencies(*progs);

    // Set initial conditions.
    // FIXME

    // Now compute the tendencies by running the process.
    Real t = 0.0, dt = 0.01;
    process->run(modal_aero_config, t, dt, *progs, *atm, *diags, *tends);

    // --------------------------------------------------
    // Check the tendencies to make sure they make sense.
    // --------------------------------------------------
    // TODO

    // Clean up.
    delete progs;
    delete diags;
    delete tends;
    delete process;
  }

  delete model;
}

