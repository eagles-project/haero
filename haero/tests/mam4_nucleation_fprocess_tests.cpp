#include "haero/model.hpp"
#include "haero/floating_point.hpp"
#include "haero/processes/mam4_nucleation_fprocess.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

// These tests demonstrate our minimal Fortran-backed prognostic process stub.
TEST_CASE("mam4_nucleation_fprocess", "") {

  // We create a phony model to be used for these tests.
  auto modes = create_mam4_modes();
  auto aero_species = create_mam4_aerosol_species(),
       gas_species = create_mam4_gas_species();
  auto mode_species = create_mam4_mode_species();
  int num_levels = 72;
  SelectedProcesses selected_processes;
  selected_processes.nucleation = SelectedProcesses::MAM4FNucleation;
  auto* model = new Model(selected_processes, modes, aero_species,
                          mode_species, gas_species, num_levels);

  // Set up some prognosics aerosol data views‥
  int num_aero_populations = model->num_aerosol_populations();
  Kokkos::View<PackType**> int_aerosols("interstitial aerosols",
                                        num_aero_populations, num_levels);
  Kokkos::View<PackType**> cld_aerosols("cloudborne aerosols",
                                        num_aero_populations, num_levels);
  int num_gases = gas_species.size();
  Kokkos::View<PackType**> gases("gases", num_gases, num_levels);
  int num_modes = modes.size();
  Kokkos::View<PackType**> modal_concs("modal number concs", num_modes,
                                       num_levels);

  // Set up atmospheric data and populate it with some views.
  Kokkos::View<PackType*> temp("temperature", num_levels);
  Kokkos::View<PackType*> press("pressure", num_levels);
  Kokkos::View<PackType*> rel_hum("relative humidity", num_levels);
  Kokkos::View<PackType*> ht("height", num_levels+1);
  auto* atm = new Atmosphere(num_levels, temp, press, rel_hum, ht);

  // Test basic construction.
  SECTION("construct") {
    auto* process = new Mam4NucleationFProcess();
    REQUIRE(process->type() == haero::NucleationProcess);
    REQUIRE(process->name() == "MAM4 Nucleation (Fortran)");
    delete process;
  }

  // Test process initialization.
  SECTION("init_process") {
    auto* process = new Mam4NucleationFProcess();
    process->init(*model);
    delete process;
  }

  // Test process tendencies.
  SECTION("run_process") {
    auto* process = new Mam4NucleationFProcess();
    process->init(*model);

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
    process->run(*model, t, dt, *progs, *atm, *diags, *tends);

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

