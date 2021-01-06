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
  int num_columns = 10;
  int num_levels = 72;
  auto* model = Model::ForUnitTests(modes, aero_species, mode_species,
                                    gas_species, num_columns, num_levels);

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
    auto* progs = model->create_prognostics();
    auto* diags = model->create_diagnostics();
    auto* tends = new Tendencies(*progs);

    // Set initial conditions.
    // FIXME

    // Now compute the tendencies by running the process.
    Real t = 0.0, dt = 0.01;
    process->run(*model, t, dt, *progs, *diags, *tends);

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

