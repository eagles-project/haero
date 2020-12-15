#include "haero/model.hpp"
#include "prog_process_stub.hpp"
#include "diag_process_stub.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

// These tests demonstrate our minimal Fortran-backed prognostic process stub.
TEST_CASE("prog_process_stub", "") {

  // We create a phony model to be used for these tests.
  std::vector<Mode> modes = create_mam4_modes();
  std::vector<Species> aero_species, gas_species;
  std::map<std::string, std::vector<std::string> > mode_species;
  int num_columns = 10;
  int num_levels = 72;
  auto* model = Model::ForUnitTests(modes, aero_species, mode_species,
                  gas_species, num_columns, num_levels);

  // Rate of decay from cloudborne to interstitial aerosols.
  Real decay_rate = -0.05;

  // Test basic construction.
  SECTION("construct") {
    auto* stub = new ProgProcessStub(decay_rate);
    REQUIRE(stub->type() == haero::ActivationProcess);
    REQUIRE(stub->name() == "Prognostic process stub (Fortran)");
    delete stub;
  }

  // Test process initialization.
  SECTION("init_process") {
    auto* stub = new ProgProcessStub(decay_rate);
    stub->init(*model);
    delete stub;
  }

  // Test process tendencies.
  SECTION("tendencies") {
    auto* stub = new ProgProcessStub(decay_rate);
    stub->init(*model);

    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    auto* progs = model->create_prognostics();
    auto* diags = model->create_diagnostics();
    auto* tends = new Tendencies(*progs);

    // Set initial conditions.
    Real N0 = 1e6;
    auto& modal_num_densities = progs->modal_num_densities();
    for (int m = 0; m < progs->num_aerosol_modes(); ++m) {
      auto& cld_aerosols = progs->cloudborne_aerosols(m);
      auto& int_aerosols = progs->interstitial_aerosols(m);
      int num_aero_species = progs->num_aerosol_species(m);
      for (int i = 0; i < progs->num_columns(); ++i) {
        for (int k = 0; k < progs->num_levels(); ++k) {

          // Set aerosol mix fractions. We start with 100% clouds for simplicity.
          for (int s = 0; s < num_aero_species; ++s) {
            cld_aerosols(i, k, s) = 1.0 / num_aero_species;
            int_aerosols(i, k, s) = 0.0;
          }

          // Set modal number densities.
          modal_num_densities(m, i, k) = N0;
        }
      }
    }

    // Set gas mole fractions.
    auto& gas_mole_fracs = progs->gas_mole_fractions();
    int num_gas_species = progs->num_gas_species();
    for (int i = 0; i < progs->num_columns(); ++i) {
      for (int k = 0; k < progs->num_levels(); ++k) {
        for (int s = 0; s < num_gas_species; ++s) {
          gas_mole_fracs(i, k, s) = 1.0 / num_gas_species;
        }
      }
    }

    // Now compute the tendencies by running the process.
    Real t = 0.0, dt = 0.01;
    stub->run(*model, t, dt, *progs, *diags, *tends);

    // Clean up.
    delete progs;
    delete diags;
    delete tends;
    delete stub;
  }

  delete model;
}

// These tests demonstrate our minimal Fortran-backed diagnostic process stub.
TEST_CASE("diag_process_stub", "") {

  // We create a phony model to be used for these tests.
  Parameterizations params;
  std::vector<Mode> modes;
  std::vector<Species> aero_species, gas_species;
  std::map<std::string, std::vector<std::string> > mode_species;
  int num_columns = 10;
  int num_levels = 72;
  auto* model = new Model(params, modes, aero_species, mode_species,
                          gas_species, num_columns, num_levels);

  // Test basic construction.
  SECTION("construct") {
    auto* stub = new DiagProcessStub();
    REQUIRE(stub->type() == haero::WaterUptakeProcess);
    REQUIRE(stub->name() == "Diagnostic process stub (Fortran)");
    delete stub;
  }

  // Test preparation of diagnostic variables.
  SECTION("prepare_diagnostic_vars") {
    Diagnostics diags(10, 72, {1, 2}, 1);
    auto* stub = new DiagProcessStub();
    stub->prepare(diags);
    REQUIRE(diags.has_modal_var("pressure"));
    delete stub;
  }

  // Test process initialization.
  SECTION("init_process") {
    // We create a phony model to be passed to the init method.
    auto* stub = new DiagProcessStub();
    stub->init(*model);
    delete stub;
  }

  delete model;
}

