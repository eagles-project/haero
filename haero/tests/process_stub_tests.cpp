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
  Parameterizations params;
  std::vector<Mode> modes = create_mam4_modes();
  std::vector<Species> aero_species, gas_species;
  std::map<std::string, std::vector<std::string> > mode_species;
  int num_columns = 10;
  int num_levels = 72;
  auto* model = new Model(params, modes, aero_species, mode_species,
                          gas_species, num_columns, num_levels);

  // Test basic construction.
  SECTION("construct") {
    auto* stub = new ProgProcessStub();
    REQUIRE(stub->type() == haero::ActivationProcess);
    REQUIRE(stub->name() == "Prognostic process stub (Fortran)");
    delete stub;
  }

  // Test process initialization.
  SECTION("init_process") {
    auto* stub = new ProgProcessStub();
    stub->init(*model);
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

