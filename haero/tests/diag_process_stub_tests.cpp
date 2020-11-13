#include "haero/model.hpp"
#include "haero/processes/diag_process_stub.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

// These tests demonstrate our minimal Fortran-backed diagnostic process stub.
TEST_CASE("diag_process_stub", "") {

  // Test basic construction.
  SECTION("construct") {
    auto* stub = new DiagProcessStub();
    REQUIRE(stub->type() == haero::WaterUptakeProcess);
    REQUIRE(stub->name() == "Diagnostic process stub (Fortran)");
    delete stub;
  }

  // Test basic construction.
  SECTION("prepare_diagnostic_vars") {
    Diagnostics diags(10, 72, {1, 2}, 1);
    auto* stub = new DiagProcessStub();
    stub->prepare(diags);
    REQUIRE(diags.has_var("diag_stub_var"));
    REQUIRE(diags.has_modal_var("diag_stub_modal_var"));
    delete stub;
  }

  // Test process initialization.
  SECTION("init_process") {
    // We create a phony model to be passed to the init method.
    Parameterizations params;
    std::vector<Mode> modes;
    std::vector<Species> aero_species, gas_species;
    std::map<std::string, std::vector<std::string> > mode_species;
    auto* model = new Model(params, modes, aero_species, mode_species,
                            gas_species, 10, 72);
    auto* stub = new DiagProcessStub();
    stub->init(*model);
    delete stub;
    delete model;
  }
}

