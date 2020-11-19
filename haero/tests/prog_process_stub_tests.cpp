#include "haero/model.hpp"
#include "haero/processes/prog_process_stub.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

// These tests demonstrate our minimal Fortran-backed prognostic process stub.
TEST_CASE("prog_process_stub", "") {

  // We create a phony model to be used for these tests.
  Parameterizations params;
  std::vector<Mode> modes;
  std::vector<Species> aero_species, gas_species;
  std::map<std::string, std::vector<std::string> > mode_species;
  auto* model = new Model(params, modes, aero_species, mode_species,
                          gas_species, 10, 72);

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

