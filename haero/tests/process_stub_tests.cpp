#include "haero/model.hpp"
#include "haero/floating_point.hpp"
#include "prog_process_stub.hpp"
#include "diag_process_stub.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

// These tests demonstrate our minimal Fortran-backed prognostic process stub.
TEST_CASE("prog_process_stub", "") {

  // We create a phony model to be used for these tests.
  auto modes = create_mam4_modes();
  auto aero_species = create_mam4_aerosol_species(),
       gas_species = create_mam4_gas_species();
  auto mode_species = create_mam4_mode_species();
  int num_columns = 10;
  int num_levels = 72;
  auto* model = Model::ForUnitTests(modes, aero_species, mode_species,
                                    gas_species, num_columns, num_levels);

  // Set up atmospheric data and populate it with some views. It's not
  // important for this data to be valid, since it's unused by these stubs.
  using PackType = Atmosphere::PackType;
  Kokkos::View<PackType**> temp, press, rel_hum, ht;
  auto* atm = new Atmosphere(model->num_columns(), model->num_levels(),
                             temp, press, rel_hum, ht);

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
    stub->run(*model, t, dt, *progs, *atm, *diags, *tends);

    // --------------------------------------------------
    // Check the tendencies to make sure they make sense.
    // --------------------------------------------------

    // Cloudborne aerosol mix fractions have negative tendencies, interstitial
    // mix fractions have positive tendencies, and their sums are zero.
    for (int m = 0; m < progs->num_aerosol_modes(); ++m) {
      auto& dqdt_c = tends->cloudborne_aerosols(m);
      auto& dqdt_i = tends->interstitial_aerosols(m);
      int num_aero_species = progs->num_aerosol_species(m);
      for (int i = 0; i < progs->num_columns(); ++i) {
        for (int k = 0; k < progs->num_levels(); ++k) {
          for (int s = 0; s < num_aero_species; ++s) {
            REQUIRE(dqdt_c(i, k, s)[0] < 0.0);
            REQUIRE(dqdt_i(i, k, s)[0] > 0.0);
            Real sum = dqdt_c(i, k, s)[0] + dqdt_i(i, k, s)[0];
            REQUIRE(FloatingPoint<Real>::equiv(sum, 0.0));
          }
        }
      }
    }

    // Aerosol modal number densities are unchanged.
    auto& dndt = tends->modal_num_densities();
    for (int m = 0; m < progs->num_aerosol_modes(); ++m) {
      for (int i = 0; i < progs->num_columns(); ++i) {
        for (int k = 0; k < progs->num_levels(); ++k) {
          REQUIRE(FloatingPoint<Real>::equiv(dndt(m, i, k)[0], 0.0));
        }
      }
    }

    // Gas mole fractions are unchanged.
    const auto& dqdt_g = tends->gas_mole_fractions();
    for (int i = 0; i < progs->num_columns(); ++i) {
      for (int k = 0; k < progs->num_levels(); ++k) {
        for (int s = 0; s < num_gas_species; ++s) {
          REQUIRE(FloatingPoint<Real>::equiv(dqdt_g(i, k, s)[0], 0.0));
        }
      }
    }


    // Clean up.
    delete progs;
    delete diags;
    delete tends;
    delete stub;
  }

  delete atm;
  delete model;
}

// These tests demonstrate our minimal Fortran-backed diagnostic process stub.
TEST_CASE("diag_process_stub", "") {

  // We create a phony model to be used for these tests.
  auto modes = create_mam4_modes();
  auto aero_species = create_mam4_aerosol_species(),
       gas_species = create_mam4_gas_species();
  auto mode_species = create_mam4_mode_species();
  int num_columns = 10;
  int num_levels = 72;
  auto* model = Model::ForUnitTests(modes, aero_species, mode_species,
                                    gas_species, num_columns, num_levels);

  // Set up atmospheric data and populate it with some views. It's not
  // important for this data to be valid, since it's unused by these stubs.
  using PackType = Atmosphere::PackType;
  Kokkos::View<PackType**> temp, press, rel_hum, ht;
  auto* atm = new Atmosphere(model->num_columns(), model->num_levels(),
                             temp, press, rel_hum, ht);

  // Test basic construction.
  SECTION("construct") {
    auto* stub = new DiagProcessStub();
    REQUIRE(stub->type() == haero::WaterUptakeProcess);
    REQUIRE(stub->name() == "Diagnostic process stub (Fortran)");
    delete stub;
  }

  // Test preparation of diagnostic variables.
  SECTION("prepare_diagnostic_vars") {
    std::vector<int> num_aero_species(modes.size());
    for (int m = 0; m < modes.size(); ++m) {
      num_aero_species[m] = mode_species[modes[m].name].size();
    }
    int num_gas_species = gas_species.size();
    Diagnostics diags(num_columns, num_levels, num_aero_species, num_gas_species);
    auto* stub = new DiagProcessStub();
    stub->prepare(diags);
    REQUIRE(diags.has_var("temperature"));
    REQUIRE(diags.has_gas_var("pressure"));
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

  SECTION("update_diagnostics") {
    auto* stub = new DiagProcessStub();
    stub->init(*model);

    // Initialize prognostic and diagnostic variables.
    auto* progs = model->create_prognostics();
    auto* diags = model->create_diagnostics();

    // Set initial conditions.
    Real N0 = 1e16;
    auto& modal_num_densities = progs->modal_num_densities();
    for (int m = 0; m < progs->num_aerosol_modes(); ++m) {
      auto& cld_aerosols = progs->cloudborne_aerosols(m);
      auto& int_aerosols = progs->interstitial_aerosols(m);
      int num_aero_species = progs->num_aerosol_species(m);
      for (int i = 0; i < progs->num_columns(); ++i) {
        for (int k = 0; k < progs->num_levels(); ++k) {

          // Set aerosol mix fractions (50% cloudborne, 50% interstitial)
          for (int s = 0; s < num_aero_species; ++s) {
            cld_aerosols(i, k, s) = 0.5 / num_aero_species;
            int_aerosols(i, k, s) = 0.5 / num_aero_species;
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

    // Set up required diagnostic variables.
    stub->prepare(*diags);

    // Set the atmospheric temperature.
    {
      auto& T = diags->var("temperature");
      for (int i = 0; i < progs->num_columns(); ++i) {
        for (int k = 0; k < progs->num_levels(); ++k) {
          T(i, k) = 273.15;
        }
      }
    }

    // Make a copy of the temperature field.
    auto T0 = diags->var("temperature");

    // Now update diagnostics at time 0.
    Real t = 0.0;
    stub->update(*model, t, *progs, *atm, *diags);

    // -------------------------------------------------
    // Make sure the temperature field was not affected.
    // -------------------------------------------------
    const auto& T = diags->var("temperature");
    for (int i = 0; i < progs->num_columns(); ++i) {
      for (int k = 0; k < progs->num_levels(); ++k) {
        REQUIRE(FloatingPoint<Real>::equiv(T(i, k)[0], T0(i, k)[0]));
      }
    }

    // ---------------------------------------
    // Check the diagnostic partial pressures.
    // ---------------------------------------
    const auto& p_m = diags->modal_var("pressure");
    for (int m = 0; m < progs->num_aerosol_modes(); ++m) {
      for (int i = 0; i < progs->num_columns(); ++i) {
        for (int k = 0; k < progs->num_levels(); ++k) {
          REQUIRE(p_m(m, i, k)[0] > 0.0);
        }
      }
    }

    const auto& p_g = diags->gas_var("pressure");
    for (int i = 0; i < progs->num_columns(); ++i) {
      for (int k = 0; k < progs->num_levels(); ++k) {
        for (int s = 0; s < num_gas_species; ++s) {
          REQUIRE(p_g(i, k, s)[0] > 0.0);
        }
      }
    }

    // Clean up.
    delete progs;
    delete diags;
    delete stub;
  }

  delete atm;
  delete model;
}

