#include "haero/model.hpp"
#include "haero/floating_point.hpp"
#include "haero/processes/mam_nucleation_process.hpp"
#include "mam_nucleation_test_bridge.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace haero;

TEST_CASE("ternary_nuc_merik2007", "mam_nucleation_process") {
  // Define a pseudo-random generator [0-1] that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  const unsigned p0  = 987659;
  const unsigned p1  =  12373;
  long unsigned seed =  54319;
  auto random = [&]() {
    seed =  (seed * p1)%p0;
    return double(seed)/p0;
  };
  for (int i=0; i<1000; ++i) {
    const double  t=  235 +   60*random();  // range 235-295
    const double rh= 0.05 +   .9*random();  // range .05-.95
    const double c2= 5.e4 + 1.e8*random();  // range 5x10^4 - 10^9 
    const double c3=  0.1 +  999*random();  // range 0.1 - 1000
    double j_log_cpp;
    double ntot_cpp; 
    double naci_cpp;
    double namm_cpp; 
    double r_cpp;
    double j_log_f90;
    double ntot_f90; 
    double naci_f90;
    double namm_f90; 
    double r_f90;
    MAMNucleationProcess::ternary_nuc_merik2007(t, rh, c2, c3, j_log_cpp, ntot_cpp, naci_cpp, namm_cpp, r_cpp);
    ternary_nuc_merik2007_bridge(t, rh, c2, c3, j_log_f90, ntot_f90, naci_f90, namm_f90, r_f90);
    REQUIRE(j_log_cpp == j_log_f90);
    REQUIRE(ntot_cpp  == ntot_f90);
    REQUIRE(naci_cpp  == naci_f90);
    REQUIRE(namm_cpp  == namm_f90);
    REQUIRE(r_cpp     == r_f90);
  }
}

// These tests exercise our transplant of the MAM nucleation process.
TEST_CASE("virtual_process_test", "mam_nucleation_process") {

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
    auto* process = new MAMNucleationProcess();
    REQUIRE(process->type() == haero::NucleationProcess);
    REQUIRE(process->name() == "MAMNucleationProcess");
    delete process;
  }

  // Test process initialization.
  SECTION("init_process") {
    auto* process = new MAMNucleationProcess();
    process->init(aero_config);
    delete process;
  }

  // Test process tendencies.
  SECTION("nucleate_without_existing_aerosols") {
    auto* process = new MAMNucleationProcess();
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

    // Clean up.
    delete progs;
    delete diags;
    delete tends;
    delete process;
  }

  delete model;
}

