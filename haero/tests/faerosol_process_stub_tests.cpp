#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "faerosol_process_stub.hpp"
#include "haero/diagnostics.hpp"
#include "haero/floating_point.hpp"

using namespace haero;

// These tests demonstrate our minimal Fortran-backed prognostic process stub.
TEST_CASE("faerosol_process_stub", "") {
  static_assert(HAERO_PACK_SIZE == 1,
                "Fortran not supported for HAERO_PACK_SIZE != 1.");

  // We create a phony model to be used for these tests.
  auto aero_config = ModalAerosolConfig::create_mam4_config();
  int num_gases = aero_config.num_gases();
  int num_modes = aero_config.num_modes();
  int num_aero_populations = aero_config.num_aerosol_populations;
  int num_levels = 72;
  int num_vert_packs = PackInfo::num_packs(num_levels);
  int num_iface_packs = PackInfo::num_packs(num_levels+1);

  // Set up atmospheric data and populate it with some views. It's not
  // important for this data to be valid, since it's unused by these stubs.
  Kokkos::View<PackType*> temp("temperature", num_vert_packs);
  Kokkos::View<PackType*> press("pressure", num_vert_packs);
  Kokkos::View<PackType*> qv("vapor mixing ratio", num_vert_packs);
  Kokkos::View<PackType*> pdel("hydrostatic_dp", num_vert_packs);
  Kokkos::View<PackType*> ht("height", num_iface_packs);
  Real pblh = 100.0;
  auto* atm = new Atmosphere(num_levels, temp, press, qv, ht, pdel, pblh);

  // Rate of decay from cloudborne to interstitial aerosols.
  Real decay_rate = -0.05;

  // Test basic construction.
  SECTION("construct") {
    auto* stub = new FAerosolProcessStub();
    stub->set_param("decay_rate", decay_rate);
    REQUIRE(stub->name() == "Aerosol process stub (Fortran)");
    delete stub;
  }

  // Test process initialization.
  SECTION("init_process") {
    auto* stub = new FAerosolProcessStub();
    stub->set_param("decay_rate", decay_rate);
    stub->init(aero_config);
    delete stub;
  }

  // Test process tendencies.
  SECTION("tendencies") {
    auto* stub = new FAerosolProcessStub();
    stub->set_param("decay_rate", decay_rate);
    stub->init(aero_config);

    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    auto* progs = new Prognostics(aero_config, num_levels);
    auto* diags = new HostDiagnostics(aero_config, num_levels);
    auto* tends = new Tendencies(*progs);

    // Set initial conditions.

    // Set aerosol mass mix fractions. We start with 100% clouds
    // for simplicity.
    for (int p = 0; p < num_aero_populations; ++p) {
      for (int k = 0; k < num_levels; ++k) {
        progs->cloud_aerosols(p, k) = 1.0 / num_aero_populations;
        progs->interstitial_aerosols(p, k) = 0.0;
      }
    }

    // Set modal number concentrations.
    Real n0 = 1e6;
    for (int m = 0; m < num_modes; ++m) {
      for (int k = 0; k < num_levels; ++k) {
        progs->interstitial_num_mix_ratios(m, k) = n0;
      }
    }

    // Set gas mass mixing rations.
    for (int g = 0; g < num_gases; ++g) {
      for (int k = 0; k < num_levels; ++k) {
        progs->gases(g, k) = 1.0 / num_gases;
      }
    }

    // Now compute the tendencies by running the process.
    Real t = 0.0, dt = 0.01;
    auto team_policy = haero::TeamPolicy(1u, Kokkos::AUTO);
    auto d_stub = stub->copy_to_device();
    const auto& p = *progs;
    const auto& a = *atm;
    const auto& d = *diags;
    auto& te = *tends;
    Kokkos::parallel_for(
        team_policy, KOKKOS_LAMBDA(const TeamType& team) {
          d_stub->run(team, t, dt, p, a, d, te);
        });
    AerosolProcess::delete_on_device(d_stub);

    // --------------------------------------------------
    // Check the tendencies to make sure they make sense.
    // --------------------------------------------------

    // Cloudborne aerosol mix fractions have negative tendencies, interstitial
    // mix fractions have positive tendencies, and their sums are zero.
    auto dqdt_c = tends->cloud_aerosols;
    auto dqdt_i = tends->interstitial_aerosols;
    for (int p = 0; p < num_aero_populations; ++p) {
      for (int k = 0; k < num_levels; ++k) {
        REQUIRE(dqdt_c(p, k)[0] < 0.0);
        REQUIRE(dqdt_i(p, k)[0] > 0.0);
        Real sum = dqdt_c(p, k)[0] + dqdt_i(p, k)[0];
        REQUIRE(FloatingPoint<Real>::equiv(sum, 0.0));
      }
    }

    // Aerosol modal number concentrations are unchanged.
    auto dndt = tends->interstitial_num_mix_ratios;
    for (int m = 0; m < num_modes; ++m) {
      for (int k = 0; k < num_levels; ++k) {
        REQUIRE(FloatingPoint<Real>::equiv(dndt(m, k)[0], 0.0));
      }
    }

    // Gas mix ratios are unchanged.
    const auto& dqdt_g = tends->gases;
    for (int g = 0; g < num_gases; ++g) {
      for (int k = 0; k < num_levels; ++k) {
        REQUIRE(FloatingPoint<Real>::equiv(dqdt_g(g, k)[0], 0.0));
      }
    }

    // Clean up.
    delete progs;
    delete diags;
    delete tends;
    delete stub;
  }

  delete atm;
}
