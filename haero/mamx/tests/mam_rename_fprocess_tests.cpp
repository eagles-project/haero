#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/processes/mam_rename_fprocess.hpp"
#include "haero/processes/mam_rename_process.hpp"

using namespace haero;

TEST_CASE("mam_rename_run", "") {
  auto aero_config = ModalAerosolConfig::create_mam4_config();
  int num_modes = aero_config.num_aerosol_modes();
  int num_aero_populations = aero_config.num_aerosol_populations;
  int num_levels = 72;

  // Set up atmospheric data
  Real pblh = 100.0;  // planetary BL height [m]
  auto atm = new Atmosphere(num_levels, pblh);

  SECTION("rename_run") {
    auto* process = new MAMRenameFProcess();
    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    auto* progs = new Prognostics(aero_config, num_levels);
    auto* diags = new HostDiagnostics(aero_config, num_levels);
    auto* tends = new Tendencies(*progs);

    // Define a pseudo-random generator [0-1) that is consistent across
    // platforms. Manually checked the first 100,000 values to be unique.
    static constexpr unsigned p0{987659};
    static constexpr unsigned p1{12373};
    long unsigned seed{54319};
    auto random = [&]() {
      seed = (seed * p1) % p0;
      return Real(seed) / p0;
    };

    // Set initial conditions
    // aerosols mass mixing ratios
    for (std::size_t p = 0; p < num_aero_populations; ++p) {
      for (std::size_t k = 0; k < num_levels; ++k) {
        progs->interstitial_aerosols(p, k) = random() * 10e-10;
        progs->cloud_aerosols(p, k) = random() * 10e-10;
      }
    }

    // aerosols number mixing ratios
    for (std::size_t imode = 0; imode < num_modes; ++imode) {
      for (std::size_t k = 0; k < num_levels; ++k) {
        progs->interstitial_num_mix_ratios(imode, k) = 1e8 + random();
        progs->cloud_num_mix_ratios(imode, k) = 1e8 + random();
      }
    }

    // Initialize the process
    process->init(aero_config);

    // Now compute the tendencies by running the process.
    Real t = 0.0, dt = 30.0;
    auto team_policy = haero::TeamPolicy(1u, Kokkos::AUTO);
    auto d_process = process->copy_to_device();
    const auto& p = *progs;
    const auto& a = *atm;
    const auto& d = *diags;
    auto& te = *tends;
    Kokkos::parallel_for(
        team_policy, KOKKOS_LAMBDA(const TeamType& team) {
          d_process->run(team, t, dt, p, a, d, te);
        });
    AerosolProcess::delete_on_device(d_process);

    delete tends;
    delete progs;
    delete diags;
    delete process;
  }

  delete atm;
}
