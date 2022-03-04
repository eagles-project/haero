#include <ekat/ekat.hpp>
#include <ekat/logging/ekat_logger.hpp>
#include <memory>

#include "catch2/catch.hpp"
#include "haero/processes/mam_calcsize_hostcxx_process.hpp"

using namespace haero;

// ADD COMMENTS
TEST_CASE("mam_calcsize_hostcxx_run", "") {
  /*-----------------------------------------------------------------------------
    We construct a phony model and initialize only the parts which we need to
    drive the test. Some of the fields defined below are created to satisfy the
    arguments needed to create the model and call the "run" method
   -----------------------------------------------------------------------------*/

  auto aero_config =
      ModalAerosolConfig::create_mam4_config();  // create MAM4 configuration

  int num_gases = aero_config.num_gases();
  int num_modes = aero_config.num_aerosol_modes();
  int num_aero_populations = aero_config.num_aerosol_populations;
  int num_levels = 72;
  int num_vert_packs = PackInfo::num_packs(num_levels);
  int num_iface_packs = PackInfo::num_packs(num_levels + 1);

  // Set up atmospheric data.
  Real pblh = 100.0;  // planetary BL height [m]
  auto atm = new Atmosphere(num_levels, pblh);

  // This will drive the "run" method of calcsize_hostcxx
  SECTION("calcsize_hostcxx_run") {
    auto process = new MAMCalcsizeHostCXXProcess();

    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    auto progs = new Prognostics(aero_config, num_levels);
    auto diags = new HostDiagnostics(aero_config, num_levels);
    auto tends = new Tendencies(*progs);

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
    for (int p = 0; p < num_aero_populations; ++p) {
      for (int k = 0; k < num_vert_packs; ++k) {
        progs->interstitial_aerosols(p, k) = random() * 10e-10;
        progs->cloud_aerosols(p, k) = random() * 10e-10;
      }
    }

    // aerosols number mixing ratios
    for (int imode = 0; imode < num_modes; ++imode) {
      for (int k = 0; k < num_vert_packs; ++k) {
        progs->interstitial_num_mix_ratios(imode, k) = 1e8 + random();
        progs->cloud_num_mix_ratios(imode, k) = 1e8 + random();
      }
    }

    ekat::logger::Log::set_level(ekat::logger::Log::level::debug);

    // Initialize the process
    process->init(aero_config);

    // Now compute the tendencies by running the process.
    Real t = 0.0, dt = 30.0;
    auto team_policy = haero::TeamPolicy(1u, Kokkos::AUTO);
    const auto& p = *progs;
    const auto& a = *atm;
    const auto& d = *diags;
    auto& te = *tends;
    Kokkos::parallel_for(
        team_policy, KOKKOS_LAMBDA(const TeamType& team) {
          process->run(team, t, dt, p, a, d, te);
        });

    // Clean up.
    delete tends;
    delete atm;
    delete progs;
    delete diags;
    delete process;
  }  // section:calcsize_hostcxx_run

}  // TEST_CASE mam_calcsize_hostcxx_run
