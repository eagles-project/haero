#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/processes/mam_rename_process.hpp"

using namespace haero;

TEST_CASE("mam_rename_run", "") {
  using View1D = Kokkos::View<PackType*>;

  auto aero_config = ModalAerosolConfig::create_mam4_config();
  std::size_t num_levels = 72;  // number of levels
  std::size_t num_vert_packs = PackInfo::num_packs(num_levels);

  const std::size_t num_modes =
      aero_config.aerosol_modes.size();  // number of modes

  // Set up some prognostics aerosol data views
  const int num_aero_populations = aero_config.num_aerosol_populations;

  Real pblh = 100.0;  // planetary BL height [m]
  auto atm = new Atmosphere(num_levels, pblh);

  SECTION("rename_run") {
    auto process = new MAMRenameProcess();
    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    auto* progs = new Prognostics(aero_config, num_levels);
    auto* diags = new HostDiagnostics(aero_config, num_levels);
    auto* tends = new Tendencies(*progs);

    SpeciesColumnView int_aerosols = progs->interstitial_aerosols;
    SpeciesColumnView cld_aerosols = progs->cloud_aerosols;
    SpeciesColumnView gases = progs->gases;
    ModeColumnView int_num_ratios = progs->interstitial_num_mix_ratios;
    ModeColumnView cld_num_ratios = progs->cloud_num_mix_ratios;

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
    auto h_int_aerosols = Kokkos::create_mirror_view(int_aerosols);
    auto h_cld_aerosols = Kokkos::create_mirror_view(cld_aerosols);
    for (std::size_t p = 0; p < num_aero_populations; ++p) {
      for (std::size_t k = 0; k < num_vert_packs; ++k) {
        h_int_aerosols(p, k) = random() * 10e-10;
        h_cld_aerosols(p, k) = random() * 10e-10;
      }
    }
    Kokkos::deep_copy(int_aerosols, h_int_aerosols);
    Kokkos::deep_copy(cld_aerosols, h_cld_aerosols);

    // aerosols number mixing ratios
    auto h_int_num_ratios = Kokkos::create_mirror_view(int_num_ratios);
    auto h_cld_num_ratios = Kokkos::create_mirror_view(cld_num_ratios);
    for (std::size_t imode = 0; imode < num_modes; ++imode) {
      for (std::size_t k = 0; k < num_vert_packs; ++k) {
        h_int_num_ratios(imode, k) = 1e8 + random();
        h_cld_num_ratios(imode, k) = 1e8 + random();
      }
    }
    Kokkos::deep_copy(int_num_ratios, h_int_num_ratios);
    Kokkos::deep_copy(cld_num_ratios, h_cld_num_ratios);

    // Initialize the process
    process->init(aero_config);

    // Now compute the tendencies by running the process.
    Real t = 0.0, dt = 30.0;
    auto team_policy = haero::TeamPolicy(1u, Kokkos::AUTO);
    auto d_process = process->copy_to_device();
    const auto& p = *progs;
    const auto& a = *atm;
    const Diagnostics& d = *diags;
    auto& te = *tends;
    Kokkos::parallel_for(
        team_policy, KOKKOS_LAMBDA(const TeamType& team) {
          d_process->run(team, t, dt, p, a, d, te);
        });
    AerosolProcess::delete_on_device(d_process);

    delete progs;
    delete diags;
    delete tends;
    delete process;
  }

  // Clean up.
  delete atm;
}
