#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/model.hpp"
#include "haero/processes/mam_rename_process.hpp"

using namespace haero;

Model* get_model_for_unit_tests(const ModalAerosolConfig& aero_config,
                                const std::size_t num_levels) {
  static Model* model(Model::ForUnitTests(aero_config, num_levels));
  return model;
}

TEST_CASE("mam_rename_run", "") {
  using View1D = Kokkos::View<PackType*>;

  auto aero_config = create_mam4_modal_aerosol_config();
  static constexpr std::size_t num_levels{72};  // number of levels
  std::size_t num_vert_packs = num_levels / HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels) {
    num_vert_packs++;
  }
  std::size_t num_iface_packs = (num_levels + 1) / HAERO_PACK_SIZE;
  if (num_iface_packs * HAERO_PACK_SIZE < (num_levels + 1)) {
    num_iface_packs++;
  }

  auto* model = get_model_for_unit_tests(aero_config, num_levels);
  const std::size_t num_gases{
      aero_config.gas_species.size()};  // number of gases
  const std::size_t num_modes{
      aero_config.aerosol_modes.size()};  // number of modes

  // Set up some prognostics aerosol data views
  const int num_aero_populations = model->num_aerosol_populations();

  SpeciesColumnView int_aerosols(
      "interstitial aerosols", num_aero_populations,
      num_vert_packs);  // interstitial aerosols mmr [kg/kg(of air)]
  SpeciesColumnView cld_aerosols(
      "cloudborne aerosols", num_aero_populations,
      num_vert_packs);  // cloud borne aerosols mmr [kg/kg(of air)]
  SpeciesColumnView gases("gases", num_gases, num_vert_packs);
  ModeColumnView int_num_concs("interstitial number concs", num_modes,
                               num_vert_packs);  // interstitial aerosols number
                                                 // mixing ratios [#/kg(of air)]
  ModeColumnView cld_num_concs("cloud borne number concs", num_modes,
                               num_vert_packs);  // cloud borne aerosols number
                                                 // mixing ratios [#/kg(of air)]

  // Set up atmospheric data and populate it with some views.
  View1D temp("temperature", num_vert_packs);  //[K]
  View1D press("pressure", num_vert_packs);    //[Pa]
  View1D rel_hum("relative humidity", num_vert_packs);
  View1D qv("vapor mixing ratio", num_vert_packs);
  View1D pdel("hydrostatic_dp", num_vert_packs);  //[Pa]
  View1D ht("height", num_iface_packs);
  Real pblh{100.0};  // planetary BL height [m]
  auto atm = new Atmosphere(num_levels, temp, press, qv, ht, pdel, pblh);

  SECTION("rename_run") {
    auto process = new MAMRenameProcess();
    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    auto* progs = model->create_prognostics(
        int_aerosols, cld_aerosols, int_num_concs, cld_num_concs, gases);
    auto* diags = model->create_diagnostics();
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
    auto h_int_num_concs = Kokkos::create_mirror_view(int_num_concs);
    auto h_cld_num_concs = Kokkos::create_mirror_view(cld_num_concs);
    for (std::size_t imode = 0; imode < num_modes; ++imode) {
      for (std::size_t k = 0; k < num_vert_packs; ++k) {
        h_int_num_concs(imode, k) = 1e8 + random();
        h_cld_num_concs(imode, k) = 1e8 + random();
      }
    }
    Kokkos::deep_copy(int_num_concs, h_int_num_concs);
    Kokkos::deep_copy(cld_num_concs, h_cld_num_concs);

    // Initialize the process
    process->init(aero_config);

    // Now compute the tendencies by running the process.
    Real t = 0.0, dt = 30.0;
    auto team_policy = haero::TeamPolicy(1u, Kokkos::AUTO);
    auto d_process = process->copy_to_device();
    Kokkos::parallel_for(
        team_policy, KOKKOS_LAMBDA(const TeamType& team) {
          d_process->run(team, t, dt, *progs, *atm, *diags, *tends);
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
