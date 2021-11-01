#include <ekat/ekat.hpp>
#include <ekat/logging/ekat_logger.hpp>
#include <memory>

#include "catch2/catch.hpp"
#include "haero/model.hpp"
#include "haero/processes/mam_calcsize_hostcxx_process.hpp"

using namespace haero;

Model* get_model_for_unit_tests(const ModalAerosolConfig& aero_config,
                                const int num_levels) {
  static Model* model(Model::ForUnitTests(aero_config, num_levels));
  return model;
}

// ADD COMMENTS
TEST_CASE("mam_calcsize_hostcxx_run", "") {
  /*-----------------------------------------------------------------------------
    We construct a phony model and initialize only the parts which we need to
    drive the test. Some of the fields defined below are created to satisfy the
    arguments needed to create the model and call the "run" method
   -----------------------------------------------------------------------------*/

  auto aero_config =
      create_mam4_modal_aerosol_config();  // create MAM4 configuration

  static constexpr int num_levels{72};  // number of levels

  std::size_t num_vert_packs = num_levels / HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels) {
    num_vert_packs++;
  }
  std::size_t num_iface_packs = (num_levels + 1) / HAERO_PACK_SIZE;
  if (num_iface_packs * HAERO_PACK_SIZE < (num_levels + 1)) {
    num_iface_packs++;
  }

  auto* model = get_model_for_unit_tests(
      aero_config, num_levels);  // get an instance of "model"

  const size_t num_gases = aero_config.gas_species.size();    // number of gases
  const size_t num_modes = aero_config.aerosol_modes.size();  // number of modes

  // Set up some prognostics aerosol data views
  const int num_aero_populations{
      model->num_aerosol_populations()};  // total number of aerosol species
  Kokkos::View<PackType**> int_aerosols(
      "interstitial aerosols", num_aero_populations,
      num_vert_packs);  // interstitial aerosols mmr [kg/kg(of air)]
  Kokkos::View<PackType**> cld_aerosols(
      "cloudborne aerosols", num_aero_populations,
      num_vert_packs);  // cloud borne aerosols mmr [kg/kg(of air)]
  Kokkos::View<PackType**> gases("gases", num_gases, num_vert_packs);
  Kokkos::View<PackType**> int_num_mix_ratios(
      "interstitial number mix ratios", num_modes,
      num_vert_packs);  // interstitial aerosols number mixing ratios [#/kg(of
                        // air)]
  Kokkos::View<PackType**> cld_num_mix_ratios(
      "cloud borne number mix ratios", num_modes,
      num_vert_packs);  // cloud borne aerosols number mixing ratios [#/kg(of
                        // air)]

  // Set up atmospheric data and populate it with some views.
  Kokkos::View<PackType*> temp("temperature", num_vert_packs);  //[K]
  Kokkos::View<PackType*> press("pressure", num_vert_packs);    //[Pa]
  Kokkos::View<PackType*> qv("vapor mixing ratio", num_vert_packs);
  Kokkos::View<PackType*> pdel("hydrostatic_dp", num_vert_packs);  //[Pa]
  Kokkos::View<PackType*> ht("height", num_iface_packs);           //[m]
  Real pblh{100.0};  // planetary BL height [m]
  auto atm = new Atmosphere(num_levels, temp, press, qv, ht, pdel, pblh);

  // This will drive the "run" method of calcsize_hostcxx
  SECTION("calcsize_hostcxx_run") {
    auto process = new MAMCalcsizeHostCXXProcess();

    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    auto progs = model->create_prognostics(int_aerosols, cld_aerosols,
                                           int_num_mix_ratios,
                                           cld_num_mix_ratios, gases);
    auto diags = model->create_diagnostics();
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
        int_aerosols(p, k) = random() * 10e-10;
        cld_aerosols(p, k) = random() * 10e-10;
      }
    }

    // aerosols number mixing ratios
    for (int imode = 0; imode < num_modes; ++imode) {
      for (int k = 0; k < num_vert_packs; ++k) {
        int_num_mix_ratios(imode, k) = 1e8 + random();
        cld_num_mix_ratios(imode, k) = 1e8 + random();
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