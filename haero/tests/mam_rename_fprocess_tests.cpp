#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/model.hpp"
#include "haero/processes/mam_rename_fprocess.hpp"
#include "haero/processes/mam_rename_process.hpp"

using namespace haero;

Model* get_model_for_unit_tests(const ModalAerosolConfig& aero_config,
                                const std::size_t num_levels) {
  static Model* model(Model::ForUnitTests(aero_config, num_levels));
  return model;
}

TEST_CASE("mam_rename_run", "") {
  auto aero_config = create_mam4_modal_aerosol_config();
  static constexpr std::size_t num_levels{72};  // number of levels
  auto* model = get_model_for_unit_tests(aero_config, num_levels);
  const std::size_t num_gases{
      aero_config.h_gas_species.size()};  // number of gases
  const std::size_t num_modes{
      aero_config.h_aerosol_modes.size()};  // number of modes

  // Set up some prognostics aerosol data views
  const int num_aero_populations = model->num_aerosol_populations();

  Kokkos::View<PackType**> int_aerosols(
      "interstitial aerosols", num_aero_populations,
      num_levels);  // interstitial aerosols mmr [kg/kg(of air)]
  Kokkos::View<PackType**> cld_aerosols(
      "cloudborne aerosols", num_aero_populations,
      num_levels);  // cloud borne aerosols mmr [kg/kg(of air)]
  Kokkos::View<PackType**> gases("gases", num_gases, num_levels);
  Kokkos::View<PackType**> int_num_mix_ratios(
      "interstitial number mix ratios", num_modes,
      num_levels);  // interstitial aerosols number mixing ratios [#/kg(of air)]
  Kokkos::View<PackType**> cld_num_mix_ratios(
      "cloud borne number mix ratios", num_modes,
      num_levels);  // cloud borne aerosols number mixing ratios [#/kg(of air)]

  // Set up atmospheric data and populate it with some views.
  Kokkos::View<PackType*> temp("temperature", num_levels);  //[K]
  Kokkos::View<PackType*> press("pressure", num_levels);    //[Pa]
  Kokkos::View<PackType*> rel_hum("relative humidity", num_levels);
  Kokkos::View<PackType*> pdel("hydrostatic_dp", num_levels);  //[Pa]
  Kokkos::View<PackType*> ht("height", num_levels + 1);        //[m]
  Real pblh{100.0};  // planetary BL height [m]
  auto atm =
      std::make_unique<Atmosphere>(num_levels, temp, press, rel_hum, ht, pdel,
                                   pblh);  // create atmosphere object

  SECTION("rename_run") {
    auto process = std::make_unique<MAMRenameFProcess>();
    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    auto* progs = model->create_prognostics(int_aerosols, cld_aerosols,
                                            int_num_mix_ratios, cld_num_mix_ratios,
                                            gases);
    auto* diags = model->create_diagnostics();
    auto tends = std::make_unique<Tendencies>(*progs);

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
        int_aerosols(p, k) = random() * 10e-10;
        cld_aerosols(p, k) = random() * 10e-10;
      }
    }

    // aerosols number mixing ratios
    for (std::size_t imode = 0; imode < num_modes; ++imode) {
      for (std::size_t k = 0; k < num_levels; ++k) {
        int_num_mix_ratios(imode, k) = 1e8 + random();
        cld_num_mix_ratios(imode, k) = 1e8 + random();
      }
    }

    // Initialize the process
    process->init(aero_config);

    // Now compute the tendencies by running the process.
    Real t = 0.0, dt = 30.0;
    process->run(aero_config, t, dt, *progs, *atm, *diags, *tends);
  }
}
