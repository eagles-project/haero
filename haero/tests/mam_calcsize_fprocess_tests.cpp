#include <yaml-cpp/yaml.h>

#include <memory>

#include "catch2/catch.hpp"
#include "haero/processes/mam_calcsize_fprocess.hpp"
#include "haero/processes/mam_calcsize_process.hpp"

using namespace haero;

// ADD COMMENTS
TEST_CASE("mam_calcsize_run", "") {
  // Ensures that Fortran implementation is ALWAYS run with a packsize of 1
  static_assert(HAERO_PACK_SIZE == 1,
                "Fortran not supported for HAERO_PACK_SIZE != 1.");

  /*-----------------------------------------------------------------------------
    We construct a phony configuration and initialize only the parts which we
   need to drive the test. Some of the fields defined below are created to
   satisfy the arguments needed to create the model and call the "run" method
   -----------------------------------------------------------------------------*/

  auto aero_config =
      ModalAerosolConfig::create_mam4_config();  // create MAM4 configuration

  static constexpr int num_levels{72};  // number of levels

  const size_t num_gases = aero_config.gas_species.size();    // number of gases
  const size_t num_modes = aero_config.aerosol_modes.size();  // number of modes
  const int num_aero_populations =
      aero_config.num_aerosol_populations;  // total number of aerosol species

  // Set up atmospheric data
  Real pblh = 100.0;  // planetary BL height [m]
  auto atm = new Atmosphere(num_levels, pblh);

  // This will drive the "run" method of calcsize
  SECTION("calcsize_run") {
    auto process = new MAMCalcsizeFProcess();

    // Initialize the process
    process->init(aero_config);

    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    auto* progs = new Prognostics(aero_config, num_levels);
    auto* diags = new HostDiagnostics(aero_config, num_levels);
    auto* tends = new Tendencies(*progs);

    // open and read calcsize data from a YAML file
    std::string datafile = HAERO_TEST_DATA_DIR;
    datafile += "/calcsize_input.yaml";

    auto calcsize_data = YAML::LoadFile(datafile);
    for (YAML::const_iterator it = calcsize_data.begin();
         it != calcsize_data.end(); ++it) {
      // read input collection
      const std::string& key = it->first.as<std::string>();

      //--read all its attributes
      auto attributes = it->second;
      //----read contents of each attribute
      auto intermmr = attributes["interstitial"];
      auto intermmr_num = attributes["interstitial_num"];
      auto cldbrnmmr = attributes["cldbrn"];
      auto cldbrnmmr_num = attributes["cldbrn_num"];

      // mmrs
      for (int p = 0; p < num_aero_populations; ++p) {
        for (int k = 0; k < num_levels; ++k) {
          progs->interstitial_aerosols(p, k) = intermmr[p].as<float>();
          progs->cloud_aerosols(p, k) = cldbrnmmr[p].as<float>();
        }
      }

      // number mmrs
      for (int imode = 0; imode < num_modes; ++imode) {
        for (int k = 0; k < num_levels; ++k) {
          progs->interstitial_num_mix_ratios(imode, k) =
              intermmr_num[imode].as<float>();
          progs->cloud_num_mix_ratios(imode, k) =
              cldbrnmmr_num[imode].as<float>();
        }
      }

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
    }

    // Clean up.
    delete tends;
    delete atm;
    delete progs;
    delete diags;
    delete process;
  }  // section:calcsize_run

}  // TEST_CASE mam_calcsize_run
