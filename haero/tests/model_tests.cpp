#include <cmath>
#include <iostream>
#include <memory>

#include "catch2/catch.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "haero/model.hpp"

using namespace haero;

TEST_CASE("model_ctor", "") {
  // define aerosol modes
  const auto modes = create_mam4_modes();
  std::cout << "created mam4 modes: ";
  for (int m = 0; m < modes.size(); ++m) {
    std::cout << modes[m].name() << " ";
  }
  std::cout << "\n";
  // define aerosol species
  const auto aero_species = create_mam4_aerosol_species();
  std::cout << "created mam4 aerosol species: ";
  for (int s = 0; s < aero_species.size(); ++s) {
    std::cout << aero_species[s].symbol() << " ";
  }
  std::cout << "\n";
  // map aerosols to modes
  const auto mode_spec_map = create_mam4_mode_species();
  std::cout << "mode_spec_map.size() = " << mode_spec_map.size() << "\n";
  // define gas species
  const auto gas_species = create_mam4_gas_species();
  std::cout << "created mam4 gas species: ";
  for (int g = 0; g < gas_species.size(); ++g) {
    std::cout << gas_species[g].symbol() << " ";
  }
  std::cout << "\n";
  // configure the aerosol model
  const auto aero_config =
      ModalAerosolConfig(modes, aero_species, mode_spec_map, gas_species);

  {
    // All of the above steps can be combined:
    const auto short_config = create_mam4_modal_aerosol_config();
    std::cout << short_config.info_string();
  }

  // select processes (for now, just nucleation)
  SelectedProcesses aerosol_processes;
  aerosol_processes.nucleation = SelectedProcesses::Nucleation::MAMNucleation;
  {
    // later:
    // const auto aerosol_processes = create_mam4_aerosol_processes();
  }

  // build the model
  const int nlev = 72;
  Model haero_model(aero_config, aerosol_processes, nlev);
  const int num_vert_packs = PackInfo::num_packs(nlev);

  /** Host model jobs start here */
  // create views of tracer data (species/mode mass mixing ratios)
  auto aero_tracers = SpeciesColumnView(
      "aerosol_tracers", haero_model.num_aerosol_populations(), num_vert_packs);
  auto mode_tracers = ModeColumnView("aerosol_mode_tracers",
                                     haero_model.num_modes(), num_vert_packs);
  auto gas_tracers =
      SpeciesColumnView("gas_tracers", haero_model.num_gases(), num_vert_packs);
  /** end host model jobs */

  auto progs = std::unique_ptr<Prognostics>(haero_model.create_prognostics(
      aero_tracers, aero_tracers, gas_tracers, mode_tracers, mode_tracers));
}
