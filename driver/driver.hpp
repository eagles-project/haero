#ifndef HAERO_DRIVER_HPP
#define HAERO_DRIVER_HPP

#include <map>
#include <string>
#include <vector>

#include "haero/aerosol_species.hpp"
#include "haero/gas_species.hpp"
#include "haero/haero_config.hpp"
#include "haero/mode.hpp"
#include "haero/physical_constants.hpp"
#include "host_params.hpp"

namespace haero {
namespace driver {

/// This type contains switches for activating/deactivating aerosol-related
/// physical processes.
struct AerosolProcesses {
  /// Model growth of aerosol particles?
  bool growth;
  /// Model gas chemistry?
  bool gas_chemistry;
  /// Model cloud chemistry?
  bool cloud_chemistry;
  /// Model gas-aerosol exchange processes?
  bool gas_aerosol_exchange;
  /// Allow the merging of modes?
  bool mode_merging;
  /// Model aerosol nucleation/formation processes?
  bool nucleation;
  /// Model aerosol coagulation?
  bool coagulation;
};

/// This type defines parameters related to atmospheric columns.
struct GridParams {
  /// Number of columns.
  int num_columns;
  /// Number of vertical levels per column.
  int num_levels;
};

/// This holds initial conditions for aerosol/gas species. Below, modes are
/// indexed in the same order as given in SimulationInput::modes (below).
struct InitialConditions {
  /// Aerosol species column profiles, stored from the top of the atmosphere
  /// downward. Specified in modal mass fractions for each mode in
  /// which they appear. So aerosols[mode_index][aerosol_name][k] gives the
  /// mass fraction of the given aerosol within the given mode at elevation k.
  /// All mass fractions for a given mode must sum to 1.
  /// (aerosols[mode index][aerosol species symbol] -> mass fraction).
  std::vector<std::map<std::string, std::vector<Real> > > aerosols;
  /// Gas species. Specified in mole fractions [kmol gas / kmol air].
  /// (gases[gas species symbol][k] -> mole fraction of the gas at elevation k).
  std::map<std::string, std::vector<Real> > gases;
  /// Modal number densities [#/m**3].
  /// (modes[mode_index][k] -> modal number density at elevation k.
  std::vector<std::vector<Real> > modes;
};

/// This is a simplified surrogate for aerosol chemistry. It'll go away when
/// we adopt a real solution.
struct SimpleChemistryModel {
  /// Uniform production rates for gases
  /// (gas symbol -> molar production rate [mol gas / mol air / s]
  std::map<std::string, Real> production_rates;
};

/// This type defines time-stepping-related parameters for a simulation.
struct SimulationParams {
  /// Fixed time step size [s]
  Real dt;
  /// Simulation duration [s]
  Real duration;
  /// Prefix for NetCDF output file.
  std::string output_prefix;
  /// Directory to which output file is written (absolute path).
  std::string output_dir;
  /// Output frequency in steps, or -1 if undefined.
  int output_freq;
};

/// This type holds essential simulation input data for the haero stand-alone
/// driver.
struct SimulationInput {
  /// List of aerosol modes.
  std::vector<Mode> modes;
  /// List of aerosol species.
  std::vector<AerosolSpecies> aerosols;
  /// List of gas species.
  std::vector<GasSpecies> gases;
  /// Settings for modeling specific processes.
  AerosolProcesses physics;
  /// Ambient atmosphere conditions.
  AtmosphericConditions atmosphere;
  /// Column grid parameters.
  GridParams grid;
  /// Initial conditions.
  InitialConditions initial_conditions;
  /// Chemistry modeling.
  SimpleChemistryModel chemistry;
  /// Simulation parameters (timestepping, duration, output)
  SimulationParams simulation;
};

/// This is the entry point for the haero driver, which executes the simulations
/// in the given ensemble.
void haero_driver(const std::vector<SimulationInput>& ensemble);

}  // namespace driver
}  // namespace haero

#endif
