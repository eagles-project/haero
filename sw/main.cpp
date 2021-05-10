#include "parse_yaml.hpp"
#include "haero/model.hpp"
#include "haero/available_processes.hpp"

#include "ekat/ekat_session.hpp"
#include "ekat/ekat_pack_kokkos.hpp"

namespace {

// Print driver usage information and exit.
void usage(const char* exe)
{
  fprintf(stderr, "%s: usage:\n", exe);
  fprintf(stderr, "%s <input.yml>\n", exe);
  exit(1);
}

// List of overridden parameter (name, value)s.
using OverriddenParameterList = std::vector<std::pair<std::string, haero::Real>>;

// Here's a container that associates input parameters with output data.
struct OutputData {
  // Names and values of overridden input parameters.
  OverriddenParameterList params;
  // Modal aerosol number concentrations [# aero molecules / kg air]
  std::vector<haero::Real> number_concs;
  // Aerosol mass mixing ratios [kg aerosol / kg air]
  std::vector<std::vector<haero::Real>> aero_mmrs;
  // Gas mass mixing ratios [kg gas / kg air]
  std::vector<haero::Real> gas_mmrs;
};

// Writes simulation output data to a Python module.
void write_py_module(const char* py_module_name,
                     const std::vector<OutputData>& data) {

}

// Overrides the parameter with the given name using the given value within
// the aerosol prognostics or atmospheric state.
void override_parameter(const haero::ModalAerosolConfig& aero_config,
                        haero::Prognostics& prognostics,
                        haero::Atmosphere& atmosphere,
                        int level,
                        const std::string& param_name,
                        haero::Real param_value) {
  size_t colon = param_name.find(':');
  if (colon != std::string::npos) { // mode:aerosol
    auto mode_name = param_name.substr(0, colon);
    auto aero_name = param_name.substr(colon+1, param_name.length());
    int mode_index = aero_config.aerosol_mode_index(mode_name);
    int aero_index = aero_config.aerosol_species_index(mode_index, aero_name, false);
    int pop_index = aero_config.population_index(mode_index, aero_index);
    auto int_aerosols = ekat::scalarize(prognostics.interstitial_aerosols);
    int_aerosols(pop_index, level) = param_value;
  } else {
    // gas or mode?
    int mode_index = aero_config.aerosol_mode_index(param_name);
    int gas_index = aero_config.gas_index(param_name);
    if (mode_index != -1) {
      auto num_concs = ekat::scalarize(prognostics.interstitial_num_concs);
      num_concs(mode_index, level) = param_value;
    } else if (gas_index != -1) {
      auto gases = ekat::scalarize(prognostics.gases);
      gases(gas_index, level) = param_value;
    } else {
      // Atmospheric state variables?
      if (param_name == "temperature") {
        auto T = ekat::scalarize(atmosphere.temperature);
        T(level) = param_value;
      } else if (param_name == "pressure") {
        auto p = ekat::scalarize(atmosphere.pressure);
        p(level) = param_value;
      } else if (param_name == "relative_humidity") {
        auto relhum = ekat::scalarize(atmosphere.relative_humidity);
        relhum(level) = param_value;
      } else if (param_name == "height") {
        auto h = ekat::scalarize(atmosphere.height);
        h(level) = param_value;
      } else if (param_name == "hydroѕtatic_dp") {
        auto dp = ekat::scalarize(atmosphere.hydrostatic_dp);
        dp(level) = param_value;
      } else { // other
      }
    }
  }
}

// Initializes prognostic and atmosphere input data according to the
// (non-plbh) parameters given in the walker. Returns an exhaustive list of
// ŧhe sets of overridden input parameters used to initialize the containers.
std::vector<OverriddenParameterList>
initialize_input(const haero::ModalAerosolConfig& aero_config,
                 const ParameterWalk& param_walk,
                 haero::Atmosphere& atmosphere,
                 haero::Prognostics& prognostics) {
  // Copy in reference data.
  int num_levels = prognostics.num_levels();
  int num_modes = prognostics.num_aerosol_modes();
  int num_gases = prognostics.num_gases();
  {
    auto T = ekat::scalarize(atmosphere.temperature);
    auto p = ekat::scalarize(atmosphere.pressure);
    auto relhum = ekat::scalarize(atmosphere.relative_humidity);
    auto h = ekat::scalarize(atmosphere.height);
    auto dp = ekat::scalarize(atmosphere.hydrostatic_dp);
    auto int_aero = ekat::scalarize(prognostics.interstitial_aerosols);
    auto cld_aero = ekat::scalarize(prognostics.cloud_aerosols);
    auto gases = ekat::scalarize(prognostics.gases);
    auto num_concs = ekat::scalarize(prognostics.interstitial_num_concs);
    for (int l = 0; l < num_levels; ++l) {
      // Atmospheric state
      T(l) = param_walk.temperature;
      p(l) = param_walk.pressure;
      relhum(l) = param_walk.relative_humidity;
      h(l) = param_walk.height;
      dp(l) = param_walk.hydrostatic_dp;

      // Aerosol prognostics.
      for (int m = 0; m < num_modes; ++m) {
        num_concs(m, l) = param_walk.number_concs[m];
        for (size_t s = 0; s < param_walk.aero_mmrs[m].size(); ++s) {
          int p = aero_config.population_index(m, static_cast<int>(s));
          int_aero(p, l) = param_walk.aero_mmrs[m][s];
        }
        // TODO: Cloud aerosols not yet treated!
      }
      for (int g = 0; g < num_gases; ++g) {
        gases(g, l) = param_walk.gas_mmrs[g];
      }
    }
  }

  // How many non-plbh parameters are we overriding?
  int num_params = param_walk.ensemble.size();
  if (param_walk.ensemble.find("planetary_boundary_layer_height") != param_walk.ensemble.end()) {
    num_params--;
  }
  EKAT_REQUIRE_MSG(((num_params >= 1) and (num_params <= 5)),
    "Invalid number of overridden parameters (" << num_params << ", must be 1-5).");

  // Override the requested parameters at each level.
  // This involves some ugly index magic based on the number of parameters.
  std::vector<OverriddenParameterList> overridden_params;
  for (int l = 0; l < num_levels; ++l) {
    if (num_params == 1) {
      auto iter = param_walk.ensemble.begin();
      auto name = iter->first;
      const auto& vals = iter->second;
      override_parameter(aero_config, prognostics, atmosphere, l, name, vals[l]);
      overridden_params.push_back({{name, vals[l]}});
    } else if (num_params == 2) {
      auto iter = param_walk.ensemble.begin();
      auto name1 = iter->first;
      const auto& vals1 = iter->second;
      iter++;
      auto name2 = iter->first;
      const auto& vals2 = iter->second;
      size_t n2 = vals2.size();
      size_t j1 = l/n2;
      size_t j2 = l - n2*j1;
      override_parameter(aero_config, prognostics, atmosphere, l, name1, vals1[j1]);
      override_parameter(aero_config, prognostics, atmosphere, l, name2, vals2[j2]);
      overridden_params.push_back({{name1, vals1[j1]},
                                   {name2, vals2[j2]}});
    } else if (num_params == 3) {
      auto iter = param_walk.ensemble.begin();
      auto name1 = iter->first;
      const auto& vals1 = iter->second;
      iter++;
      auto name2 = iter->first;
      const auto& vals2 = iter->second;
      iter++;
      auto name3 = iter->first;
      const auto& vals3 = iter->second;
      size_t n2 = vals2.size();
      size_t n3 = vals3.size();
      size_t j1 = l/(n2*n3);
      size_t j2 = (l - n2*n3*j1) / n3;
      size_t j3 = l - n2*n3*j1 - n3*j2;
      override_parameter(aero_config, prognostics, atmosphere, l, name1, vals1[j1]);
      override_parameter(aero_config, prognostics, atmosphere, l, name2, vals2[j2]);
      override_parameter(aero_config, prognostics, atmosphere, l, name3, vals3[j3]);
      overridden_params.push_back({{name1, vals1[j1]},
                                   {name2, vals2[j2]},
                                   {name3, vals3[j3]}});
    } else if (num_params == 4) {
      auto iter = param_walk.ensemble.begin();
      auto name1 = iter->first;
      const auto& vals1 = iter->second;
      iter++;
      auto name2 = iter->first;
      const auto& vals2 = iter->second;
      iter++;
      auto name3 = iter->first;
      const auto& vals3 = iter->second;
      iter++;
      auto name4 = iter->first;
      const auto& vals4 = iter->second;
      size_t n2 = vals2.size();
      size_t n3 = vals3.size();
      size_t n4 = vals4.size();
      size_t j1 = l/(n2*n3*n4);
      size_t j2 = (l - n2*n3*n4*j1) / (n3*n4);
      size_t j3 = (l - n2*n3*n4*j1 - n3*n4*j2) / n4;
      size_t j4 = l - n2*n3*n4*j1 - n3*n4*j2 - n4*j3;
      override_parameter(aero_config, prognostics, atmosphere, l, name1, vals1[j1]);
      override_parameter(aero_config, prognostics, atmosphere, l, name2, vals2[j2]);
      override_parameter(aero_config, prognostics, atmosphere, l, name3, vals3[j3]);
      override_parameter(aero_config, prognostics, atmosphere, l, name4, vals4[j4]);
      overridden_params.push_back({{name1, vals1[j1]},
                                   {name2, vals2[j2]},
                                   {name3, vals3[j3]},
                                   {name4, vals4[j4]}});
    } else { // if (num_params == 5)
      auto iter = param_walk.ensemble.begin();
      auto name1 = iter->first;
      const auto& vals1 = iter->second;
      iter++;
      auto name2 = iter->first;
      const auto& vals2 = iter->second;
      iter++;
      auto name3 = iter->first;
      const auto& vals3 = iter->second;
      iter++;
      auto name4 = iter->first;
      const auto& vals4 = iter->second;
      iter++;
      auto name5 = iter->first;
      const auto& vals5 = iter->second;
      size_t n2 = vals2.size();
      size_t n3 = vals3.size();
      size_t n4 = vals4.size();
      size_t n5 = vals5.size();
      size_t j1 = l/(n2*n3*n4*n5);
      size_t j2 = (l - n2*n3*n4*n5*j1) / (n3*n4*n5);
      size_t j3 = (l - n2*n3*n4*n5*j1 - n3*n4*n5*j2) / (n4*n5);
      size_t j4 = (l - n2*n3*n4*n5*j1 - n3*n4*n5*j2 - n4*n5*j3) / n5;
      size_t j5 = l - n2*n3*n4*n5*j1 - n3*n4*n5*j2 - n4*n5*j3 - n5*j4;
      override_parameter(aero_config, prognostics, atmosphere, l, name1, vals1[j1]);
      override_parameter(aero_config, prognostics, atmosphere, l, name2, vals2[j2]);
      override_parameter(aero_config, prognostics, atmosphere, l, name3, vals3[j3]);
      override_parameter(aero_config, prognostics, atmosphere, l, name4, vals4[j4]);
      override_parameter(aero_config, prognostics, atmosphere, l, name5, vals5[j5]);
      overridden_params.push_back({{name1, vals1[j1]},
                                   {name2, vals2[j2]},
                                   {name3, vals3[j3]},
                                   {name4, vals4[j4]},
                                   {name5, vals5[j5]}});
    }
  }

  return overridden_params;
}

// Runs an aerosol process using the parameters in param_walk, writing output
// ŧo a Python module with the given name.
void run_process(const haero::ModalAerosolConfig& aero_config,
                 const ParameterWalk& param_walk,
                 const char* py_module_name) {
  // Count up the number of simulations we need (excluding the planetary
  // boundary layer parameter). We can run all simulations simultaneously
  // by setting data for each simulation at a specific vertical level.
  std::vector<haero::Real> pblhs, dts;
  int num_levels = 1;
  for (auto iter = param_walk.ensemble.begin();
       iter != param_walk.ensemble.end(); ++iter) {
    if (iter->first == "planetary_boundary_layer_height") {
      pblhs = iter->second;
    } else if (iter->first == "dt") {
      dts = iter->second;
    } else {
      num_levels *= static_cast<int>(iter->second.size());
    }
  }
  if (pblhs.empty()) { // pblh is not a walked parameter!
    pblhs.push_back(param_walk.planetary_boundary_layer_height);
  }
  if (dts.empty()) { // dt is not a walked parameter!
    dts.push_back(param_walk.dt);
  }

  // Create a model initialized for a number of vertical levels equal to the
  // number of (0D) simulations we need for our parameter walk.
  haero::Model* model = haero::Model::ForUnitTests(aero_config, num_levels);

  // Create column views containing reference data.
  int num_modes = aero_config.h_aerosol_modes.size();
  int num_aero_populations = model->num_aerosol_populations();
  int num_gases = aero_config.h_gas_species.size();
  haero::SpeciesColumnView int_aerosols("interstitial aerosols",
                                        num_aero_populations, num_levels);
  haero::SpeciesColumnView cld_aerosols("cloud aerosols",
                                        num_aero_populations, num_levels);
  haero::SpeciesColumnView gases("gases", num_gases, num_levels);
  haero::ModalColumnView int_num_concs("interstitial number concs", num_modes,
                                       num_levels);
  haero::ModalColumnView cld_num_concs(" cloudborne number concs", num_modes,
                                       num_levels);

  auto* prognostics = model->create_prognostics(int_aerosols, cld_aerosols,
                                                gases, int_num_concs, cld_num_concs);
  auto* diagnostics = model->create_diagnostics();

  // Set up an atmospheric state and initialize it with reference data.
  haero::ColumnView temp("temperature", num_levels);
  haero::ColumnView press("pressure", num_levels);
  haero::ColumnView rel_hum("relative humidity", num_levels);
  haero::ColumnView ht("height", num_levels+1);
  haero::ColumnView dp("hydrostatic pressure thickness", num_levels);
  auto* atmosphere = new haero::Atmosphere(num_levels, temp, press, rel_hum, ht,
                       dp, param_walk.planetary_boundary_layer_height);

  // Create tendencies for the given prognostics.
  auto* tendencies = new haero::Tendencies(*prognostics);

  // Create the specified process.
  haero::AerosolProcess* process = nullptr;
  if (param_walk.process == "MAMNucleationProcess") { // C++ nucleation
    process = new haero::MAMNucleationProcess();
#if HAERO_FORTRAN
  } else if (param_walk.process == "MAMNucleationFProcess") { // fortran nucleation
    process = new haero::MAMNucleationFProcess();
#endif
  } else { // unknown
    fprintf(stderr, "Unknown aerosol process: %s\n", param_walk.process.c_str());
    return;
  }

  // Initialize it for the given aerosol configuration.
  process->init(aero_config);

  // Run a set of simulations, each of which computes tendencies for aerosols
  // and gases over different vertical levels, to maximize parallelism. We do
  // include two loops here to accommodate different values of the planetary
  // boundary layer height and the use of different time steps.
  std::vector<OutputData> output_data;
  for (auto pblh: pblhs) {
    // Initialize the input.
    auto overridden_params = initialize_input(aero_config, param_walk,
      *atmosphere, *prognostics);
    atmosphere->set_planetary_boundary_height(pblh);

    for (auto dt: dts) {
      // Run the thing.
      haero::Real t = 0.0, t_end = param_walk.total_time;
      while (t < t_end) {
        process->run(aero_config, t, dt, *prognostics, *atmosphere, *diagnostics,
                     *tendencies);

        // Advance the time and prognostic state.
        t += dt;
        prognostics->scale_and_add(dt, *tendencies);
      }

      // If the planetary boundary layer height is actually a walked parameter,
      // add it each one of our list of overridden parameters.
      if (pblhs.size() > 1) {
        for (int l = 0; l < num_levels; ++l) {
          overridden_params[l].push_back({"planetary_boundary_layer_height", pblh});
        }
      }

      // Same for time steps.
      if (dts.size() > 1) {
        for (int l = 0; l < num_levels; ++l) {
          overridden_params[l].push_back({"dt", dt});
        }
      }

      // Stash output data.
      for (auto params: overridden_params) {
        OutputData output;
        output.params = params;
        output_data.push_back(output);
      }
    }
  }

  // Write the output data to a Python module.
  write_py_module(py_module_name, output_data);

  // Clean up.
  delete tendencies;
  delete atmosphere;
  delete diagnostics;
  delete prognostics;
  delete process;
}

} // anonymous namespace

int main(int argc, const char* argv[]) {
  ekat::initialize_ekat_session(argc, const_cast<char**>(argv), false);

  if (argc < 2)
    usage(argv[0]);

  // Set up an aerosol configuration. For now, we just use MAM4.
  auto aero_config = haero::create_mam4_modal_aerosol_config();

  // Read the input file and extract input.
  std::string input_file(argv[1]);
  try {
    auto param_walk = parse_yaml(aero_config, input_file);

    // Set up the desired aerosol process and run it, dumping output to
    // "haero_skywalker.py".
    run_process(aero_config, param_walk, "haero_skywalker");
  } catch (std::exception& e) {
    printf("%s: error: %s\n", argv[0], e.what());
    exit(1);
  }
}
