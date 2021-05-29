#include "haero/modal_aerosol_config.hpp"
#include "skywalker.hpp"

namespace {

// Provides a reference to a zero value.
skywalker::Real zero_value;

bool is_aerosol(const std::string& param_name) {
  return (param_name.find("aerosols.") != std::string::npos);
}

bool is_number_conc(const std::string& param_name) {
  return (param_name.find("number_conc") != std::string::npos);
}

bool is_gas(const std::string& param_name) {
  return (param_name.find("gases.") != std::string::npos);
}

bool is_atmosphere(const std::string& param_name) {
  return (param_name.find("atmosphere.") != std::string::npos);
}

void parse_aerosol(const haero::ModalAerosolConfig& aero_config,
                   const std::string& param_name, bool& cloudy,
                   int& pop_index) {
  size_t last_dot = param_name.rfind('.');
  size_t penultimate_dot = param_name.rfind('.', last_dot - 1);
  auto mode_name =
      param_name.substr(penultimate_dot + 1, last_dot - penultimate_dot - 1);
  cloudy = (param_name.find("cloudy.") != std::string::npos);
  auto aero_name = param_name.substr(last_dot + 1);
  int mode_index = aero_config.aerosol_mode_index(mode_name, false);
  int aero_index =
      aero_config.aerosol_species_index(mode_index, aero_name, false);
  pop_index = aero_config.population_index(mode_index, aero_index);
}

void parse_number_conc(const haero::ModalAerosolConfig& aero_config,
                       const std::string& param_name, bool& cloudy,
                       int& mode_index) {
  size_t last_dot = param_name.rfind('.');
  size_t penultimate_dot = param_name.rfind('.', last_dot - 1);
  auto mode_name =
      param_name.substr(penultimate_dot + 1, last_dot - penultimate_dot - 1);
  cloudy = (param_name.find("cloudy.") != std::string::npos);
  mode_index = aero_config.aerosol_mode_index(mode_name, false);
}

void parse_gas(const haero::ModalAerosolConfig& aero_config,
               const std::string& param_name, int& gas_index) {
  size_t last_dot = param_name.rfind('.');
  auto gas_name = param_name.substr(last_dot + 1);
  gas_index = aero_config.gas_index(gas_name, false);
}

}  // namespace

namespace skywalker {

Real InputData::operator[](const std::string& param_name) const {
  if (is_number_conc(param_name)) {
    bool cloud;
    int mode_index;
    parse_number_conc(aero_config, param_name, cloud, mode_index);
    if (cloud) {
      return cloud_number_concs[mode_index];
    } else {
      return interstitial_number_concs[mode_index];
    }
  } else if (is_aerosol(param_name)) {
    bool cloud;
    int pop_index;
    parse_aerosol(aero_config, param_name, cloud, pop_index);
    if (cloud) {
      return cloud_aero_mmrs[pop_index];
    } else {
      return interstitial_aero_mmrs[pop_index];
    }
  } else if (is_gas(param_name)) {
    int gas_index;
    parse_gas(aero_config, param_name, gas_index);
    return gas_mmrs[gas_index];
  } else if (is_atmosphere(param_name)) {
    if (param_name.find("temperature") != std::string::npos) {
      return temperature;
    } else if (param_name.find("pressure") != std::string::npos) {
      return pressure;
    } else if (param_name.find("relative_humidity") != std::string::npos) {
      return relative_humidity;
    } else if (param_name.find("height") != std::string::npos) {
      return height;
    } else if (param_name.find("hydrostatic_dp") != std::string::npos) {
      return hydrostatic_dp;
    } else if (param_name.find("planetary_boundary_layer_height") !=
               std::string::npos) {
      return planetary_boundary_layer_height;
    } else {
      return 0.0;
    }
  } else {
    return 0.0;
  }
}

Real& InputData::operator[](const std::string& param_name) {
  if (is_number_conc(param_name)) {
    bool cloud;
    int mode_index;
    parse_number_conc(aero_config, param_name, cloud, mode_index);
    if (cloud) {
      if (cloud_number_concs.size() <= mode_index) {
        cloud_number_concs.resize(mode_index + 1);
      }
      return cloud_number_concs[mode_index];
    } else {
      if (interstitial_number_concs.size() <= mode_index) {
        interstitial_number_concs.resize(mode_index + 1);
      }
      return interstitial_number_concs[mode_index];
    }
  } else if (is_aerosol(param_name)) {
    bool cloud;
    int pop_index;
    parse_aerosol(aero_config, param_name, cloud, pop_index);
    if (cloud) {
      if (cloud_aero_mmrs.size() <= pop_index) {
        cloud_aero_mmrs.resize(pop_index + 1);
      }
      return cloud_aero_mmrs[pop_index];
    } else {
      if (interstitial_aero_mmrs.size() <= pop_index) {
        interstitial_aero_mmrs.resize(pop_index + 1);
      }
      return interstitial_aero_mmrs[pop_index];
    }
  } else if (is_gas(param_name)) {
    int gas_index;
    parse_gas(aero_config, param_name, gas_index);
    if (gas_mmrs.size() <= gas_index) {
      gas_mmrs.resize(gas_index + 1);
    }
    return gas_mmrs[gas_index];
  } else if (is_atmosphere(param_name)) {
    if (param_name.find("temperature") != std::string::npos) {
      return temperature;
    } else if (param_name.find("pressure") != std::string::npos) {
      return pressure;
    } else if (param_name.find("relative_humidity") != std::string::npos) {
      return relative_humidity;
    } else if (param_name.find("height") != std::string::npos) {
      return height;
    } else if (param_name.find("hydrostatic_dp") != std::string::npos) {
      return hydrostatic_dp;
    } else if (param_name.find("planetary_boundary_layer_height") !=
               std::string::npos) {
      return planetary_boundary_layer_height;
    } else {
      zero_value = 0.0;
      return zero_value;
    }
  } else {
    zero_value = 0.0;
    return zero_value;
  }
}

Real OutputData::operator[](const std::string& param_name) const {
  if (is_number_conc(param_name)) {
    bool cloud;
    int mode_index;
    parse_number_conc(aero_config, param_name, cloud, mode_index);
    if (cloud) {
      return cloud_number_concs[mode_index];
    } else {
      return interstitial_number_concs[mode_index];
    }
  } else if (is_aerosol(param_name)) {
    bool cloud;
    int pop_index;
    parse_aerosol(aero_config, param_name, cloud, pop_index);
    if (cloud) {
      return cloud_aero_mmrs[pop_index];
    } else {
      return interstitial_aero_mmrs[pop_index];
    }
  } else if (is_gas(param_name)) {
    int gas_index;
    parse_gas(aero_config, param_name, gas_index);
    return gas_mmrs[gas_index];
  } else {
    return 0.0;
  }
}

std::vector<InputData> ParameterWalk::gather_inputs(const std::set<std::string>& excluded_params) const {
  // Count up the number of inputs defined by the parameter walk thingy,
  // excluding those parameters specified.
  size_t num_inputs = 1, num_params = 0;
  for (auto param: ensemble) {
    if (excluded_params.find(param.first) == excluded_params.end()) { // not excluded
      num_inputs *= param.second.size(); // set of parameter values
      num_params++;
    }
  }
  EKAT_REQUIRE_MSG(((num_params >= 1) and (num_params <= 5)),
                   "Invalid number of overridden parameters ("
                       << num_params << ", must be 1-5).");

  // Start from reference data and build a list of inputs corresponding to all
  // the overridden parameters. This involves some ugly index magic based on the
  // number of parameters.
  std::vector<InputData> inputs(num_inputs, ref_input);
  for (size_t l = 0; l < num_inputs; ++l) {
    if (num_params == 1) {
      auto iter = ensemble.begin();
      auto name = iter->first;
      const auto& vals = iter->second;
      inputs[l][name] = vals[l];
    } else if (num_params == 2) {
      auto iter = ensemble.begin();
      auto name1 = iter->first;
      const auto& vals1 = iter->second;
      iter++;
      auto name2 = iter->first;
      const auto& vals2 = iter->second;
      size_t n2 = vals2.size();
      size_t j1 = l / n2;
      size_t j2 = l - n2 * j1;
      inputs[l][name1] = vals1[j1];
      inputs[l][name2] = vals2[j2];
    } else if (num_params == 3) {
      auto iter = ensemble.begin();
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
      size_t j1 = l / (n2 * n3);
      size_t j2 = (l - n2 * n3 * j1) / n3;
      size_t j3 = l - n2 * n3 * j1 - n3 * j2;
      inputs[l][name1] = vals1[j1];
      inputs[l][name2] = vals2[j2];
      inputs[l][name3] = vals3[j3];
    } else if (num_params == 4) {
      auto iter = ensemble.begin();
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
      size_t j1 = l / (n2 * n3 * n4);
      size_t j2 = (l - n2 * n3 * n4 * j1) / (n3 * n4);
      size_t j3 = (l - n2 * n3 * n4 * j1 - n3 * n4 * j2) / n4;
      size_t j4 = l - n2 * n3 * n4 * j1 - n3 * n4 * j2 - n4 * j3;
      inputs[l][name1] = vals1[j1];
      inputs[l][name2] = vals2[j2];
      inputs[l][name3] = vals3[j3];
      inputs[l][name4] = vals4[j4];
    } else {  // if (num_params == 5)
      auto iter = ensemble.begin();
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
      size_t j1 = l / (n2 * n3 * n4 * n5);
      size_t j2 = (l - n2 * n3 * n4 * n5 * j1) / (n3 * n4 * n5);
      size_t j3 = (l - n2 * n3 * n4 * n5 * j1 - n3 * n4 * n5 * j2) / (n4 * n5);
      size_t j4 =
          (l - n2 * n3 * n4 * n5 * j1 - n3 * n4 * n5 * j2 - n4 * n5 * j3) / n5;
      size_t j5 = l - n2 * n3 * n4 * n5 * j1 - n3 * n4 * n5 * j2 -
                  n4 * n5 * j3 - n5 * j4;
      inputs[l][name1] = vals1[j1];
      inputs[l][name2] = vals2[j2];
      inputs[l][name3] = vals3[j3];
      inputs[l][name4] = vals4[j4];
      inputs[l][name5] = vals5[j5];
    }
  }

  return inputs;
}

}  // namespace skywalker

//----------------------------
// Skywalker Fortran bindings
//----------------------------

extern "C" {

typedef skywalker::Real Real;

/// Parses the given file, assuming the given named aerosol configuration,
/// returning an opaque pointer to the ensemble data.
/// @param [in] aerosol_config The named aerosol configuration. The only valid
///                            configuration at this time is "mam4".
/// @param [in] filename The name of the YAML file containing ensemble data.
void* sw_load_ensemble(const char* aerosol_config, const char* filename) {
  return nullptr;
}

/// Returns the number of inputs (members) for the given ensemble data.
int sw_ensemble_size(void* ensemble) {
  return 0;
}

/// Fetches an opaque pointer to the ith set of input data from the given
/// ensemble.
void* sw_ensemble_input(void* ensemble, int i) {
  return nullptr;
}

/// Fetches timestepping data from the given ensemble input data pointer.
void sw_input_get_timestepping(void* input, Real* dt, Real* total_time) {
}

/// Fetches atmosphere data from the given ensemble input data pointer.
void sw_input_get_atmosphere(void* input, Real* temperature, Real* pressure,
                             Real* relative_humidty, Real* height,
                             Real* hydrostatic_dp,
                             Real* planetary_boundary_layer_height) {
}

/// Fetches aerosol data from the given ensemble input data pointer.
void sw_input_get_aerosols(void* input, Real* interstitial_number_concs,
                           Real* cloud_number_concs, Real* interstitial_aero_mmrs,
                           Real* cloud_aero_mmrs) {
}

/// Fetches gas data from the given ensemble input data pointer.
void sw_input_get_gases(void* input, Real* gas_mmrs) {
}

/// Fetches an opaque pointer to the ith set of output data from the given
/// ensemble.
void* sw_ensemble_output(void* ensemble, int i) {
  return nullptr;
}

/// Sets aerosol data for the given ensemble output data pointer.
void sw_output_set_aerosols(void* output, Real* interstitial_number_concs,
                            Real* cloud_number_concs, Real* interstitial_aero_mmrs,
                            Real* cloud_aero_mmrs) {
}

/// Sets gas data for the given ensemble output data pointer.
void sw_output_set_gases(void* output, Real* gas_mmrs) {
}


/// Frees all memory associated with the ensemble, including input and output
/// data.
void sw_enѕemble_free(void* ensemble) {
}

}
