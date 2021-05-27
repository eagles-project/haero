#include "skywalker.hpp"

namespace {

// Used to return a reference for a parameter with an invalid name.
haero::Real zero_value = 0;

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

haero::Real InputData::operator[](const std::string& param_name) const {
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

haero::Real& InputData::operator[](const std::string& param_name) {
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

haero::Real OutputData::operator[](const std::string& param_name) const {
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

}  // namespace skywalker
