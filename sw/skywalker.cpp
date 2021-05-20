#include "skywalker.hpp"

namespace {

bool is_aerosol(const std::string& param_name) {
  return (param_name.find("aerosol:") != std::string::npos);
}

bool is_number_conc(const std::string& param_name) {
  return (param_name.find("number_conc") != std::string::npos);
}

bool is_gas(const std::string& param_name) {
  return (param_name.find("gases:") != std::string::npos);
}

void parse_aerosol(const haero::ModalAerosolConfig& aero_config,
  const std::string& param_name, bool& cloud, int& pop_index) {
  size_t last_colon = param_name.rfind(':');
  size_t penultimate_colon = param_name.rfind(':', last_colon-1);
  auto mode_name = param_name.substr(penultimate_colon+1, last_colon-penultimate_colon-1);
  cloud = (param_name.find("cloud:") != std::string::npos);
  auto aero_name = param_name.substr(last_colon+1);
  int mode_index = aero_config.aerosol_mode_index(mode_name);
  int aero_index = aero_config.aerosol_species_index(mode_index, aero_name, false);
  pop_index = aero_config.population_index(mode_index, aero_index);
}

void parse_number_conc(const haero::ModalAerosolConfig& aero_config,
  const std::string& param_name, bool& cloud, int& mode_index) {
  size_t last_colon = param_name.rfind(':');
  size_t penultimate_colon = param_name.rfind(':', last_colon-1);
  auto mode_name = param_name.substr(penultimate_colon+1, last_colon-penultimate_colon-1);
  cloud = (param_name.find("cloud:") != std::string::npos);
  mode_index = aero_config.aerosol_mode_index(param_name);
}

void parse_gas(const haero::ModalAerosolConfig& aero_config,
  const std::string& param_name, int& gas_index) {
  size_t last_colon = param_name.rfind(':');
  auto gas_name = param_name.substr(last_colon+1);
  gas_index = aero_config.gas_index(gas_name, false);
}

}

namespace skywalker {

haero::Real InputData::operator[](const std::string& param_name) const {
  if (is_aerosol(param_name)) {
    bool cloud;
    int pop_index;
    parse_aerosol(aero_config, param_name, cloud, pop_index);
    if (cloud) {
      return cloud_aero_mmrs[pop_index];
    } else {
      return interstitial_aero_mmrs[pop_index];
    }
  } else if (is_number_conc(param_name)) {
    bool cloud;
    int mode_index;
    parse_number_conc(aero_config, param_name, cloud, mode_index);
    if (cloud) {
      return cloud_number_concs[mode_index];
    } else {
      return interstitial_number_concs[mode_index];
    }
  } else if (is_gas(param_name)) {
    int gas_index;
    parse_gas(aero_config, param_name, gas_index);
    return gas_mmrs[gas_index];
  } else {
    if (param_name == "temperature") {
    return temperature;
    } else if (param_name == "pressure") {
    return pressure;
    } else if (param_name == "relative_humidity") {
    return relative_humidity;
    } else if (param_name == "height") {
      return height;
    } else if (param_name == "hydroѕtatic_dp") {
      return hydrostatic_dp;
    } else if (param_name == "planetary_boundary_layer_height") {
      return planetary_boundary_layer_height;
    } else {
      return 0.0;
    }
  }
}

haero::Real& InputData::operator[](const std::string& param_name) {
  if (is_aerosol(param_name)) {
    bool cloud;
    int pop_index;
    parse_aerosol(aero_config, param_name, cloud, pop_index);
    if (cloud) {
      return cloud_aero_mmrs[pop_index];
    } else {
      return interstitial_aero_mmrs[pop_index];
    }
  } else if (is_number_conc(param_name)) {
    bool cloud;
    int mode_index;
    parse_number_conc(aero_config, param_name, cloud, mode_index);
    if (cloud) {
      return cloud_number_concs[mode_index];
    } else {
      return interstitial_number_concs[mode_index];
    }
  } else if (is_gas(param_name)) {
    int gas_index;
    parse_gas(aero_config, param_name, gas_index);
    return gas_mmrs[gas_index];
  } else {
    if (param_name == "temperature") {
    return temperature;
    } else if (param_name == "pressure") {
    return pressure;
    } else if (param_name == "relative_humidity") {
    return relative_humidity;
    } else if (param_name == "height") {
      return height;
    } else if (param_name == "hydroѕtatic_dp") {
      return hydrostatic_dp;
    } else if (param_name == "planetary_boundary_layer_height") {
      return planetary_boundary_layer_height;
    } else {
      return temperature; // Not great...
    }
  }
}

haero::Real OutputData::operator[](const std::string& param_name) const {
  if (is_aerosol(param_name)) {
    bool cloud;
    int pop_index;
    parse_aerosol(aero_config, param_name, cloud, pop_index);
    if (cloud) {
      return cloud_aero_mmrs[pop_index];
    } else {
      return interstitial_aero_mmrs[pop_index];
    }
  } else if (is_number_conc(param_name)) {
    bool cloud;
    int mode_index;
    parse_number_conc(aero_config, param_name, cloud, mode_index);
    if (cloud) {
      return cloud_number_concs[mode_index];
    } else {
      return interstitial_number_concs[mode_index];
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
