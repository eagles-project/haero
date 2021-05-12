#include "skywalker.hpp"

namespace skywalker {

haero::Real InputData::operator[](const std::string& param_name) const {
  size_t colon = param_name.find(':');
  if (colon != std::string::npos) {  // mode:aerosol
    auto mode_name = param_name.substr(0, colon);
    auto aero_name = param_name.substr(colon + 1, param_name.length());
    int mode_index = aero_config.aerosol_mode_index(mode_name);
    int aero_index =
        aero_config.aerosol_species_index(mode_index, aero_name, false);
    int pop_index = aero_config.population_index(mode_index, aero_index);
    return interstitial_aero_mmrs[pop_index];  // TODO: what about cloud
                                               // aerosols?
  } else {
    // gas or mode?
    int mode_index = aero_config.aerosol_mode_index(param_name);
    int gas_index = aero_config.gas_index(param_name);
    if (mode_index != -1) {
      return number_concs[mode_index];
    } else if (gas_index != -1) {
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
}

haero::Real& InputData::operator[](const std::string& param_name) {
  size_t colon = param_name.find(':');
  if (colon != std::string::npos) {  // mode:aerosol
    auto mode_name = param_name.substr(0, colon);
    auto aero_name = param_name.substr(colon + 1, param_name.length());
    int mode_index = aero_config.aerosol_mode_index(mode_name);
    int aero_index =
        aero_config.aerosol_species_index(mode_index, aero_name, false);
    int pop_index = aero_config.population_index(mode_index, aero_index);
    return interstitial_aero_mmrs[pop_index];  // TODO: what about cloud
                                               // aerosols?
  } else {
    // gas or mode?
    int mode_index = aero_config.aerosol_mode_index(param_name);
    int gas_index = aero_config.gas_index(param_name);
    if (mode_index != -1) {
      return number_concs[mode_index];
    } else if (gas_index != -1) {
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
        return planetary_boundary_layer_height;  // TODO
      }
    }
  }
}

haero::Real OutputData::operator[](const std::string& param_name) const {
  size_t colon = param_name.find(':');
  if (colon != std::string::npos) {  // mode:aerosol
    auto mode_name = param_name.substr(0, colon);
    auto aero_name = param_name.substr(colon + 1, param_name.length());
    int mode_index = aero_config.aerosol_mode_index(mode_name);
    int aero_index =
        aero_config.aerosol_species_index(mode_index, aero_name, false);
    int pop_index = aero_config.population_index(mode_index, aero_index);
    return interstitial_aero_mmrs[pop_index];  // TODO: what about cloud
                                               // aerosols?
  } else {
    // gas or mode?
    int mode_index = aero_config.aerosol_mode_index(param_name);
    int gas_index = aero_config.gas_index(param_name);
    if (mode_index != -1) {
      return number_concs[mode_index];
    } else if (gas_index != -1) {
      return gas_mmrs[gas_index];
    } else {
      return 0.0;
    }
  }
}

}  // namespace skywalker
