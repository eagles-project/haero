#include "skywalker.hpp"

namespace skywalker {

void InputData::override_parameter(const haero::ModalAerosolConfig& aero_config,
                                   const std::string& param_name,
                                   haero::Real param_value) {
  size_t colon = param_name.find(':');
  if (colon != std::string::npos) { // mode:aerosol
    auto mode_name = param_name.substr(0, colon);
    auto aero_name = param_name.substr(colon+1, param_name.length());
    int mode_index = aero_config.aerosol_mode_index(mode_name);
    int aero_index = aero_config.aerosol_species_index(mode_index, aero_name, false);
    int pop_index = aero_config.population_index(mode_index, aero_index);
    interstitial_aero_mmrs[pop_index] = param_value;
    cloud_aero_mmrs[pop_index] = param_value;
  } else {
    // gas or mode?
    int mode_index = aero_config.aerosol_mode_index(param_name);
    int gas_index = aero_config.gas_index(param_name);
    if (mode_index != -1) {
      number_concs[mode_index] = param_value;
    } else if (gas_index != -1) {
      gas_mmrs[gas_index] = param_value;
    } else {
      // Atmospheric state variables?
      if (param_name == "temperature") {
        temperature = param_value;
      } else if (param_name == "pressure") {
        pressure = param_value;
      } else if (param_name == "relative_humidity") {
        relative_humidity = param_value;
      } else if (param_name == "height") {
        height = param_value;
      } else if (param_name == "hydroѕtatic_dp") {
        hydrostatic_dp = param_value;
      } else if (param_name == "planetary_boundary_layer_height") {
        planetary_boundary_layer_height = param_value;
      } else { // other
      }
    }
  }
}

}
