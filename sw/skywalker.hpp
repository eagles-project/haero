#ifndef HAERO_SKYWALKER_HPP
#define HAERO_SKYWALKER_HPP

#include "haero/haero.hpp"
#include "haero/modal_aerosol_config.hpp"

namespace skywalker {

// Reference input data for simulations.
struct InputData {

  explicit InputData(const haero::ModalAerosolConfig& aero_config):
    aero_config(aero_config) {}
  InputData() = delete;

  const haero::ModalAerosolConfig& aero_config;

  // timestepping parameters
  haero::Real dt, total_time;

  // atmospheric state parameters
  haero::Real temperature, pressure, relative_humidity, height,
              hydrostatic_dp, planetary_boundary_layer_height;

  // aerosol initial data

  // Modal aerosol number concentrations [# aero molecules / kg air]
  std::vector<haero::Real> number_concs;
  // Aerosol mass mixing ratios [kg aerosol / kg air]
  std::vector<haero::Real> interstitial_aero_mmrs, cloud_aero_mmrs;
  // Gas mass mixing ratios [kg gas / kg air]
  std::vector<haero::Real> gas_mmrs;

  // Fetches the parameter with the given name.
  haero::Real operator[](const std::string& param_name) const;
  haero::Real& operator[](const std::string& param_name);
};


// Data structure that stores parameter walking information.
struct ParameterWalk {
  explicit ParameterWalk(const haero::ModalAerosolConfig& aero_config):
    aero_config(aero_config), process(), ref_input(aero_config) {}
  ParameterWalk() = delete;

  // aerosol configuration
  const haero::ModalAerosolConfig& aero_config;

  // process name
  std::string process;

  // reference input data
  InputData ref_input;

  // ensemble: parameters to walk (name -> vector of values)
  std::map<std::string, std::vector<haero::Real>> ensemble;
};

// Here's a container that associates input parameters with output data.
struct OutputData {
  explicit OutputData(const haero::ModalAerosolConfig& aero_config):
    aero_config(aero_config) {}
  OutputData() = delete;

  const haero::ModalAerosolConfig& aero_config;

  // Modal aerosol number concentrations [# aero molecules / kg air]
  std::vector<haero::Real> number_concs;
  // Aerosol mass mixing ratios [kg aerosol / kg air]
  std::vector<haero::Real> interstitial_aero_mmrs, cloud_aero_mmrs;
  // Gas mass mixing ratios [kg gas / kg air]
  std::vector<haero::Real> gas_mmrs;

  // Fetches the parameter with the given name.
  haero::Real operator[](const std::string& param_name) const;

};

}

#endif
