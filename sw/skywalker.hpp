#ifndef HAERO_SKYWALKER_HPP
#define HAERO_SKYWALKER_HPP

#include "haero/haero.hpp"
#include "haero/modal_aerosol_config.hpp"

namespace skywalker {

// Reference input data for simulations.
struct InputData {
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

  // Override a parameter value in the given InputData given its name.
  void override_parameter(const haero::ModalAerosolConfig& aero_config,
                          const std::string& param_name,
                          haero::Real param_value);
};


// Data structure that stores parameter walking information.
struct ParameterWalk {
  // aerosol configuration
  haero::ModalAerosolConfig aero_config;

  // process name
  std::string process;

  // reference input data
  InputData ref_input;

  // ensemble: parameters to walk (name -> vector of values)
  std::map<std::string, std::vector<haero::Real>> ensemble;
};

// Here's a container that associates input parameters with output data.
struct OutputData {
  // Modal aerosol number concentrations [# aero molecules / kg air]
  std::vector<haero::Real> number_concs;
  // Aerosol mass mixing ratios [kg aerosol / kg air]
  std::vector<haero::Real> interstitial_aero_mmrs, cloud_aero_mmrs;
  // Gas mass mixing ratios [kg gas / kg air]
  std::vector<haero::Real> gas_mmrs;
};

}

#endif
