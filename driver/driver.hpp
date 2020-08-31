#ifndef HAERO_DRIVER_HPP
#define HAERO_DRIVER_HPP

#include <string>
#include <map>
#include <vector>

#include "haero/haero_config.hpp"

namespace haero {

// This type holds essential simulation input data.
struct Sim_input_data {
  // Simulation parameters.
  int do_gaschem, do_cloudchem, do_gasaerexch, do_rename, do_newnuc,
      do_coag, do_calcsize, num_unit, frac_unit, gas_unit;
  // "MET" input.
  Real temp, press, RH_CLEA, hgt, cld_frac;
  // Table of initial conditions.
  std::vector<std::map<std::string, Real> > initial_conditions;
  // Uniform perturbation factor.
  Real perturb_factor;
  // Time step and simulation duration.
  std::vector<Real> dts;
  Real duration;
  // Output parameters.
  std::string output_prefix, output_dir;
  int output_freq;
};

// This is the entry point for our C++ driver.
void haero_driver(const Sim_input_data& data);

} // namespace haero

#endif
