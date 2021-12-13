#include <haero/processes/kerminen2002.hpp>

#include <validation.hpp>
#include <skywalker.hpp>

#include <iostream>

// This driver evaluates the particle growth parameterization for the nucleation
// process, using input generated from the MAM box model.

void usage(const std::string& prog_name) {
  std::cerr << prog_name << ": usage:" << std::endl;
  std::cerr << prog_name << " <input.yaml>" << std::endl;
  exit(0);
}

using namespace skywalker;
using namespace haero;

int main(int argc, char** argv) {
  if (argc == 1) {
    usage((const char*)argv[0]);
  }
  std::string input_file = argv[1];
  std::string output_file = validation::output_name(input_file);
  std::cout << argv[0] << ": reading " << input_file << std::endl;

  // Load the ensemble. Any error encountered is fatal.
  Ensemble* ensemble = skywalker::load_ensemble(input_file, "haero");

  // We don't need any settings for this particular test.
  // Settings settings = ensemble->settings();

  // Run the ensemble.
  try {
    ensemble->process([](const Input& input, Output& output) {
      // Ensemble parameters
      PackType c_so4 = input.get("c_so4");
      PackType c_nh4 = input.get("c_nh4");
      PackType nh4_to_so4_molar_ratio = input.get("nh4_to_so4_molar_ratio");
      PackType temp = input.get("temperature");
      PackType rel_hum = input.get("relative_humidity");
      PackType d_dry_crit = input.get("dry_critical_diameter");
      PackType d_wet_crit = input.get("wet_critical_diameter");
      PackType d_dry_grown = input.get("dry_grown_diameter");
      PackType rho_grown = input.get("grown_mass_density");
      PackType rho_air = input.get("air_mass_density");
      Real mw_h2so4 = input.get("mw_h2so4");
      if (input.has("h2so4_uptake_rate")) {
      }

      // Compute the growth rate GR, corrected for NH3/H2O uptake
      const auto bounded_rel_hum = max(0.10, min(0.95, rel_hum));
      const auto wet_dry_vol_ratio = 1.0 - 0.56 / log(bounded_rel_hum);
      PackType V_frac_wet_so4 =
          1.0 / (wet_dry_vol_ratio *
          (1.0 + nh4_to_so4_molar_ratio * 17.0 / 98.0));
      PackType GR = kerminen2002::growth_rate(c_so4, rho_grown, mw_h2so4, temp);
      GR /= V_frac_wet_so4;

      // Compute the condensation sink CS'.
      PackType d_wet_grown = 1e9 * d_dry_grown *
          pow(wet_dry_vol_ratio, 1.0 / 3.0);
      PackType CS1 = kerminen2002::condensation_sink(rho_air, d_wet_grown,
                                                     c_so4 + c_nh4);

      // Compute the growth parameter eta.
      PackType eta = kerminen2002::growth_parameter(c_so4, c_nh4,
          nh4_to_so4_molar_ratio, temp, rel_hum, d_dry_crit, d_wet_crit,
          d_dry_grown, rho_grown, rho_air, mw_h2so4);

      // Jot everything down
      output.set("GR", GR[0]);
      output.set("CS1", CS1[0]);
      output.set("eta", eta[0]);
    });

  } catch (Exception& e) {
    std::cerr << argv[0] << ": Error: " << e.what() << std::endl;
  }

  // Write out a Python module.
  std::cout << argv[0] << ": writing " << output_file << std::endl;
  ensemble->write(output_file);

  // Clean up.
  delete ensemble;
}
