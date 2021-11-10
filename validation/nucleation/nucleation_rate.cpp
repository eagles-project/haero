#include <haero/merikanto2007.hpp>
#include <haero/vehkamaki2002.hpp>
#include <haero/wang2008.hpp>

#include <skywalker.hpp>

#include <iostream>

// This driver computes the binary or ternary nucleation rate for the given
// input.

void usage(const std::string& prog_name) {
  std::cerr << prog_name << ": usage:" << std::endl;
  std::cerr << prog_name << " <input.yaml>" << std::endl;
  exit(0);
}

using namespace skywalker;

void run_vehkamaki2002(Ensemble* ensemble, int pbl_method) {
  // Ensemble data
  ensemble->process([](const Input& input, Output& output) {
    // Ensemble parameters
    Real c_h2so4 = input.get("c_h2so4");
    Real rel_hum = input.get("relative_humidity");
    Real temp = input.get("temperature");

    Real J = 0.0;

    // Write the computed nucleation rate.
    output.set("nucleation_rate", J);
  });
}

void run_merikanto2007(Ensemble* ensemble, int pbl_method) {
}

int main(int argc, char **argv) {

  if (argc == 1) {
    usage((const char*)argv[0]);
  }
  std::string input_file = argv[1];
  std::string output_file;
  {
    size_t slash = input_file.find_last_of('/');
    size_t dot = input_file.find_last_of('.');
    if ((dot == std::string::npos) and (slash == std::string::npos)) {
      dot = input_file.length();
    }
    if (slash == std::string::npos) {
      slash = 0;
    } else {
      slash += 1;
      dot -= slash;
    }
    output_file = std::string("haero_") + input_file.substr(slash, dot) +
                  std::string(".py");
  }

  // Load the ensemble. Any error encountered is fatal.
  Ensemble* ensemble = skywalker::load_ensemble(input_file, "settings");

  // Figure out settings for binary/ternary nucleation and planetary boundary
  // layer treatment
  Settings settings = ensemble->settings();
  int nuc_method = std::stoi(settings.get("nucleation_method"));
  if ((nuc_method != 2) and (nuc_method != 3)) {
    std::cerr << "Invalid nucleation method: " << nuc_method << std::endl;
    exit(0);
  }

  int pbl_method = std::stoi(settings.get("pbl_method"));
  if ((pbl_method < 0) or (pbl_method > 2)) {
    std::cerr << "Invalid planetary boundary layer method: " << pbl_method
              << std::endl;
    exit(0);
  }

  // Run the ensemble.
  if (nuc_method == 2) { // binary nucleation
    run_vehkamaki2002(ensemble, pbl_method);
  } else { // ternary nucleation
    run_merikanto2007(ensemble, pbl_method);
  }

  // Write out a Python module.
  ensemble->write(output_file);

  // Clean up.
  delete ensemble;
}
