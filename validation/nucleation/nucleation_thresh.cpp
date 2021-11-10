#include <haero/processes/vehkamaki2002.hpp>
#include <iostream>
#include <skywalker.hpp>

// This driver computes the threshold concentration of H2SO4 for binary
// nucleation.

void usage(const std::string& prog_name) {
  std::cerr << prog_name << ": usage:" << std::endl;
  std::cerr << prog_name << " <input.yaml>" << std::endl;
  exit(0);
}

using namespace skywalker;
using namespace haero;

void run(Ensemble* ensemble, int pbl_method) {}

int main(int argc, char** argv) {
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
  std::cout << argv[0] << ": reading " << input_file << std::endl;

  // Load the ensemble. Any error encountered is fatal.
  Ensemble* ensemble = skywalker::load_ensemble(input_file, "settings");

  Settings settings = ensemble->settings();
  if (settings.has("nucleation_method")) {
    int nuc_method = std::stoi(settings.get("nucleation_method"));
    if ((nuc_method != 2)) {
      std::cerr << "Invalid nucleation method: " << nuc_method << std::endl;
      exit(0);
    }
  }

  // Run the ensemble.
  try {
    ensemble->process([](const Input& input, Output& output) {
      // Ensemble parameters
      PackType rel_hum = input.get("relative_humidity");
      PackType temp = input.get("temperature");

      // Compute the mole fraction of H2SO4 in a critical cluster, and from it
      // the nucleation rate.
      PackType c_thresh =
          vehkamaki2002::h2so4_nucleation_threshold(temp, rel_hum);

      // Write the computed nucleation rate.
      output.set("nucleation_threshold", c_thresh[0]);
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
