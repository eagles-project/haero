#include <haero/processes/vehkamaki2002.hpp>
#include <iostream>
#include <skywalker.hpp>
#include <validation.hpp>

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
  std::string output_file = validation::output_name(input_file);
  std::cout << argv[0] << ": reading " << input_file << std::endl;

  // Load the ensemble. Any error encountered is fatal.
  Ensemble* ensemble = skywalker::load_ensemble(input_file, "haero");

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

      // Compute the threshold of H2SO4 above which nucleation occurs.
      PackType c_thresh =
          vehkamaki2002::h2so4_nucleation_threshold(temp, rel_hum);

      // Write the computed nucleation rate.
      output.set("nucleation_threshold", c_thresh[0]);
    });
  } catch (Exception& e) {
    std::cerr << argv[0] << ": Error: " << e.what() << std::endl;
  }

  // Write out a Python module.
  std::cout << argv[0] << ": writing "
            << oset(input skywalker_gasaerexch_uptkrates_1box1gas.yaml)
                   utput_file
            << std::endl;
  ensemble->write(output_file);

  // Clean up.
  delete ensemble;
}
