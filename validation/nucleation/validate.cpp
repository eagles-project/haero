#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_session.hpp"
#include "haero/processes/vehkamaki2002.hpp"
#include "skywalker/skywalker.hpp"

namespace {

using namespace haero;

// Print driver usage information and exit.
void usage(const char* exe) {
  fprintf(stderr, "%s: usage:\n", exe);
  fprintf(stderr, "%s <input.yaml>\n", exe);
  exit(1);
}

// Runs an aerosol process using the parameters in param_walk, writing output
// ŧo a Python module with the given name.
void run_process(const ModalAerosolConfig& aero_config,
                 const skywalker::ParameterWalk& param_walk,
                 const std::string& py_module_name) {
  // Binary or ternary nucleation?
  auto program_params = param_walk.program_params;
  int nucleation_method = atoi(program_params["nucleation_method"].c_str());
  if ((nucleation_method != 2) and (nucleation_method != 3)) {
    fprintf(stderr, "Invalid nucleation method: %d (must be 2 or 3)\n",
            nucleation_method);
    return;
  }

  // Line up all the ensemble members so we can run them.
  auto inputs = param_walk.gather_inputs();
  std::vector<skywalker::OutputData> outputs;

  for (auto& input : inputs) {
    haero::PackType c_h2so4(input.user_params["c_h2so4"]);
    haero::PackType T(input.user_params["temperature"]);
    haero::PackType RH(input.user_params["relative_humidity"]);

    haero::PackType J;  // nucleation rate
    if (nucleation_method == 2) {
      auto x_crit = vehkamaki2002::h2so4_critical_mole_fraction(c_h2so4, T, RH);
      J = vehkamaki2002::nucleation_rate(c_h2so4, T, RH, x_crit);
    } else {
      // TODO
    }
    skywalker::OutputData output(aero_config);
    output.metrics["nucleation_rate"] = J[0];
    outputs.push_back(output);
  }

  // Write the output data to a Python module.
  skywalker::write_py_module(inputs, outputs, py_module_name);
}

}  // anonymous namespace

int main(int argc, const char* argv[]) {
  ekat::initialize_ekat_session(argc, const_cast<char**>(argv), false);

  if (argc < 2) usage(argv[0]);

  // Compute input and output filenames.
  std::string input_file(argv[1]);
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

  // Read the input file and extract input.
  try {
    ModalAerosolConfig no_config;  // user-defined parameter study
    auto param_walk = skywalker::load_ensemble(no_config, input_file, "haero");

    // Set up the desired aerosol process and run it, dumping output to
    // "haero_skywalker.py".
    run_process(no_config, param_walk, output_file);
  } catch (std::exception& e) {
    printf("%s: error: %s\n", argv[0], e.what());
    exit(1);
  }
  return 0;
}
