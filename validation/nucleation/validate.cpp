#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_session.hpp"
#include "haero/processes/vehkamaki2002.hpp"
#include "skywalker/skywalker.hpp"

namespace {

using namespace haero;

// Print driver usage information and exit.
void usage(const char* exe) {
  fprintf(stderr, "%s: usage:\n", exe);
  fprintf(stderr, "%s <input.yml>\n", exe);
  exit(1);
}

// Runs an aerosol process using the parameters in param_walk, writing output
// ŧo a Python module with the given name.
void run_process(const ModalAerosolConfig& aero_config,
                 const skywalker::ParameterWalk& param_walk,
                 const char* py_module_name) {

  // Line up all the ensemble members so we can run them.
  auto inputs = param_walk.gather_inputs();
  std::vector<skywalker::OutputData> outputs;

  for (auto& input: inputs) {
    haero::PackType c_h2so4(input.user_params["c_h2so4"]);
    haero::PackType T(input.user_params["temperature"]);
    haero::PackType RH(input.user_params["relative_humidity"]);
    auto x_crit = vehkamaki2002::h2so4_critical_mole_fraction(c_h2so4, T, RH);
    auto J = vehkamaki2002::nucleation_rate(c_h2so4, T, RH, x_crit);
    skywalker::OutputData output(aero_config);
    output.metrics["nucleation_rate"] = J[0];
    outputs.push_back(output);
  }

  // Write the output data to a Python module.
  write_py_module(inputs, outputs, py_module_name);
}

} // anonymous namespace

int main(int argc, const char* argv[]) {
  ekat::initialize_ekat_session(argc, const_cast<char**>(argv), false);

  if (argc < 2) usage(argv[0]);

  // Read the input file and extract input.
  std::string input_file(argv[1]);
  try {
    ModalAerosolConfig no_config; // user-defined parameter study
    auto param_walk = skywalker::load_ensemble(no_config, input_file, "haero");

    // Set up the desired aerosol process and run it, dumping output to
    // "haero_skywalker.py".
    run_process(no_config, param_walk, "haero_nucleation.py");
  } catch (std::exception& e) {
    printf("%s: error: %s\n", argv[0], e.what());
    exit(1);
  }
  return 0;
}
