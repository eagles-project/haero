#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_session.hpp"
#include "ekat/logging/ekat_logger.hpp"
#include "haero/processes/merikanto2007.hpp"
#include "haero/processes/vehkamaki2002.hpp"
#include "skywalker/skywalker.hpp"

namespace {

using namespace haero;
using namespace ekat::logger;

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

  if (nucleation_method == 2) {  // binary nucleation
    Log::info("Validating binary nucleation (Vehkamaki et al 2002)");
  } else {
    Log::info("Validating ternary nucleation (Merikanto et al 2007)");
  }

  for (auto& input : inputs) {
    haero::PackType T(input.user_params["temperature"]);
    haero::PackType RH(input.user_params["relative_humidity"]);
    haero::PackType c_h2so4;

    haero::PackType J;             // nucleation rate
    if (nucleation_method == 2) {  // binary nucleation
      // Figure out the number concentration for H2SO4.
      auto iter = input.user_params.find("c_h2so4");
      if (iter != input.user_params.end()) {  // given as input
        c_h2so4 = haero::PackType(iter->second);
      } else {  // not given: compute nucleation threshold
        Log::info("c_h2so4 not given... computing nucleation threshold.");
        c_h2so4 = vehkamaki2002::h2so4_nucleation_threshold(T, RH);
      }
      auto x_crit = vehkamaki2002::h2so4_critical_mole_fraction(c_h2so4, T, RH);
      J = vehkamaki2002::nucleation_rate(c_h2so4, T, RH, x_crit);
    } else {  // ternary nucleation
      // Fetch the number concentration for H2SO4 gas.
      auto iter = input.user_params.find("c_h2so4");
      EKAT_REQUIRE_MSG(
          iter != input.user_params.end(),
          "c_h2so4 was not given and is required for ternary nucleation.");
      c_h2so4 = haero::PackType(iter->second);

      // Fetch the mole fraction for NH3 gas.
      haero::PackType xi_nh3;
      iter = input.user_params.find("xi_nh3");
      EKAT_REQUIRE_MSG(
          iter != input.user_params.end(),
          "xi_nh3 was not given and is required for ternary nucleation.");
      xi_nh3 = haero::PackType(iter->second);

      // Check the onset temperature, above which J = 0.
      J = 0.0;
      auto T_onset = merikanto2007::onset_temperature(RH, c_h2so4, xi_nh3);
      auto below_T_onset = (T < T_onset);
      if (below_T_onset.any()) {
        haero::PackType log_J;
        log_J.set(below_T_onset,
                  merikanto2007::log_nucleation_rate(T, RH, c_h2so4, xi_nh3));
        J.set(below_T_onset, exp(log_J));
      }
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
