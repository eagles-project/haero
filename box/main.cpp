// This is a box model used to validate Haero against MAM's box model. It
// behaves identically to that box model.

#include "box_input.hpp"
#include "box_output.hpp"
#include <haero/aerosol_process.hpp>
#include <haero/atmosphere.hpp>
#include <haero/diagnostics.hpp>
#include <haero/haero.hpp>
#include <haero/modal_aerosol_config.hpp>
#include <haero/prognostics.hpp>

#include <ekat/ekat_session.hpp>

#include <cstdio>
#include <cstdlib>
#include <limits>

using namespace std;
using namespace haero;

namespace {

// Print version banner.
void print_banner() {
  string version = haero::version();
  string revision = haero::revision();
  bool not_for_science = haero::has_uncommitted_changes();
  printf("haero v%s (git revision %s)\n", version.c_str(), revision.c_str());
  if (not_for_science) {
    fprintf(stderr, "WARNING: haero was built with uncommitted changeѕ!\n");
    fprintf(stderr, "WARNING: DO NOT USE FOR SCIENCE!\n");
  }
}

// Print driver usage information and exit.
void usage(const char* exe) {
  fprintf(stderr, "%s: usage:\n", exe);
  fprintf(stderr, "%s <input.nml>\n", exe);
  exit(1);
}

// Create the processes selected within the given input.
std::vector<AerosolProcess*> create_processes(const ModalAerosolConfig& config,
                                              const BoxInput& input) {
  std::vector<AerosolProcess*> processes;
  return processes;
}

// Create prognostic variables from input.
Prognostics* create_prognostics(const ModalAerosolConfig& config,
                                const BoxInput& input) {
  return nullptr;
}

// Create atmospheric variables from input.
Atmosphere* create_atmosphere(const BoxInput& input) {
  return nullptr;
}

// Create diagnostic variables from input.
HostDiagnostics* create_diagnostics(const ModalAerosolConfig& config,
                                    const BoxInput& input) {
  return nullptr;
}

// Advance the selected processes using the same sequential splitting approach
// used in MAM4.
void advance(const BoxInput& input,
             const std::vector<AerosolProcess*>& processes,
             Prognostics& prognostics,
             Atmosphere& atmosphere,
             HostDiagnostics& diagnostics) {
}

}  // anonymous namespace

int main(int argc, const char** argv) {
  ekat::initialize_ekat_session(argc, const_cast<char**>(argv), false);
  print_banner();

  if (argc < 2) usage(argv[0]);

  // Read the input file and extract input.
  std::string input_file(argv[1]);
  auto input = Box::read_namelist(input_file);

  // Initialize a set of processes with the given input.
  auto config = ModalAerosolConfig::create_mam4_config();
  std::vector<AerosolProcess*> processes = create_processes(config, input);

  // Initialize aerosol prognostics and atmospheric conditions with input.
  Prognostics* prognostics = create_prognostics(config, input);
  Atmosphere* atmosphere = create_atmosphere(input);
  HostDiagnostics* diagnostics = create_diagnostics(config, input);

  // Open the output file and define quantities.
  Box::BoxOutput output(config);

  for (int n = 0; n < input.mam_nstep; ++n) {
    advance(input, processes, *prognostics, *atmosphere, *diagnostics);
    if ((n % input.mam_output_intvl) == 0) {
      output.append(*prognostics, *diagnostics);
    }
  }

  // Write output.
  output.write("haero_output.nc");

  // Clean up.
  delete diagnostics;
  delete atmosphere;
  delete prognostics;
  for (size_t i = 0; i < processes.size(); ++i) {
    delete processes[i];
  }
  ekat::finalize_ekat_session();
  return 0;
}
