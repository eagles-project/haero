// This is Haero's standalone C++ driver. It reads a YAML file containing input
// parameters for a simulation (or an ensemble of simulations) and dispatches
// the runs to a Fortran driver subroutine.

#include <cstdio>
#include <cstdlib>
#include <limits>
#include "haero/haero.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/ekat_session.hpp"
#include "driver/driver.hpp"

using namespace std;
using namespace haero;

namespace {

// Print version banner.
void print_banner()
{
  string version = haero::version();
  string revision = haero::revision();
  bool not_for_science = haero::has_uncommitted_changes();
  printf("haero v%s (git revision %s)\n", version.c_str(), revision.c_str());
  if (not_for_science)
  {
    fprintf(stderr, "WARNING: haero was built with uncommitted changeѕ!\n");
    fprintf(stderr, "WARNING: DO NOT USE FOR SCIENCE!\n");
  }
}

// Print driver usage information and exit.
void usage(const char* exe)
{
  fprintf(stderr, "%s: usage:\n", exe);
  fprintf(stderr, "%s <input.yml>\n", exe);
  exit(1);
}

// Extract simulation input from a parsed YAML file.
std::vector<SimulationInput> extract_input(const ekat::ParameterList& yaml)
{
  std::vector<SimulationInput> ensemble;
  return ensemble;
}

} // anonymous namespace

int main(int argc, const char** argv)
{
  ekat::initialize_ekat_session(argc, const_cast<char**>(argv));
  print_banner();

  if (argc < 2)
    usage(argv[0]);

  // Read the input file.
  std::string input_file(argv[1]);
  auto yaml = ekat::parse_yaml_file(input_file);

  // Extract input. This produces one or more simulations belonging to an
  // ensemble, depending on the initial conditions specified in the yaml file.
  std::vector<SimulationInput> ensemble = extract_input(yaml);

  // Now run the driver with our input.
  haero_driver(ensemble);

  // Clean up.
  ekat::finalize_ekat_session();
  return 0;
}
