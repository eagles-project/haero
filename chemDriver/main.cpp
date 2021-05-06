// This is Haero's chemistry batch reactor driver. It reads a YAML file
// containing input parameters for a simulation (or an ensemble of simulations)
// and solves the chemical system using TChem.

#include <cstdio>
// #include <cstdlib>
// #include <limits>
#include "haero/haero.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/ekat_session.hpp"
#include "chemUtil.hpp"
#include "read_chem_input.hpp"

// ***TESTING***
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;
using namespace haero;
using namespace haero::chemDriver;

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
    fprintf(stderr, "WARNING: haero was built with uncommitted changes!\n");
    fprintf(stderr, "WARNING: DO NOT USE FOR SCIENCE!\n");
  }
}

// Print driver usage information and exit.
void usage(const char* exe)
{
  fprintf(stderr, "Too few inputs to %s--usage:\n", exe);
  fprintf(stderr, "%s <input.yml>\n", exe);
  exit(1);
}

} // anonymous namespace

int main(int argc, const char** argv)
{
  ekat::initialize_ekat_session(argc, const_cast<char**>(argv), true);
  print_banner();

  if (argc < 2)
    usage(argv[0]);

    // Write out some test data to our current working directory.
    const char* chem_inp = R"INPUT(ELEMENTS
X /1/
END
SPECIES
X X2
END
THERM ALL
    300.000  1000.000  5000.000
X                        X  1               G   200.000  6000.000 1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 4.37967000E+00 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 4.37967000E+00                   4
X2                       X  2               G   200.000  6000.000 1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 4.37967000E+00 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 4.37967000E+00                   4
END
REACTIONS
X2=>2X      1E+0    1   1
X+X=>X2      1E+0    1   1
END
    )INPUT";
    mkdir("tests", 0777);
    mkdir("tests/toy-problem", 0777);
    FILE* f = fopen("tests/toy-problem/chem.inp", "w");
    fprintf(f, "%s", chem_inp);
    fclose(f);
    f = fopen("tests/toy-problem/therm.dat", "w");
    fclose(f);

  // Read the input file and extract input.
  std::string input_file(argv[1]);
  SimulationInput sim_input = read_chem_input(input_file);

  // lat/lon for sun zenith
  const Real latz(20.0); // ! degrees
  const Real lonz(300.0);// ! degrees

  // lat/lon for column position
  Real lat = 20.0;
  Real lon = 37.5;
  // initial concentrations for the two species
  Real initX = 1.0e-6;
  Real initX2 = 1.0e-6;

  // calculate k1 reaction rate, based on position of column and sun's zenith
  Real k1 = sin(lat*PI()/180) * sin(latz*PI()/180) +
                       cos(lat*PI()/180) * cos(latz*PI()/180) *
                       cos(lon*PI()/180 - lonz*PI()/180);
  k1 = k1 > 0 ? k1 : 0;
  Real k2 = 1;
  // Kokkos::initialize();
  // provide arbitrarily chosen inputs
  chemUtil::chemSolver chem_solver("toy-problem/", false, 1, false, k1, k2,
                                   initX, initX2);

  // run the problem on device
  real_type_2d_view results = chem_solver.get_results();

  // create mirror view and deep copy to host
  auto results_host = Kokkos::create_mirror_view(results);
  Kokkos::deep_copy(results_host, results);

  // eq (4) tendency should be positive, eq (5) negative
  const Real val00 = results_host(0, 0);
  const Real val01 = results_host(0, 1);
  std::cout << "val00, val01 = " << val00 << ", " << val01 << "\n";
  // Kokkos::finalize();

  // // Now run the driver with our input.
  // haero_driver(ensemble);

  // FIXME: this causes the following error:
    // libc++abi.dylib: terminating with uncaught exception of type
    // std::runtime_error: Kokkos allocation "KMD::sPhase" is being deallocated
    // after Kokkos::finalize was called
  // // Clean up.
  // ekat::finalize_ekat_session();
  return 0;
}
