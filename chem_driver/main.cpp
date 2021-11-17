// This is Haero's chemistry batch reactor driver. It reads a YAML file
// containing input parameters for a simulation (or an ensemble of simulations)
// and solves the chemical system using TChem.

#include <cstdio>

// #include "chem_driver.hpp"
// #include "ekat/ekat_parameter_list.hpp"
// #include "ekat/ekat_session.hpp"
// #include "haero/haero.hpp"
// #include "haero/physical_constants.hpp"
// #include "read_chem_input.hpp"
// #include "tests/toy-problem/toy_problem.hpp"

// using namespace std;
// using namespace haero;
// using namespace haero::chem_driver;

// namespace {

// // Print version banner.
// void print_banner() {
//   string version = haero::version();
//   string revision = haero::revision();
//   bool not_for_science = haero::has_uncommitted_changes();
//   printf("haero v%s (git revision %s)\n", version.c_str(), revision.c_str());
//   if (not_for_science) {
//     fprintf(stderr, "WARNING: haero was built with uncommitted changes!\n");
//     fprintf(stderr, "WARNING: DO NOT USE FOR SCIENCE!\n");
//   }
// }

// // Print driver usage information and exit.
// void usage(const char* exe) {
//   fprintf(stderr, "Too few inputs to %s--usage:\n", exe);
//   fprintf(stderr, "%s <input.yml>\n", exe);
//   exit(1);
// }

// }  // end anonymous namespace

// // this is a work in progress--for now we just run the toy problem and print
// // the results to screen, as a prototype for usage
int main(int argc, const char** argv) {
//   print_banner();

//   if (argc < 2) usage(argv[0]);

//   // Read the input file and extract input.
//   std::string input_file(argv[1]);
//   SimulationInput sim_input = read_chem_input(input_file);

//   ekat::initialize_ekat_session(argc, const_cast<char**>(argv), true);
//   {
//     // this creates the toy-problem-specific input file for TChem
//     // FIXME: learn more about this and whether they've switched to the yaml
//     // input spec
//     create_chem_files();
//     ChemSolver chem_solver(sim_input);

//     // run the problem on device
//     real_type_2d_view results = chem_solver.get_tendencies();

//     // create mirror view and deep copy to host
//     auto results_host = Kokkos::create_mirror_view(results);
//     Kokkos::deep_copy(results_host, results);

//     // get the values and print to screen
//     const Real val00 = results_host(0, 0);
//     const Real val01 = results_host(0, 1);
//     const Real val10 = results_host(1, 0);
//     const Real val11 = results_host(1, 1);
//     std::cout << "val00, val01 = " << val00 << ", " << val01 << "\n";
//     std::cout << "val10, val11 = " << val10 << ", " << val11 << "\n";
//   }
//   ekat::finalize_ekat_session();
  return 0;
}
