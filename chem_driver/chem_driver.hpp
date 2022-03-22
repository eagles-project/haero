#ifndef HAERO_CHEM_DRIVER_HPP
#define HAERO_CHEM_DRIVER_HPP

#include <yaml-cpp/yaml.h>

#include <cstdarg>
#include <map>
#include <string>
#include <vector>

#include "haero/haero.hpp"
#include "tchem/TChem_AtmosphericChemistry.hpp"
#include "tchem/TChem_KineticModelData.hpp"
#include "tchem/TChem_Util.hpp"

namespace haero {
namespace chem_driver {

using Real = haero::Real;
using ordinal_type = TChem::ordinal_type;
using time_advance_type = TChem::time_advance_type;

using Real_1d_view = TChem::real_type_1d_view;
using Real_2d_view = TChem::real_type_2d_view;

using time_advance_type_1d_view = TChem::time_advance_type_1d_view;

using Real_1d_view_host = TChem::real_type_1d_view_host;
using Real_2d_view_host = TChem::real_type_2d_view_host;

using time_advance_type_1d_view_host = TChem::time_advance_type_1d_view_host;

using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

// these are various parameters that are required by the TChem solver
struct SolverParams {
  // Real dtmin, dtmax;
  // Real tbeg, tend;
  // int num_time_iterations_per_interval;
  int max_time_iterations;
  // int max_newton_iterations;
  Real atol_newton, rtol_newton;
  Real atol_time, tol_time;
  // int jacobian_interval;
  std::string outputfile;
  void set_params(const std::string& filename, const bool& verbose);
  SolverParams() = default;
  // tchem's struct for holding time-stepping info
  time_advance_type tadv_default;
};

class ChemSolver {
 private:
  // detail output flag
  bool verbose_;
  // flag for printing qoi's--useful for debugging
  bool print_qoi_;
  // the output file for the time series of system states
  FILE* fout_;
  // number of chemical batches/samples that will have the given composition
  int nbatch_;
  // struct containing parameters that may be passed to the TChem time
  // integrator
  SolverParams solver_params_;
  // 2D views containing chemical state on device and host
  Real_2d_view state_;
  Real_2d_view_host state_host_;
  // kokkos team policy
  policy_type policy_;
  int team_size_, vector_size_;
  // TChem kinetic model data
  TChem::KineticModelData kmd_;
  // a const version of the object that contains the data describing the
  // kinetic model
  TChem::KineticModelNCAR_ConstData<Tines::UseThisDevice<exec_space>::type>
      kmcd_;
  // parses the yaml file for user-provided inputs required by TChem
  void parse_tchem_inputs_(const std::string& input_file);

 public:
  // constructor
  explicit ChemSolver(std::string input_file);
  // destructor
  ~ChemSolver();
  // different time integrator options to get begin and end times from input
  // file, or specify them (allowing external time stepping)
  void time_integrate();
  void time_integrate(const Real& tbeg, const Real& tend);
  // returns the current state of the simulation
  Real_2d_view& get_state() { return state_; };
};

}  // end namespace chem_driver
}  // end namespace haero

#endif
