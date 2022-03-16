#ifndef HAERO_CHEM_DRIVER_HPP
#define HAERO_CHEM_DRIVER_HPP

#include <yaml-cpp/yaml.h>

#include <map>
#include <string>
#include <vector>

#include "haero/haero.hpp"
#include "tchem/TChem_KineticModelData.hpp"
#include "tchem/TChem_Util.hpp"
#include "tchem/TChem_AtmosphericChemistry.hpp"

namespace haero {
namespace chem_driver {

// using Real = haero::Real;
using ordinal_type = TChem::ordinal_type;
using time_advance_type = TChem::time_advance_type;

using Real_0d_view = TChem::real_type_0d_view;
using Real_1d_view = TChem::real_type_1d_view;
using Real_2d_view = TChem::real_type_2d_view;

using time_advance_type_0d_view = TChem::time_advance_type_0d_view;
using time_advance_type_1d_view = TChem::time_advance_type_1d_view;

using Real_0d_view_host = TChem::real_type_0d_view_host;
using Real_1d_view_host = TChem::real_type_1d_view_host;
using Real_2d_view_host = TChem::real_type_2d_view_host;

using time_advance_type_0d_view_host = TChem::time_advance_type_0d_view_host;
using time_advance_type_1d_view_host = TChem::time_advance_type_1d_view_host;

using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

struct SolverParams {
  double dtmin, dtmax;
  double tbeg, tend;
  int num_time_iterations_per_interval;
  int max_time_iterations, max_newton_iterations;
  double atol_newton, rtol_newton;
  double atol_time, tol_time;
  int jacobian_interval;
  std::string outputfile;
  // constructor that accepts the name of the yaml input file
  void set_params(const std::string& filename, const bool& verbose);
  SolverParams() = default;
};

class ChemSolver {
 private:
  /// detail output flag
  bool verbose_;
  /// flag for printing qoi's--useful for debugging
  bool print_qoi_;
  FILE* fout_;
  /// number of chemical batches/samples that will have the given composition
  int nbatch_;
  // struct containing parameters that may be passed to the TChem time integrator
  SolverParams solver_params_;
  /// 2D views containing chemical state on device and host
  Real_2d_view state_;
  Real_2d_view_host state_host_;
  /// kokkos team policy
  policy_type policy_;
  int team_size_, vector_size_;
  // /// TChem kinetic model data
  TChem::KineticModelData kmd_;
  // /// a const version of the object that contains the data describing the
  // /// kinetic model
  TChem::KineticModelNCAR_ConstData<Tines::UseThisDevice<exec_space>::type> kmcd_;
  /// parses the yaml file for user-provided inputs required by TChem
  void parse_tchem_inputs_(const std::string& input_file);
 public:
  /// constructor
  explicit ChemSolver(std::string input_file);
  /// destructor
  ~ChemSolver();
  void time_integrate();
  void time_integrate(const double& tbeg, const double& tend);
  Real_2d_view& get_state() { return state_; };
};

}  // end namespace chem_driver
}  // end namespace haero

#endif
