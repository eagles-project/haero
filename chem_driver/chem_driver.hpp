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
// #include "tchem/TChem_Impl_RateOfProgress.hpp"
// #include "tchem/TChem_KineticModelData.hpp"
// #include "tchem/TChem_Util.hpp"

namespace haero {
namespace chem_driver {

using Real = haero::Real;
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

// /// This type contains the chemical species relevant to the simulation
// struct ChemicalSpecies {
//   /// name of reactant
//   std::string name;
//   /// initial value (in the specified units)
//   Real initial_value;
//   /// units
//   std::string units;
//   /// constructor
//   ChemicalSpecies(const std::string& name, Real initial_value,
//                   const std::string& units);
// };

// /// This type contains the environmental conditions for the simulation
// struct EnvironmentalConditions {
//   /// initial temperature and pressure
//   Real initial_temp, initial_pressure;
//   /// units
//   std::string units_temp, units_pressure;
//   /// constructors: the one with arguments is currently unused but here in case
//   EnvironmentalConditions(Real T0, const std::string& T_units, Real P0,
//                           const std::string& P_units);
//   EnvironmentalConditions() = default;
// };

// /// enumeration of reaction type: used for determining how rates are calculated
// enum ReactionType { arrhenius, troe };

// /// This type contains the chemical species relevant to the simulation
// struct Reaction {
//   /// reaction type, enumerated
//   ReactionType type;
//   /// name of reaction, comes from input file and used to create the enum
//   std::string type_str;
//   /// reactants and products that are vectors of maps between a string (name)
//   /// and the integer-valued coefficient for that species
//   std::map<std::string, Real> reactants;
//   std::map<std::string, Real> products;
//   /// coefficients that are used to calculate the reaction rate
//   std::map<std::string, Real> rate_coefficients;
//   /// constructor
//   Reaction(const std::string& type_str,
//            const std::map<std::string, Real>& reactants,
//            const std::map<std::string, Real>& products,
//            const std::map<std::string, Real>& rate_coefficients);
//   /// copy constructor: used when copying the reactions from SimulationInput
//   /// to ChemSolver
//   Reaction(const Reaction& rxn)
//       : type{rxn.type},
//         type_str{rxn.type_str},
//         reactants{rxn.reactants},
//         products{rxn.products},
//         rate_coefficients{rxn.rate_coefficients} {}
//   /// method for either getting the coefficient from SimulationInput
//   /// or using defaults
//   void get_or_default(const std::map<std::string, Real>& mrate_coefficients,
//                       const std::string& name);
// };

// /// This type contains the data required to initialize and run a simulation
// struct SimulationInput {
//   /// name of the input file used to create the input
//   /// used in ChemSolver constructor to get the TChem inputs
//   std::string input_file;
//   /// vector containing the species participating in the simulation
//   std::vector<ChemicalSpecies> species;
//   /// struct containing the environmental conditions (temp and pressure)
//   EnvironmentalConditions env_conditions;
//   /// vector containing the reactions occurring in the simulation
//   std::vector<Reaction> reactions;
// };

// /// these are the file paths for the inputs and output Tchem requires
// /// Note: these are currently created by the ChemSolver constructor and deleted
// /// by the destructor
// class ChemFiles {
//  public:
//   std::string chemFile = "chem.inp";
//   std::string thermFile = "therm.dat";
//   std::string outputFile = "omega.dat";
//   ChemFiles() = default;
// };

struct SolverParams {
  Real dtmin, dtmax;
  Real tbeg, tend;
  int num_time_iterations_per_interval;
  int max_time_iterations, max_newton_iterations;
  Real atol_newton, rtol_newton;
  Real atol_time, tol_time;
  int jacobian_interval;
  std::string outputfile;
  // constructor that accepts the name of the yaml input file
  void set_params(const std::string& filename, const bool& verbose);
  SolverParams() = default;
};

class ChemSolver {
 private:
  /// detail output flag
  bool verbose;
  /// flag for printing qoi's for debugging
  // FIXME: maybe can this when things are working
  bool print_qoi;
  FILE* fout;
  /// number of chemical batches/samples that will have the given composition
  int nbatch;
  // struct containing parameters that may be passed to the TChem time integrator
  SolverParams solver_params;
  /// 2D views containing chemical state on device and host
  Real_2d_view state;
  Real_2d_view_host state_host;
  // /// timer variables
  // Real t_device_batch;
  /// kokkos team policy
  policy_type policy;
  int team_size, vector_size;
  // /// temperature and pressure, used for calculating reaction rates
  // Real temperature, pressure;
  // std::string units_temp, units_pressure;
  // /// timer object for kokkos
  // Kokkos::Impl::Timer timer;
  // /// TChem kinetic model data
  TChem::KineticModelData kmd;
  // /// a const version of the object that contains the data describing the
  // /// kinetic model
  TChem::KineticModelNCAR_ConstData<Tines::UseThisDevice<exec_space>::type> kmcd;
  // /// vector of reactions, copied directly from the SimulationInput class
  // std::vector<Reaction> reactions;
  // /// *-----------------------------------------------------------------------*
  // /// prints a summary of the TChem run
  // void print_summary(const ChemFiles& cfiles);
  /// parses the yaml file for user-provided inputs required by TChem
  void parse_tchem_inputs(const std::string& input_file);
  // /// sets the reaction rates: currently called prior to get_tendencies() and
  // /// uses current temp and pressure, along with the user-provided rate
  // /// coefficients
  // void set_reaction_rates();

 public:
 //  /// ChemFiles object containing file paths for inputs/output
 //  ChemFiles cfiles;
  /// constructor
  ChemSolver(std::string input_file);
  /// destructor
  ~ChemSolver();
  void time_integrate();
 //  /// runs chemical model and returns the resulting tendencies
 //  Real_2d_view get_tendencies();
};

}  // end namespace chem_driver
}  // end namespace haero

#endif
