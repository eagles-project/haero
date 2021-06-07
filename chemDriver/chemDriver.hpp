#ifndef HAERO_CHEM_DRIVER_HPP
#define HAERO_CHEM_DRIVER_HPP

#include <yaml-cpp/yaml.h>

#include <map>
#include <string>
#include <vector>

#include "haero/haero.hpp"
#include "tchem/TChem_Impl_RateOfProgress.hpp"
#include "tchem/TChem_KineticModelData.hpp"
#include "tchem/TChem_Util.hpp"

namespace haero {
namespace chemDriver {

// some aliases
using ordinal_type = TChem::ordinal_type;
using Real = haero::Real;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

/// This type contains the chemical species relevant to the simulation
struct ChemicalSpecies {
  /// name of reactant
  /// FIXME(maybe): should there be a long name and short name here?
  std::string name;
  /// initial value (in the specified units)
  Real initial_value;
  /// units
  std::string units;
  ChemicalSpecies(std::string name, Real initial_value, std::string units);
};

/// This type contains the environmental conditions for the simulation
struct EnvironmentalConditions {
  Real initial_temp, initial_pressure;
  std::string units_temp, units_pressure;
  EnvironmentalConditions(Real T0, std::string T_units, Real P0,
                          std::string P_units);
  EnvironmentalConditions() = default;
};

enum ReactionType { arrhenius, troe };

/// This type contains the chemical species relevant to the simulation
struct Reaction {
  /// reaction type, enumerated
  ReactionType type;
  std::string type_str;
  /// reactants and products that are vectors of maps between a string (name)
  /// and the integer-valued coefficient for that species
  std::map<std::string, Real> reactants;
  std::map<std::string, Real> products;
  // these are the coefficients that are used to calculate the reaction rate
  std::map<std::string, Real> rate_coefficients;
  Reaction(std::string type_str, std::map<std::string, Real> reactants,
           std::map<std::string, Real> products,
           std::map<std::string, Real> rate_coefficients);
  // copy constructor
  Reaction(const Reaction& rxn) {
    type = rxn.type;
    type_str = rxn.type_str;
    reactants = rxn.reactants;
    products = rxn.products;
    rate_coefficients = rxn.rate_coefficients;
  };
};

/// This type contains the data required to initialize and run a simulation
struct SimulationInput {
  std::string input_file;
  /// vector containing the species participating in the simulation
  std::vector<ChemicalSpecies> species;
  /// struct containing the environmental conditions
  EnvironmentalConditions env_conditions;
  /// vector containing the reactions occurring in the simulation
  std::vector<Reaction> reactions;
};

// these are the file paths for the inputs and output Tchem requires
// Note: these are currently created by the ChemSolver constructor and deleted
// by the destructor
class ChemFiles {
 public:
  std::string chemFile = "chem.inp";
  std::string thermFile = "therm.dat";
  std::string outputFile = "omega.dat";
  ChemFiles() = default;
};

class ChemSolver {
 private:
  // verbose/detail output flags
  bool verbose, detail;
  // timer variables
  Real t_device_batch;
  // kokkos team policy
  policy_type policy;
  // temperature and pressure, used for calculating reaction rates
  Real temperature, pressure;
  std::string units_temp, units_pressure;
  // timer object for kokkos
  Kokkos::Impl::Timer timer;
  // TChem kinetic model data
  TChem::KineticModelData kmd;
  // a const version of the object that contains the data describing the
  // kinetic model
  TChem::KineticModelConstData<TChem::exec_space> kmcd;
  // vector of reactions, copied directly from the SimulationInput class
  std::vector<Reaction> reactions;
  // *-----------------------------------------------------------------------*
  // prints a summary of the TChem run
  void print_summary(const ChemFiles& cfiles);
  void parse_tchem_inputs(SimulationInput& sim_inp);
  void set_reaction_rates();

 public:
  // ChemFiles object containing file paths for inputs/output
  ChemFiles cfiles;
  // number of chemical batches/samples that will have the given composition
  int nBatch;
  // 2D views containing chemical state, reaction rates, and omega (results)
  real_type_2d_view state, kfor, krev, omega;
  // constructor
  ChemSolver(SimulationInput& sim_inp);
  ~ChemSolver();
  // runs chemical model and saves the results (tendencies) to the output file
  real_type_2d_view get_results();
};

}  // namespace chemDriver
}  // namespace haero

// NOTE: everything in this namespace is copied directly from TChem
namespace from_tchem {

struct SourceTermToyProblem {
  template <typename KineticModelConstDataType>
  static inline ordinal_type getWorkSpaceSize(
      const KineticModelConstDataType& kmcd) {
    return SourceTermToyProblem::getWorkSpaceSize(kmcd);
  }

  //
  static void runDeviceBatch(  /// input
      typename UseThisTeamPolicy<exec_space>::type& policy,
      const real_type_2d_view& kfor, const real_type_2d_view& krev,
      const real_type_2d_view& state,
      /// output
      const real_type_2d_view& SourceTermToyProblem,
      /// const data from kinetic model
      const KineticModelConstDataDevice& kmcd);

};  // end STTP struct

template <typename KineticModelConstDataType>
KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize(
    const KineticModelConstDataType& kmcd);

}  // namespace from_tchem

#endif
