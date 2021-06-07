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
  std::string name;
  /// initial value (in the specified units)
  Real initial_value;
  /// units
  std::string units;
  /// constructor
  ChemicalSpecies(std::string name, Real initial_value, std::string units);
};

/// This type contains the environmental conditions for the simulation
struct EnvironmentalConditions {
  /// initial temperature and pressure
  Real initial_temp, initial_pressure;
  /// units
  std::string units_temp, units_pressure;
  /// constructors: the one with arguments is currently unused but here in case
  EnvironmentalConditions(Real T0, std::string T_units, Real P0,
                          std::string P_units);
  EnvironmentalConditions() = default;
};

/// enumeration of reaction type: used for determining how rates are calculated
enum ReactionType { arrhenius, troe };

/// This type contains the chemical species relevant to the simulation
struct Reaction {
  /// reaction type, enumerated
  ReactionType type;
  /// name of reaction, comes from input file and used to create the enum
  std::string type_str;
  /// reactants and products that are vectors of maps between a string (name)
  /// and the integer-valued coefficient for that species
  std::map<std::string, Real> reactants;
  std::map<std::string, Real> products;
  /// coefficients that are used to calculate the reaction rate
  std::map<std::string, Real> rate_coefficients;
  /// constructor
  Reaction(std::string type_str, std::map<std::string, Real> reactants,
           std::map<std::string, Real> products,
           std::map<std::string, Real> rate_coefficients);
  /// copy constructor: used when copying the reactions from SimulationInput
  /// to ChemSolver
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
  /// name of the input file used to create the input
  /// used in ChemSolver constructor to get the TChem inputs
  std::string input_file;
  /// vector containing the species participating in the simulation
  std::vector<ChemicalSpecies> species;
  /// struct containing the environmental conditions (temp and pressure)
  EnvironmentalConditions env_conditions;
  /// vector containing the reactions occurring in the simulation
  std::vector<Reaction> reactions;
};

/// these are the file paths for the inputs and output Tchem requires
/// Note: these are currently created by the ChemSolver constructor and deleted
/// by the destructor
class ChemFiles {
 public:
  std::string chemFile = "chem.inp";
  std::string thermFile = "therm.dat";
  std::string outputFile = "omega.dat";
  ChemFiles() = default;
};

class ChemSolver {
 private:
  /// verbose output flag
  bool verbose;
  /// number of chemical batches/samples that will have the given composition
  int nBatch;
  /// 2D views containing chemical state, reaction rates, and omega (results)
  real_type_2d_view state, kfor, krev, omega;
  /// timer variables
  Real t_device_batch;
  /// kokkos team policy
  policy_type policy;
  /// temperature and pressure, used for calculating reaction rates
  Real temperature, pressure;
  std::string units_temp, units_pressure;
  /// timer object for kokkos
  Kokkos::Impl::Timer timer;
  /// TChem kinetic model data
  TChem::KineticModelData kmd;
  /// a const version of the object that contains the data describing the
  /// kinetic model
  TChem::KineticModelConstData<TChem::exec_space> kmcd;
  /// vector of reactions, copied directly from the SimulationInput class
  std::vector<Reaction> reactions;
  /// *-----------------------------------------------------------------------*
  /// prints a summary of the TChem run
  void print_summary(const ChemFiles& cfiles);
  /// parses the yaml file for user-provided inputs required by TChem
  void parse_tchem_inputs(SimulationInput& sim_inp);
  /// sets the reaction rates: currently called prior to get_tendencies() and uses
  /// current temp and pressure, along with the user-provided rate coefficients
  void set_reaction_rates();

 public:
  /// ChemFiles object containing file paths for inputs/output
  ChemFiles cfiles;
  /// constructor
  ChemSolver(SimulationInput& sim_inp);
  /// destructor that removes the temporary input files required by TChem
  ~ChemSolver();
  /// runs chemical model and returns the resulting tendencies
  real_type_2d_view get_tendencies();
};

}  // end namespace chemDriver
}  // end namespace haero

// NOTE: everything in this namespace is copied directly from the TChem
// implementation of the toy problem
namespace from_tchem {

struct SourceTermToyProblem {
  template <typename KineticModelConstDataType>
  static inline ordinal_type getWorkSpaceSize(
      const KineticModelConstDataType& kmcd) {
    return SourceTermToyProblem::getWorkSpaceSize(kmcd);
  }

  static void runDeviceBatch(
      /// input
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
