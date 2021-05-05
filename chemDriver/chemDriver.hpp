#ifndef HAERO_CHEM_DRIVER_HPP
#define HAERO_CHEM_DRIVER_HPP

#include <string>
#include <map>
#include <vector>

// #include "haero/mode.hpp"
// #include "haero/aerosol_species.hpp"
// #include "haero/haero_config.hpp"
// #include "haero/physical_constants.hpp"
// #include "host_params.hpp"
#include "haero/haero.hpp"

namespace haero {
namespace chemDriver {

/// This type contains the chemical species relevant to the simulation
struct ChemicalSpecies{
  /// long (descriptive) name
  std::string name;
  /// initial value (in the specified units)
  Real initial_value;
  /// units
  std::string units;
  ChemicalSpecies(std::string name, Real initial_value, std::string units);
};

/// This type contains the environmental conditions for the simulation
struct EnvironmentalConditions{
  /// long (descriptive) name
  std::string name;
  /// initial value (in the specified units)
  Real initial_value;
  /// units
  std::string units;
  EnvironmentalConditions(std::string name, Real initial_value, std::string units);
};

enum ReactionType {arrhenius, troe};

/// This type contains the chemical species relevant to the simulation
struct Reaction{
  /// reaction type, enumerated
  ReactionType type;
  std::string type_str;
  /// reactants and products that are vectors of maps between a string (name)
  /// and the integer-valued coefficient for that species
  std::map<std::string, Real> reactants;
  std::map<std::string, Real> products;
  // these are the coefficients that are used to calculate the reaction rate
  std::map<std::string, Real> rate_coefficients;
  Reaction(std::string type_str,
           std::map<std::string, Real> reactants,
           std::map<std::string, Real> products,
           std::map<std::string, Real> rate_coefficients);
};

/// This type contains the data required to initialize and run a simulation
struct SimulationInput{
  std::string input_file;
  /// vector containing the species participating in the simulation
  std::vector<ChemicalSpecies> species;
  /// struct containing the environmental conditions
  std::vector<EnvironmentalConditions> env_conditions;
  /// vector containing the reactions occurring in the simulation
  std::vector<Reaction> reactions;
};

// /// This is the entry point for the haero driver, which executes the simulations
// /// in the given ensemble.
// void haero_driver(const std::vector<SimulationInput>& ensemble);

} // namespace chemDriver
} // namespace haero

#endif
