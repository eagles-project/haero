#include "read_chem_input.hpp"

#include <yaml-cpp/yaml.h>

#include <cstdarg>
#include <map>

namespace haero {
namespace chem_driver {

YamlException::YamlException(const char* fmt, ...) {
  char ss[256];
  va_list args;
  va_start(args, fmt);
  vsnprintf(ss, 255, fmt, args);
  va_end(args);
  _message.assign(ss);
}

/// anonymous namespace containing the read_xxx() functions that are used by
/// the read_chem_input() function to build a SimulationInput object
namespace {

/// Note that all of the read_xxx() functions work in more or less the same way.
/// As such, only the first is rigorously commented.

/// function to read a number of species and extract their name, initial value,
/// and units
std::vector<ChemicalSpecies> read_species(const YAML::Node& root) {
  // create a blank vector of ChemicalSpecies
  std::vector<ChemicalSpecies> species;
  // do we have a species section in the input file?
  if (root["species"] and root["species"].IsMap()) {
    auto node = root["species"];
    // loop over the provided number of species
    for (auto iter : node) {
      std::string name = iter.first.as<std::string>();
      auto mnode = iter.second;
      if (not mnode["initial_value"]) {
        throw YamlException(
            "species entry has no initial value (initial_value).");
      } else if (not mnode["units"]) {
        throw YamlException("species entry has no units (units).");
      } else {
        // call the ChemicalSpecies constructor with the necessary info from the
        // input yaml and insert it into the vector
        species.push_back(ChemicalSpecies(name,
                                          mnode["initial_value"].as<Real>(),
                                          mnode["units"].as<std::string>()));
      }
    }
  } else {
    throw YamlException("No species section was found!");
  }
  return species;
}

/// function to read the environmental conditions (temperature and pressure),
/// and extract their initial value and units
EnvironmentalConditions read_env_conditions(const YAML::Node& root) {
  EnvironmentalConditions env_conditions;
  if (root["environmental_conditions"] and
      root["environmental_conditions"].IsMap()) {
    auto node = root["environmental_conditions"];
    if (not node["temperature"]) {
      throw YamlException(
          "environmental_conditions contains no temperature "
          "entry (temperature).");
    }
    if (not node["pressure"]) {
      throw YamlException(
          "environmental_conditions contains no pressure "
          "entry (pressure).");
    } else {
      auto mnode = node["temperature"];
      env_conditions.initial_temp = mnode["initial_value"].as<Real>();
      env_conditions.units_temp = mnode["units"].as<std::string>();
      mnode = node["pressure"];
      env_conditions.initial_pressure = mnode["initial_value"].as<Real>();
      env_conditions.units_pressure = mnode["units"].as<std::string>();
    }
  } else {
    throw YamlException("No environmental_conditions section was found!");
  }
  return env_conditions;
}

/// function to read a number of reactions and extract the reactants, products,
/// and rate coefficients
/// Note that each of these fields is a map between a string value (name) and
/// the corresponding coefficient (stoichiometric or rate)
std::vector<Reaction> read_reactions(const YAML::Node& root) {
  std::vector<Reaction> reactions;
  std::map<std::string, Real> mreactants;
  std::map<std::string, Real> mproducts;
  std::map<std::string, Real> mcoeffs;
  // do we have a reaction node?
  if (root["reactions"] and root["reactions"].IsSequence()) {
    // define the reaction node
    auto reaction_node = root["reactions"];
    // loop over reaction nodes
    for (auto mnode : reaction_node) {
      if (not mnode["type"]) {
        throw YamlException("reactions entry has no initial type (type).");
      } else if (not mnode["reactants"]) {
        throw YamlException("reactions entry has no reactants (reactants).");
      } else if (not mnode["products"]) {
        throw YamlException("reactions entry has no products (products).");
      } else if (not mnode["rate_coefficients"]) {
        throw YamlException(
            "reactions entry has no rate_coefficients "
            "(rate_coefficients).");
      } else {
        if (mnode["reactants"] and mnode["reactants"].IsMap()) {
          auto react_node = mnode["reactants"];
          for (auto iter : react_node) {
            mreactants.insert(std::pair<std::string, Real>(
                iter.first.as<std::string>(), iter.second.as<Real>()));
          }
        }
        if (mnode["products"] and mnode["products"].IsMap()) {
          auto prod_node = mnode["products"];
          for (auto iter : prod_node) {
            mproducts.insert(std::pair<std::string, Real>(
                iter.first.as<std::string>(), iter.second.as<Real>()));
          }
        }
        if (mnode["rate_coefficients"] and mnode["rate_coefficients"].IsMap()) {
          auto rc_node = mnode["rate_coefficients"];
          for (auto iter : rc_node) {
            mcoeffs.insert(std::pair<std::string, Real>(
                iter.first.as<std::string>(), iter.second.as<Real>()));
          }
        }
        reactions.push_back(Reaction(mnode["type"].as<std::string>(),
                                     mreactants, mproducts, mcoeffs));
        // clear out the maps we provide to the constructor so there's no junk
        // in there for the next constructor call
        mreactants.clear();
        mproducts.clear();
        mcoeffs.clear();
      }
    }
  } else {
    throw YamlException("No reactions section was found!");
  }
  return reactions;
}

//  function to validate that the reactions only contain species that are given
/// in the species section
void validate_reactions(SimulationInput sim_inp) {
  // construct a list of species names
  std::vector<std::string> species_list;
  for (auto species : sim_inp.species) {
    species_list.push_back(species.name);
  }
  // loop over reactions
  for (auto rxn : sim_inp.reactions) {
    // determine if reactants are in the species list vector
    for (auto rxt : rxn.reactants) {
      // Note: find() returns vec.end() if it is not found
      if (std::find(species_list.begin(), species_list.end(), rxt.first) ==
          species_list.end()) {
        std::string str = "Reactant " + rxt.first + " not in species list";
        throw YamlException(str);
      }
    }
    // determine if products are in the species list vector
    for (auto prod : rxn.products) {
      if (std::find(species_list.begin(), species_list.end(), prod.first) ==
          species_list.end()) {
        std::string str = "Product " + prod.first + " not in species list";
        throw YamlException(str);
      }
    }
  }
}

}  // end anonymous namespace

/// function to read the chemistry input yaml file and construct a
/// SimulationInput object from what is found there
SimulationInput read_chem_input(const std::string& filename) {
  // Try to load the input from the yaml file
  try {
    auto root = YAML::LoadFile(filename);

    SimulationInput sim_inp;
    // save the filename for use by the ChemSolver
    sim_inp.input_file = filename;
    // get the required info from the file, using the above functions
    sim_inp.species = read_species(root);
    sim_inp.env_conditions = read_env_conditions(root);
    sim_inp.reactions = read_reactions(root);
    // ensure that reactions contain only products/reactants in the species list
    validate_reactions(sim_inp);

    return sim_inp;
  } catch (YAML::BadFile& e) {
    throw YamlException(e.what());
  } catch (YAML::ParserException& e) {
    throw YamlException(e.what());
  }
}

}  // end namespace chem_driver
}  // end namespace haero
