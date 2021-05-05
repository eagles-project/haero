#include "chemDriver.hpp"
// #include "haero/model.hpp"

namespace haero {
namespace chemDriver {

ChemicalSpecies::ChemicalSpecies(std::string mname,
                                 Real minitial_value,
                                 std::string munits)
                                 :
                                 name(mname),
                                 initial_value(minitial_value),
                                 units(munits){}

EnvironmentalConditions::EnvironmentalConditions(std::string mname,
                                                 Real minitial_value,
                                                 std::string munits)
                                                 :
                                                 name(mname),
                                                 initial_value(minitial_value),
                                                 units(munits){}
Reaction::Reaction(std::string mtype_str,
                   std::map<std::string, Real> mreactants,
                   std::map<std::string, Real> mproducts,
                   std::map<std::string, Real> mrate_coefficients)
                   :
                   type_str(mtype_str),
                   reactants(mreactants),
                   products(mproducts){
// convert reaction type string to lowercase
transform(type_str.begin(), type_str.end(), type_str.begin(), ::tolower);
// use the type string to set the enumerated type for the reaction and assign
// the rate coefficients, either based on input or to default values
if (type_str.compare("arrhenius") == 0)
{
  type = arrhenius;
  rate_coefficients["A"]  = (mrate_coefficients["A"]) ?
                            mrate_coefficients["A"] : 1.0;
  rate_coefficients["Ea"] = (mrate_coefficients["Ea"]) ?
                            mrate_coefficients["Ea"] : 0.0;
  rate_coefficients["B"]  = (mrate_coefficients["B"]) ?
                            mrate_coefficients["B"] : 0.0;
  rate_coefficients["D"]  = (mrate_coefficients["D"]) ?
                            mrate_coefficients["D"] : 300.0;
  rate_coefficients["E"]  = (mrate_coefficients["E"]) ?
                            mrate_coefficients["E"] : 0.0;
}else if (type_str.compare("troe") == 0)
{
  type = troe;
  rate_coefficients["k0_A"]   = (mrate_coefficients["k0_A"]) ?
                                mrate_coefficients["k0_A"] : 1.0;
  rate_coefficients["k0_B"]   = (mrate_coefficients["k0_B"]) ?
                                mrate_coefficients["k0_B"] : 0.0;
  rate_coefficients["k0_C"]   = (mrate_coefficients["k0_C"]) ?
                                mrate_coefficients["k0_C"] : 0.0;
  rate_coefficients["kinf_A"] = (mrate_coefficients["kinf_A"]) ?
                                mrate_coefficients["kinf_A"] : 1.0;
  rate_coefficients["kinf_B"] = (mrate_coefficients["kinf_B"]) ?
                                mrate_coefficients["kinf_B"] : 0.0;
  rate_coefficients["kinf_C"] = (mrate_coefficients["kinf_C"]) ?
                                mrate_coefficients["kinf_C"] : 0.0;
  rate_coefficients["Fc"]     = (mrate_coefficients["Fc"]) ?
                                mrate_coefficients["Fc"] : 0.6;
  rate_coefficients["N"]      = (mrate_coefficients["N"]) ?
                                mrate_coefficients["N"] : 1.0;
}else {
  std::cout << "ERROR: reaction type currently unsupported." << "\n";
}
}

// void haero_driver(const std::vector<SimulationInput>& ensemble)
// {

// }

} // namespace chemDriver
} // namespace haero

