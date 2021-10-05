#include "read_solver_params.hpp"

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

// /// function to read the chemistry input yaml file and construct a
// /// SolverParams struct from what is found there
// SolverParams::SolverParams(const std::string& filename) {
//   // load the input from the yaml file
//   auto root = YAML::LoadFile(filename);
//   if (root["solver_parameters"] and root["solver_parameters"].IsMap()) {
//     auto node = root["solver_parameters"];
//     if (not node["dtmin"]) {
//       throw YamlException(
//           "solver_parameters contains no dtmin "
//           "entry (dtmin).");
//     }
//     if (not node["dtmax"]) {
//       throw YamlException(
//           "solver_parameters contains no dtmax "
//           "entry (dtmax).");
//     }
//     if (not node["tend"]) {
//       throw YamlException(
//           "solver_parameters contains no tend "
//           "entry (tend).");
//     }
//     if (not node["max_time_iterations"]) {
//       throw YamlException(
//           "solver_parameters contains no max_time_iterations "
//           "entry (max_time_iterations).");
//     }
//     if (not node["max_newton_iterations"]) {
//       throw YamlException(
//           "solver_parameters contains no max_newton_iterations "
//           "entry (max_newton_iterations).");
//     }
//     if (not node["atol_newton"]) {
//       throw YamlException(
//           "solver_parameters contains no atol_newton "
//           "entry (atol_newton).");
//     }
//     if (not node["rtol_newton"]) {
//       throw YamlException(
//           "solver_parameters contains no rtol_newton "
//           "entry (rtol_newton).");
//     }
//     if (not node["tol_time"]) {
//       throw YamlException(
//           "solver_parameters contains no tol_time "
//           "entry (tol_time).");
//     }
//     if (not node["outputfile"]) {
//       throw YamlException(
//           "solver_parameters contains no outputfile "
//           "entry (outputfile).");
//     } else {
//       dtmin = node["dtmin"].as<Real>();
//       dtmax = node["dtmax"].as<Real>();
//       tend = node["tend"].as<Real>();
//       max_time_iterations = node["max_time_iterations"].as<int>();
//       max_newton_iterations = node["max_newton_iterations"].as<int>();
//       atol_newton = node["atol_newton"].as<Real>();
//       rtol_newton = node["rtol_newton"].as<Real>();
//       tol_time = node["tol_time"].as<Real>();
//       outputfile = node["outputfile"].as<std::string>();
//     }
//   } else {
//     throw YamlException("No solver_parameters section was found!");
//   }
// }

}  // end namespace chem_driver
}  // end namespace haero
