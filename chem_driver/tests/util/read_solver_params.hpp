#ifndef HAERO_READ_SOLVER_PARAMS_HPP
#define HAERO_READ_SOLVER_PARAMS_HPP

#include "chem_driver/chem_driver.hpp"

namespace haero {
namespace chem_driver {

/// This exception class stores information about errors encountered in reading
/// data from a YAML file.
class YamlException : public std::exception {
 public:
  /// Constructs an exception containing the given descriptive message.
  YamlException(const std::string& message) : _message(message) {}

  /// Constructs an exception containing the given formatting string and
  /// C-style variadic arguments (a la printf).
  YamlException(const char* fmt, ...);

  const char* what() const throw() { return _message.c_str(); }

 private:
  std::string _message;
};

// struct SolverParams {
//   Real dtmin;
//   Real dtmax;
//   Real tend;
//   int max_time_iterations;
//   int max_newton_iterations;
//   Real atol_newton;
//   Real rtol_newton;
//   Real tol_time;
//   std::string outputfile;
//   // constructor that accepts the name of the yaml input file
//   SolverParams(const std::string& filename);
// };

}  // end namespace chem_driver
}  // end namespace haero
#endif
