#ifndef HAERO_READ_CHEM_INPUT_HPP
#define HAERO_READ_CHEM_INPUT_HPP

#include "chemDriver.hpp"

namespace haero {
namespace chemDriver {

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

/// This function reads simulation input from a YAML input file. If an error
/// is encountered, this throws a YamlException.
/// \param [in] filename The name of the file to be read.
/// \returns SimulationInput object that is meant to be passed to ChemSolver
SimulationInput read_chem_input(const std::string& filename);

}  // end namespace chemDriver
}  // end namespace haero
#endif
