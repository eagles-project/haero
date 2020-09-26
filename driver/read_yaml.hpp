#ifndef HAERO_READ_YAML_HPP
#define HAERO_READ_YAML_HPP

#include "driver.hpp"

namespace haero {

/// This exception class stores information about errors encountered in reading
/// data from a YAML file.
class YamlException: public std::exception {
  public:

  /// Constructs an exception containing the given descriptive message.
  YamlException(const std::string& message):
    _message(message) {}

  /// Constructs an exception containing the given formatting string and
  /// C-style variadic arguments (a la printf).
  YamlException(const char* fmt, ...);

  const char* what() const throw() {
    return _message.c_str();
  }

  private:

  std::string _message;
};

/// This function reads simulation input from a YAML input file. If an error
/// is encountered, this throws a YamlException. The input spec for the YAML
/// file is documented in Haero's design document.
/// \param [in] filename The name of the file to be read.
/// \returns A vector containing one or more sets of simulation input read from
///          the file.
std::vector<SimulationInput> read_yaml(const std::string& filename);

}

#endif
