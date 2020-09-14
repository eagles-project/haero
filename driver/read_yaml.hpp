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

  const char* what() const throw() {
    return _message.c_str();
  }

  private:

  std::string _message;
};

/// This function reads simulation input from a YAML input file. If an error
/// is encountered, this throws a YamlException.
/// TODO: We could document the input spec here for clarity. For now, see
/// TODO: tests/smoke_test.yml for an example.
/// \param [in] filename The name of the file to be read.
/// \returns A vector containing one or more sets of simulation input read from
///          the file.
std::vector<SimulationInput> read_yaml(const std::string& filename);

}

#endif