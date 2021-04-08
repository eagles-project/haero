#ifndef HAERO_SW_PARSE_YAML_HPP
#define HAERO_SW_PARSE_YAML_HPP

#include "haero/haero.hpp"
#include "haero/modal_aerosol_config.hpp"

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

// Data structure that stores parameter walking information.
struct ParameterWalk {
  // aerosol configuration
  haero::ModalAerosolConfig aero_config;

  // process name
  std::string process;

  // timestepping parameters
  haero::Real dt;
  int nsteps;

  // parameters to walk (name -> vector of values)
  std::map<std::string, std::vector<haero::Real>> parameters;

  // atmospheric state parameters
  haero::Real temperature, pressure, relative_humidity, height, cloud_fraction;

  // aerosol initial data
  std::vector<haero::Real> number_concs;
  std::vector<std::vector<haero::Real>> aero_mix_fractions;
  std::vector<haero::Real> gas_mix_fractions;
};

/// This function reads input from a YAML input file. If an error is
/// encountered, this throws a YamlException.
/// \param [in] aerosol_config A selected modal aerosol configuration that
///                            defines how the input is interpreted.
/// \param [in] filename The name of the file to be read
/// \returns A set of data and metadata used by skywalker to run a parameter
///          study
ParameterWalk parse_yaml(const haero::ModalAerosolConfig& aerosol_config,
                         const std::string& filename);

#endif
