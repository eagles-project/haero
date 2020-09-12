#ifndef HAERO_YAMLFILE_HPP
#define HAERO_YAMLFILE_HPP

#include <vector>
#include <map>
#include <yaml-cpp/yaml.h>
#include "haero/mode.hpp"
#include "haero/species.hpp"

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

/// This class defines an interface for reading data from a YAML file. The
/// interface supports an "ensemble mode" in which the parameters for several
/// related simulations may be defined succinctly within a single input file.
/// Each of the methods on this interface may throw a Yaml_exception.
class YamlFile {
  public:

  /// Given a filename, creates a YAML file from input can be read.
  /// @param [in] filename An absolute path to a YAML file to be read.
  YamlFile(const std::string& filename);

  // Disallowed machinery.
  YamlFile(const YamlFile&) = delete;
  YamlFile& operator=(const YamlFile&) = delete;

  /// Destructor.
  ~YamlFile();

  /// Reads data for all modes in this input file into a vector.
  /// @returns the vector of modes read from the file.
  std::vector<Mode> read_modes() const;

  /// Reads data for all aerosol species in this input file into a vector.
  /// @returns the vector of species read from the file.
  std::vector<Species> read_aerosol_species() const;

  /// Reads data for all gas species in this input file into a vector.
  /// @returns the vector of species read from the file.
  std::vector<Species> read_gas_species() const;

  /// Reads data for activating/deactivating various physical processes.
  /// @returns a map from process names to booleans (on/off settings).
  std::map<std::string, bool> read_physics_settings() const;

  /// Reads data for ambient atmospheric conditions. Returns the name of
  /// the model and a mapping of parameter names to their values.
  std::pair<std::string, std::map<std::string, Real> > read_atmosphere() const;

  /// Reads data for all initial conditions in this input file into a vector. A
  /// set of initial conditions is a mapping of prognostic variable names onto
  /// real-valued initial quantities.
  /// @returns the vector of initial conditions read from the file.
  std::vector<std::map<std::string, Real> > read_initial_conditions() const;

  /// Reads data for all timestep sizes for simulations into a vector.
  /// @returns the vector of timestep sizes read from the file.
  std::vector<Real> read_timesteps() const;

  /// Reads the duration for all simulations within this file.
  Real read_duration() const;

  /// Reads the directory and prefix for output files.
  std::pair<std::string, std::string> read_output_params() const;

  /// Reads the output frequency (in steps). Returns -1 if no frequency is
  /// given.
  int read_output_freq() const;

  private:

  // The root node for the YAML file.
  YAML::Node _root;
};

}

#endif
