#ifndef HAERO_YAML_FILE_HPP
#define HAERO_YAML_FILE_HPP

#include <vector>
#include <map>
#include <yaml-cpp/yaml.h>
#include "haero/mode.hpp"
#include "haero/species.hpp"
#include "haero/reaction.hpp"

namespace haero {

/// This exception class stores information about errors encountered in reading
/// data from a YAML file.
class Yaml_exception: public std::exception {
  public:

  /// Constructs an exception containing the given descriptive message.
  Yaml_exception(const std::string& message):
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
class Yaml_file {
  public:

  /// Given a filename, creates a YAML file from input can be read.
  /// @param [in] filename An absolute path to a YAML file to be read.
  Yaml_file(const std::string& filename);

  // Disallowed machinery.
  Yaml_file(const Yaml_file&) = delete;
  Yaml_file& operator=(const Yaml_file&) = delete;

  /// Destructor.
  ~Yaml_file();

  /// Reads data for all modes in this input file into a vector.
  /// @returns the vector of modes read from the file.
  std::vector<Mode> read_modes() const;

  /// Reads data for all species in this input file into a vector.
  /// @returns the vector of species read from the file.
  std::vector<Species> read_species() const;

  /// Reads data for all chemical reactions in this input file into a vector.
  /// @returns the vector of reactions read from the file.
  std::vector<Reaction> read_reactions() const;

  /// Reads data for all initial conditions in this input file into a vector. A
  /// set of initial conditions is a mapping of prognostic variable names onto
  /// real-valued initial quantities.
  /// @returns the vector of initial conditions read from the file.
  std::vector<std::map<std::string, Real> > read_initial_conditions() const;

  /// Reads the constant perturbation factor to be applied to all initial
  /// conditions in this input file (across all timesteps and members of
  /// ensembles). This factor is multipled by machine epsilon and added to
  /// each initial condition.
  /// @returns the perturbaton factor for initial conditions.
  Real read_perturb_factor() const;

  /// Reads simulation control parameters.
  /// \param [out] do_gaschem 1 if gas chemistry enabled, 0 otherwise
  /// \param [out] do_cloudchem 1 if cloud chemistry enabled, 0 otherwise
  /// \param [out] do_gasaerexch 1 if gas-aerosol exchange enabled, 0 otherwise
  /// \param [out] do_rename 1 if ??? enabled, 0 otherwise
  /// \param [out] do_newnuc 1 if nucleation(?) enabled, 0 otherwise
  /// \param [out] do_coag 1 if coagulation enabled, 0 otherwise
  /// \param [out] do_calcsize 1 if ??? enabled, 0 otherwise
  /// \param [out] num_unit ???
  /// \param [out] frac_unit ???
  /// \param [out] gas_unit ???
  void read_control_params(int& do_gaschem, int& do_cloudchem,
                           int& do_gasaerexch, int& do_rename, int& do_newnuc,
                           int& do_coag, int& do_calcsize, int& num_unit,
                           int& frac_unit, int& gas_unit);

  /// Reads "MET" (background?) atmospheric parameters.
  /// \param [out] temp air temperature [K]
  /// \param [out] press pressure [Pa]
  /// \param [out] RH_CLEA clear-sky RH [?]
  /// \param [out] hgt height above Earth's surface [m]
  /// \param [out] cld_frac cloud fraction [-]
  void read_met_input(Real& temp, Real& press, Real& RH_CLEA, Real& hgt,
                      Real& cld_frac);

  /// Reads data for all timestep sizes for simulations into a vector.
  /// @returns the vector of timestep sizes read from the file.
  std::vector<Real> read_timesteps() const;

  /// Reads the duration for all simulations within this file.
  Real read_duration() const;

  /// Reads output parameters.
  /// \param [out] dir the directory to which haero output files are written
  ///                  (default: ".")
  /// \param [out] prefix the prefix for haero output files (default: "mam_output")
  /// \param [out] freq the frequency at which output is recorded, expressed as
  ///                   the number of steps between recordings.
  void read_output_params(std::string& dir,
                          std::string& prefix,
                          int& freq) const;

  /// Reads the directory in which haero output files are written.
  std::string read_output_dir() const;

  private:

  // The root node for the YAML file.
  YAML::Node _root;
};

}

#endif
