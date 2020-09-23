#ifndef HAERO_NCFILE_HPP
#define HAERO_NCFILE_HPP

#include "haero/haero_config.hpp"
// #include "column_base.hpp" // requires column_dynamics branch
// #include "dyn_column.hpp" // requires column_dynamics branch
#include <string>
#include "ekat/util/ekat_units.hpp" // requires updated ekat
#include "netcdf.h"

namespace haero {

/** @brief This class reads NetCDF files written by haero's driver.
 */
class NcReader {
  public:
    /** @brief Constructor.  Must be called from host.

      Opens the NetCDF file with the given name for reading, and validates its
      format.

      @param filename name of NetCDF data file to open.
    */
    NcReader(const std::string& filename);

    /** @brief Reads a variable with the given name from the NetCDF file.
      @param var_name name of the variable to read from the file.
      @param data a vector into which the data is read.
     */
    void read_var(const std::string& var_name,
                  std::vector<Real>& data);

  protected:
    /** @brief `netcdf.h` defines numerous error codes (integers); this function
      decodes this integers and throws an exeption with (hopefully) a more helpful message.

      @throws
    */
    void handle_errcode(const int& ec) const;

    /** @brief Calls nc_create to make a new .nc file with appropriate
    parameters from `netcdf.h` (NETCDF4 format, NC_CLOBBER for overwrite), etc.

      Called by the constructor, should not be called by the user.
    */
    void open();

    /// filename for .nc data
    std::string fname;

    /// NetCDF file ID.  Assigned by `nc_create` during open()
    int ncid;
    /// Level dimension ID. Assigned by add_level_dims.
    int level_dimid;
    /// Interface dimension ID.  Assigned by add_level_dims.
    int interface_dimid;
    /// Mode dimension ID. Assigned by add_mode_dim.
    int mode_dimid;
    /// Time dimension ID.  Assigned by add_time_dim.
    int time_dimid;
    /// Tracks the number of NetCDF dimensions in the current file.
    int ndims;
    /// Tracks the number of NetCDF variables in the current file.
    int nvars;
};

}
#endif
