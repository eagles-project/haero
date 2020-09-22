#ifndef HAERO_NCFILE_HPP
#define HAERO_NCFILE_HPP

#include "haero/haero_config.hpp"
// #include "column_base.hpp" // requires column_dynamics branch
// #include "dyn_column.hpp" // requires column_dynamics branch
#include <string>
//#include "ekat/util/ekat_units.hpp" // requires updated ekat
#include "netcdf.h"

namespace haero {

class NcFile {
  public:
    /** @brief Constructor.  Must be called from host.

      Assumes at most 1 file open per NcFile instance.
      All values initialized to invalid data; they will be set to valid data by other
      member functions.

      Writes HAERO metadata (version, git hash) to file.

      @todo Require new_filename has extension .nc

      @throws

      @param new_filename name of NetCDF data file to create.
    */
    NcFile(const std::string& new_filename) : fname(new_filename),
     ncid(NC_EBADID), level_dimid(NC_EBADID), interface_dimid(NC_EBADID),
     mode_dimid(NC_EBADID), time_dimid(NC_EBADID), ndims(0), nvars(0) {
      open();
      add_time_dim();
    }

    /** @brief String with basic info about this NcFile.

      @param tab_level indent level for console output.
    */
    std::string info_string(const int& tab_level=0) const;

    /** @brief Adds 2 dimensions to an existing NcFile; one for level midpoints,
      one for level interfaces.

      @throws

      This is required before any column variables can be defined.
      @param nlev number of vertical levels (midpoints) in a column.
    */
    void add_level_dims(const int& nlev);


    /** @brief Adds a dimension for the number of modes.

      @throws

      This is required before any variables indexed by mode can be defined.
      @param nmodes number of modes in a HAERO aerosol computation.
    */
    void add_mode_dim(const int& nmodes);


    /** @brief Close the data file associated with this NcFile.

      @throws

      @warning Data may not be written to file until this method is called.
    */
    void close();

    /** @brief return the current number of dimensions
    */
    inline int get_ndims() const {return ndims;}

    /** @brief return the filename
    */
    inline std::string get_fname() const {return fname;}

    /** @brief return the current number of variables
    */
    inline int get_nvars() const {return nvars;}

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

    /** @brief Adds a dimension for time.

      Called by the constructor, should not be called by the user.

      @throws

      This is required before any variables indexed by time step can be defined.
    */
    void add_time_dim();

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
