#ifndef HAERO_NCFILE_HPP
#define HAERO_NCFILE_HPP

#include "haero/haero_config.hpp"
// #include "column_base.hpp" // requires column_dynamics branch
// #include "dyn_column.hpp" // requires column_dynamics branch
#include <string>
#include "ekat/util/ekat_units.hpp" // requires updated ekat
#include "netcdf.h"

namespace haero {
namespace driver {

/** @brief This class reads NetCDF files written by haero's driver.
 */
class NcReader final {
  public:
    /// Opens the NetCDF file with the given name for reading, and validates
    /// its contents.
    /// @param filename name of NetCDF data file to open.
    NcReader(const std::string& filename);

    /// Destructor.
    ~NcReader();

    /// Reads a variable defined at grid levels (column cells) with the given
    /// name from the NetCDF file.
    /// @param var_name name of the variable to read from the file.
    /// @param time_index the time index associated with the value.
    /// @param data a vector into which the data is read. The vector is resized
    ///             as necessary.
    void read_level_var(const std::string& var_name,
                        size_t time_index,
                        std::vector<Real>& data);

    /// Reads a variable defined at interfaces between grid levels (column
    /// faces) with the given name from the NetCDF file.
    /// @param var_name name of the variable to read from the file.
    /// @param time_index the time index associated with the value.
    /// @param data a vector into which the data is read. The vector is resized
    ///             as necessary.
    void read_interface_var(const std::string& var_name,
                            size_t time_index,
                            std::vector<Real>& data);

  private:
    // `netcdf.h` defines numerous error codes (integers); this function
    // decodes this integer and throws an exeption with (hopefully) a more
    // helpful message.
    // @throws
    void handle_errcode(const int& ec) const;

    // filename for .nc data
    std::string fname;

    // NetCDF file ID.
    int ncid;
    // Level dimension ID.
    int level_dimid;
    // Interface dimension ID.
    int interface_dimid;
    // Mode dimension ID.
    int mode_dimid;
    // Time dimension ID.
    int time_dimid;
};

} // namespace driver
} // namespace haero
#endif
