#ifndef HAERO_NCFILE_HPP
#define HAERO_NCFILE_HPP

#include "haero/haero_config.hpp"
#include "column_base.hpp"
#include "ekat/util/ekat_units.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/ekat_pack.hpp"
#include "netcdf.h"
#include <string>
#include <map>
#include <vector>

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

      @param [in] new_filename name of NetCDF data file to create.
    */
    NcFile(const std::string& new_filename) : fname(new_filename),
     ncid(NC_EBADID), level_dimid(NC_EBADID), interface_dimid(NC_EBADID),
     mode_dimid(NC_EBADID), time_dimid(NC_EBADID), ndims(0), view_var_map() {
      open();
      add_time_dim();
    }

    /** @brief String with basic info about this NcFile.

      @param [in] tab_level indent level for console output.
    */
    std::string info_string(const int& tab_level=0) const;

    /** @brief Adds 2 dimensions to an existing NcFile; one for level midpoints,
      one for level interfaces.

      @throws

      This is required before any column variables can be defined.
      @param [in] nlev number of vertical levels (midpoints) in a column.
    */
    void add_level_dims(const int& nlev);


    /** @brief Adds a dimension for the number of modes.

      @throws

      This is required before any variables indexed by mode can be defined.
      @param [in] nmodes number of modes in a HAERO aerosol computation.
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
    inline int get_nvars() const {return view_var_map.size();}

    /** @brief return the variable id for a view
    */
    template <typename ViewType>
    inline int get_varid(const ViewType& view) const {return view_var_map.at(view.label());}

    /** @brief returns the labels of the views associated with netcdf variables

    */
    std::vector<std::string> get_variable_view_labels() const;

    /** @brief defines a level midpoint variable (single precision real)

      @throws

      @param [in] name (name of variable to write to file)
      @param [in] units
      @param [in] view
    */
    template <typename ViewType=ColumnBase::view_1d>
    typename std::enable_if<std::is_same<typename ekat::ScalarTraits<typename ViewType::value_type>::scalar_type, float>::value,void>::type
    define_level_var(const std::string& name, const ekat::units::Units& units, const ViewType& view);


    /** @brief defines a level midpoint variable (single precision real)

      @throws

      @param [in] name (name of variable to write to file)
      @param [in] units
      @param [in] view
    */
    template <typename ViewType=ColumnBase::view_1d>
    typename std::enable_if<std::is_same<typename ekat::ScalarTraits<typename ViewType::value_type>::scalar_type, double>::value,void>::type
    define_level_var(const std::string& name, const ekat::units::Units& units, const ViewType& view);


    /** @brief defines the time variable

      @throws

      @param [in] units
    */
    template <typename RealType>
    typename std::enable_if<std::is_same<RealType,float>::value,void>::type
    define_time_var(const ekat::units::Units& units = ekat::units::s);

    /** @brief defines the time variable

      @throws

      @param [in] units
    */
    template <typename RealType>
    typename std::enable_if<std::is_same<RealType,double>::value,void>::type
    define_time_var(const ekat::units::Units& units = ekat::units::s);

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

    /** @brief adds metadata attributes to newly defined variables
    */
    template <typename ViewType>
    void add_var_atts(const int varid, const ekat::units::Units& units, const ViewType& view);

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

    std::map<std::string, int> view_var_map;
};

}
#endif
