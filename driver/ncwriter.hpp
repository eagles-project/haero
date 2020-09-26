#ifndef HAERO_NC_WRITER_HPP
#define HAERO_NC_WRITER_HPP

#include "haero/haero_config.hpp"
#include "haero/view_pack_helpers.hpp"
#include "haero/utils.hpp"
#include "column_base.hpp"
#include "ekat/util/ekat_units.hpp"
#include "ekat/ekat_assert.hpp"
#include "netcdf.h"
#include <string>
#include <map>
#include <vector>

namespace haero {

class NcWriter {
  public:
    /// key-value pairs for metadata attributes
    typedef std::pair<std::string,std::string> text_att_type;

    /// Use the correct netcdf real kind parameter
    static constexpr int NC_REAL_KIND = (std::is_same<Real,double>::value ? NC_DOUBLE : NC_FLOAT);

    /** @brief Constructor.  Must be called from host.

      Assumes at most 1 file open per NcWriter instance.
      All values initialized to invalid data; they will be set to valid data by other
      member functions.

      Data model:
        - Non-aerosol data have rank = 3, with dimensions (time, column, level)
        - Aerosol data have rank = 4, with dimensions (time, column, mode, level)

      Writes HAERO metadata (version, git revision) to file.

      @throws

      @param [in] new_filename name of NetCDF data file to create.
    */
    NcWriter(const std::string& new_filename) : fname(new_filename),
     ncid(NC_EBADID), col_dimid(NC_EBADID), level_dimid(NC_EBADID), interface_dimid(NC_EBADID),
     mode_dimid(NC_EBADID), time_dimid(NC_EBADID), ndims(0), name_varid_map() {
      const auto fext = get_filename_ext(new_filename);
      EKAT_REQUIRE_MSG(fext == ".nc",
        "NcWriter error: filename extension (" + fext + ")  must be .nc");
      open();
      add_time_dim();
    }

    /** @brief adds a metadata attribute to the file

      @param [in] att_pair
    */
    void add_file_attribute(const text_att_type& att_pair) const;

    /** @brief String with basic info about this NcWriter.

      @param [in] tab_level indent level for console output.
    */
    std::string info_string(const int& tab_level=0) const;

    /** @brief Adds a dimension for columns

    */
    void add_column_dim(const int& ncol=1);

    int ncol() const;
    int nlev() const;
    inline int ninterfaces() const {return nlev() + 1;}
    int nmodes() const;
    int ntimesteps() const;

    /** @brief Adds 2 dimensions to an existing NcWriter; one for level midpoints,
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
      @param [in] mode names (e.g., Aitken, accumulation, etc.)
      @param [in] mode sizes in min/max pairs
    */
    void add_mode_dim(const int& nmodes, const std::vector<std::string>& mode_names,
      const std::vector<std::pair<Real,Real>>& mode_max_mins);


    /** @brief Close the data file associated with this NcWriter.

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
    inline int get_nvars() const {return name_varid_map.size();}

    /** @brief return the variable id associated with a name

      @throws

      @param [in] name
      @return varid
    */
    inline int get_varid(const std::string& name) const {return name_varid_map.at(name);}

    /** @brief returns the names of the variables currently defined
    */
    std::vector<std::string> get_variable_names() const;

    /** @brief defines a level midpoint variable

      @param [in] name
      @param [in] units
      @param [in] atts
    */
    void define_level_var(const std::string& name, const ekat::units::Units& units,
      const std::vector<text_att_type>& atts=std::vector<text_att_type>());

    /** @brief defines a level midpoint variable

      @param [in] name
      @param [in] units
      @param [in] atts
    */
    void define_interface_var(const std::string& name, const ekat::units::Units& units,
      const std::vector<text_att_type>& atts=std::vector<text_att_type>());

    /** @brief defines a modal aerosol variable on level midpoints.

      @param [in] name
      @param [in] units
      @param [in] atts
    */
    void define_modal_var(const std::string& name, const ekat::units::Units& units,
      const std::vector<text_att_type>& atts=std::vector<text_att_type>());

    /** @brief defines a level midpoint variable from a view of that variable

      Here, you're allowed to define a variable name that is different than the view label.
      The view label will be recorded as an attribute of the variable.

      @throws

      @param [in] name (name of variable to write to file)
      @param [in] units
      @param [in] view
    */
    template <typename ViewType>
    void define_level_var(const std::string& name, const ekat::units::Units& units, const ViewType& view);

    /** @brief defines a modal aerosol variable on level midpoints.

      @throws

      @param [in] name
      @param [in] units
      @param [in] view
    */
    template <typename ViewType>
    void define_modal_var(const std::string& name, const ekat::units::Units& units, const ViewType& view);

    /** @brief defines an interface variable  from a view of that variable

      Here, you're allowed to define a variable name that is different than the view label.
      The view label will be recorded as an attribute of the variable.

      @throws

      @param [in] name (name of variable to write to file)
      @param [in] units
      @param [in] view
    */
    template <typename ViewType>
    void define_interface_var(const std::string& name, const ekat::units::Units& units, const ViewType& view);

    /** @brief defines the time variable

      @throws

      @param [in] units (default = seconds)
      @param [in] t0 first value of time coordinate (default = 0)
    */
    void define_time_var(const ekat::units::Units& units=ekat::units::s, const Real t0=0.0);

    /** @brief adds a time value to the time variable

    @warning This will increase the size of the time dimension; re-index your variables accordingly!

      @param [in] t
    */
    void add_time_value(const Real& t) const;

    /** @brief write one column's data from a std::vector to the nc file.

      @throws

      @param [in] varname name of variable to write
      @param [in] time_index time index associated with data
      @param [in] data
    */
    void add_level_variable_data(const std::string& varname, const size_t& time_index,
      const size_t& col_index, const std::vector<Real>& data) const;

    /** @brief write one column's variable data from a std::vector to the nc file.

      @throws

      @param [in] varname name of variable to write
      @param [in] time_index time index associated with data
      @param [in] data
    */
    void add_interface_variable_data(const std::string& varname, const size_t& time_index,
      const size_t& col_index, const std::vector<Real>& data) const;

    template <typename ViewType>
    void add_variable_data(const std::string& varname, const size_t& time_index,
      const size_t& col_index, const ViewType& v) const;

  protected:

    /** @brief `netcdf.h` defines numerous error codes (integers); this function
      decodes this integers and throws an exeption with (hopefully) a more helpful message.

      @throws

      @todo improve this function to include more info & use ekat_assert.hpp macros
    */
    void handle_errcode(const int& ec, const std::string& file="", const std::string& fn="", const int& line=-1) const;

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
    /// Column dimension ID.  Assigned by add_column_dim()
    int col_dimid;
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
    /// key = variable name, value = variable id
    std::map<std::string, int> name_varid_map;

};

}
#endif
