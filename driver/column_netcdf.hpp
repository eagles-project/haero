#ifndef HAERO_COLUMN_NETCDF_HPP
#define HAERO_COLUMN_NETCDF_HPP

#include "haero/haero_config.hpp"
// #include "column_base.hpp" // requires column_dynamics branch
// #include "dyn_column.hpp" // requires column_dynamics branch
#include <string>
// #include "ekat/util/ekat_units.hpp" // requires updated ekat
#include "netcdf.h"

namespace haero {

class NetCDFFileHandler {
  public:
    NetCDFFileHandler(const std::string& new_filename) : fname(new_filename),
     ncid(NC_EBADID), level_dimid(NC_EBADID), interface_dimid(NC_EBADID),
     mode_dimid(NC_EBADID), time_dimid(NC_EBADID), ndims(0), nvars(0) {
      open_file();
    }

    std::string info_string(const int& tab_level=0) const;

    void add_level_dims(const int& nlev);

    void add_mode_dim(const int& nmodes);

    void add_time_dim();

    void open_file();

    void close_file();

    inline int get_ndims() const {return ndims;}

    inline std::string get_fname() const {return fname;}

    inline int get_nvars() const {return nvars;}

  protected:
    void handle_errcode(const int& ec) const;

    std::string fname;
    int ncid;
    int level_dimid;
    int interface_dimid;
    int mode_dimid;
    int time_dimid;
    int ndims;
    int nvars;
};

}
#endif
