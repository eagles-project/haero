#include "ncreader.hpp"

#include <cstring>
#include <ekat/ekat_assert.hpp>
#include <exception>
#include <iostream>
#include <sstream>

#include "haero/haero.hpp"
#include "haero/utils.hpp"

namespace haero {
namespace driver {

#define CHECK_ERR(ec) \
  if (ec != NC_NOERR) handle_errcode(ec)

NcReader::NcReader(const std::string& filename)
    : fname(filename),
      ncid(NC_EBADID),
      level_dimid(NC_EBADID),
      interface_dimid(NC_EBADID),
      mode_dimid(NC_EBADID),
      time_dimid(NC_EBADID) {
  // We open the file and read it entirely into memory.
  int retval = nc_open(fname.c_str(), NC_NETCDF4 | NC_DISKLESS, &ncid);
  CHECK_ERR(retval);

  // Make sure it's actually a Haero file!
  const char* haero_str = "High-performance AEROsols standalone driver";
  char text[128];
  retval = nc_get_att_text(ncid, NC_GLOBAL, "HAERO", text);
  CHECK_ERR(retval);
  if (std::strcmp(text, haero_str) != 0) {
    throw std::runtime_error("NetCDF file is not a Haero file!");
  }

  const char* version_string = version();
  retval = nc_get_att_text(ncid, NC_GLOBAL, "HAERO_version", text);
  CHECK_ERR(retval);
  if (std::strcmp(text, version_string) != 0) {
    throw std::runtime_error(
        "NetCDF Haero file version does not match haero library version!");
  }

  const char* revision_str = revision();
  retval = nc_get_att_text(ncid, NC_GLOBAL, "HAERO_revision", text);
  if (std::strcmp(text, revision_str) != 0) {
    throw std::runtime_error(
        "NetCDF Haero file git revision does not match haero library "
        "revision!");
  }

  // Read some dimension IDs.
  retval = nc_inq_dimid(ncid, "level_midpts", &level_dimid);
  CHECK_ERR(retval);
  retval = nc_inq_dimid(ncid, "level_interfaces", &interface_dimid);
  CHECK_ERR(retval);
  retval = nc_inq_dimid(ncid, "modes", &mode_dimid);
  CHECK_ERR(retval);
  retval = nc_inq_dimid(ncid, "time", &time_dimid);
  CHECK_ERR(retval);
}

NcReader::~NcReader() { nc_close(ncid); }

void NcReader::read_level_var(const std::string& var_name, size_t time_index,
                              std::vector<Real>& data) {}

void NcReader::read_interface_var(const std::string& var_name,
                                  size_t time_index, std::vector<Real>& data) {}

void NcReader::handle_errcode(const int& ec) const {
  std::ostringstream ss;
  ss << "NcReader error: ";
  switch (ec) {
    case (NC_NOERR): {
      // no error: should not have called this routine
      assert(ec != NC_NOERR);
      return;
    }
    case (NC_EEXIST): {
      ss << "File exists, and overwrite is not permitted.";
      break;
    }
    case (NC_EPERM): {
      ss << "cannot create file; a permission issue, or trying to write to "
            "read-only file.";
      break;
    }
    case (NC_ENOMEM): {
      ss << "out of memory.";
      break;
    }
    case (NC_ENFILE): {
      ss << "too many netcdf files open.";
      break;
    }
    case (NC_EHDFERR): {
      ss << "unknown hdf5 error.";
      break;
    }
    case (NC_EFILEMETA): {
      ss << "metadata write error.";
      break;
    }
    case (NC_EDISKLESS): {
      ss << "error creating file in memory.";
      break;
    }
    case (NC_EINVAL): {
      ss << "more than one fill value defined, or trying to set global "
            "_FillValue, or invalid input.";
      break;
    }
    case (NC_ENOTVAR): {
      ss << "could not locate variable id.";
      break;
    }
    case (NC_EBADTYPE): {
      ss << "fill value and var must have same type.";
      break;
    }
    case (NC_ELATEFILL): {
      ss << "Fill values must be written while file is still in 'define' mode.";
      break;
    }
    case (NC_EMAXNAME): {
      ss << "name is too long.";
      break;
    }
    case (NC_EDIMSIZE): {
      ss << "invalid dim size.";
      break;
    }
    case (NC_ENOTINDEFINE): {
      ss << "netcdf not in define mode.";
      break;
    }
    case (NC_EUNLIMIT): {
      ss << "NC_UNLIMITED is already used.";
      break;
    }
    case (NC_EMAXDIMS): {
      ss << "ndims > NC_MAX_DIMS";
      break;
    }
    case (NC_ENAMEINUSE): {
      ss << "name already used.";
      break;
    }
    case (NC_EBADNAME): {
      ss << "name breaks NetCDF naming rules.";
      break;
    }
    default: {
      ss << "unknown netcdf error; ec = " << ec;
    }
  }
  throw std::runtime_error(ss.str());
}

}  // namespace driver
}  // namespace haero
