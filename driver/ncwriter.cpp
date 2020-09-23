#include "ncwriter.hpp"
#include "ncwriter_impl.hpp"
#include "haero/haero.hpp"
#include "haero/utils.hpp"
#include <exception>
#include <iostream>
#include <sstream>
#include <cassert>

namespace haero {

void NcWriter::open() {
  int retval = nc_create(fname.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid);
  CHECK_NCERR(retval);
  const std::string haero_str = "High-performance AEROsols standalone driver";
  retval = nc_put_att_text(ncid, NC_GLOBAL, "HAERO", haero_str.size(),
    haero_str.c_str());
  CHECK_NCERR(retval);
  const std::string version_string(version());
  retval = nc_put_att_text(ncid, NC_GLOBAL, "HAERO_version",
    version_string.size(), version_string.c_str());
  CHECK_NCERR(retval);
  const std::string revision_str(revision());
  retval = nc_put_att_text(ncid, NC_GLOBAL, "HAERO_revision",
    revision_str.size(), revision_str.c_str());
}

void NcWriter::add_level_dims(const int& nlev) {
  if (level_dimid == NC_EBADID and interface_dimid == NC_EBADID) {
    int retval = nc_def_dim(ncid, "level_midpts", nlev, &level_dimid);
    CHECK_NCERR(retval);
    retval = nc_def_dim(ncid, "level_interfaces", nlev+1, &interface_dimid);
    ndims += 2;
  }
  else {
    throw std::logic_error("level dimensions already defined.");
  }
}

void NcWriter::add_mode_dim(const int& nmodes) {
  if (mode_dimid == NC_EBADID) {
    int retval = nc_def_dim(ncid, "modes", nmodes, &mode_dimid);
    CHECK_NCERR(retval);
    ndims++;
  }
  else {
    throw std::logic_error("mode dimension already defined.");
  }
}

void NcWriter::add_time_dim() {
  if (time_dimid == NC_EBADID) {
    int retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dimid);
    CHECK_NCERR(retval);
    ndims++;
  }
  else {
    throw std::logic_error("time dimension already defined.");
  }
}

void NcWriter::close() {
  int retval = nc_close(ncid);
  CHECK_NCERR(retval);
}

std::string NcWriter::info_string(const int& tab_level) const {
  std::ostringstream ss;
  auto tabstr = indent_string(tab_level);
  ss << tabstr << "NcWriter info:\n";
  tabstr += "\t";
  ss << tabstr << "fname = " << fname << '\n';
  ss << tabstr << "ncid = " << ncid << '\n';
  ss << tabstr << "level_dimid = " << level_dimid << '\n';
  ss << tabstr << "interface_dimid = " << interface_dimid << '\n';
  ss << tabstr << "mode_dimid = " << mode_dimid << '\n';
  ss << tabstr << "time_dimid = " << time_dimid << '\n';
  ss << tabstr << "ndims = " << ndims << '\n';
  ss << tabstr << "nvars = " << get_nvars() << '\n';
  return ss.str();
}

void NcWriter::define_level_var(const std::string& name, const ekat::units::Units& units,
      const std::vector<text_att_type>& atts) {

  assert(time_dimid != NC_EBADID && level_dimid != NC_EBADID);

  int varid = NC_EBADID;
  const int m_ndims = 2;
  const int dimids[2] = {time_dimid, level_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_REAL_KIND, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  const std::string unit_str = ekat::units::to_string(units);
  retval = nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
  if (!atts.empty()) {
    for (int i=0; i<atts.size(); ++i) {
      retval = nc_put_att_text(ncid, varid, atts[i].first.c_str(),
        atts[i].second.size(), atts[i].second.c_str());
      CHECK_NCERR(retval);
    }
  }
  name_varid_map.emplace(name, varid);
}

void NcWriter::define_time_var(const ekat::units::Units& units) {

  assert(time_dimid != NC_EBADID);

  int varid = NC_EBADID;
  int retval = nc_def_var(ncid, "time", NC_REAL_KIND, 1, &time_dimid, &varid);
  CHECK_NCERR(retval);
  const std::string unit_str = ekat::units::to_string(units);
  retval  = nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
  name_varid_map.emplace("time", varid);
}

void NcWriter::define_interface_var(const std::string& name, const ekat::units::Units& units,
      const std::vector<text_att_type>& atts) {

  assert(time_dimid != NC_EBADID && interface_dimid != NC_EBADID);

  int varid = NC_EBADID;
  const int m_ndims = 2;
  const int dimids[2] = {time_dimid, interface_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_REAL_KIND, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  const std::string unit_str = ekat::units::to_string(units);
  retval = nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
  if (!atts.empty()) {
    for (int i=0; i<atts.size(); ++i) {
      retval = nc_put_att_text(ncid, varid, atts[i].first.c_str(),
        atts[i].second.size(), atts[i].second.c_str());
      CHECK_NCERR(retval);
    }
  }
  name_varid_map.emplace(name, varid);
}


std::vector<std::string> NcWriter::get_variable_names() const {
  std::vector<std::string> result;
  for (const auto& v : name_varid_map) {
    result.push_back(v.first);
  }
  return result;
}

void NcWriter::handle_errcode(const int& ec) const {
  std::ostringstream ss;
  ss << "NcWriter error: ";
  switch (ec) {
    case (NC_NOERR) : {
      // no error: should not have called this routine
      assert(ec != NC_NOERR);
      return;
    }
    case (NC_EEXIST) : {
      ss << "File exists, and overwrite is not permitted.";
      break;
    }
    case (NC_EPERM) : {
      ss << "cannot create file; a permission issue, or trying to write to read-only file.";
      break;
    }
    case (NC_ENOMEM) : {
      ss << "out of memory.";
      break;
    }
    case (NC_ENFILE) : {
      ss << "too many netcdf files open.";
      break;
    }
    case (NC_EHDFERR) : {
      ss << "unknown hdf5 error.";
      break;
    }
    case (NC_EFILEMETA) : {
      ss << "metadata write error.";
      break;
    }
    case (NC_EDISKLESS) : {
      ss << "error creating file in memory.";
      break;
    }
    case (NC_EINVAL) : {
      ss << "more than one fill value defined, or trying to set global _FillValue, or invalid input.";
      break;
    }
    case (NC_ENOTVAR) : {
      ss << "could not locate variable id.";
      break;
    }
    case (NC_EBADTYPE) : {
      ss << "fill value and var must have same type.";
      break;
    }
    case (NC_ELATEFILL) : {
      ss << "Fill values must be written while file is still in 'define' mode.";
      break;
    }
    case (NC_EMAXNAME) : {
      ss << "name is too long.";
      break;
    }
    case (NC_EDIMSIZE) : {
      ss << "invalid dim size.";
      break;
    }
    case (NC_ENOTINDEFINE) : {
      ss << "netcdf not in define mode.";
      break;
    }
    case (NC_EUNLIMIT) : {
      ss << "NC_UNLIMITED is already used.";
      break;
    }
    case (NC_EMAXDIMS) : {
      ss << "ndims > NC_MAX_DIMS";
      break;
    }
    case (NC_ENAMEINUSE) : {
      ss << "name already used.";
      break;
    }
    case (NC_EBADNAME) : {
      ss << "name breaks NetCDF naming rules.";
      break;
    }
    default : {
      ss << "unknown netcdf error; ec = " << ec;
    }
  }
  throw std::runtime_error(ss.str());
}

} // namespace haero
