#ifndef HAERO_NC_WRITER_IMPL_HPP
#define HAERO_NC_WRITER_IMPL_HPP

#include "ncwriter.hpp"
#include <cassert>

namespace haero {

#define CHECK_NCERR(ec) if (ec != NC_NOERR) handle_errcode(ec)

template <typename ViewType>
void NcWriter::add_var_atts(const int varid, const ekat::units::Units& units, const ViewType& view) {

  assert(varid != NC_EBADID);

  const std::string label_str = view.label();
  int retval = nc_put_att_text(ncid, varid, "view_label", label_str.size(), label_str.c_str());
  CHECK_NCERR(retval);
  const std::string unit_str = ekat::units::to_string(units);
  retval = nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
}

template <typename ViewType> VIEW_REAL_TYPE_IS_SP
NcWriter::define_level_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {

   assert(time_dimid != NC_EBADID && level_dimid != NC_EBADID);

  int varid = NC_EBADID;
  const int m_ndims = 2;
  const int dimids[2] = {time_dimid, level_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_FLOAT, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view);
  view_var_map.emplace(view.label(), varid);
}


template <typename ViewType> VIEW_REAL_TYPE_IS_DP
NcWriter::define_level_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {

  assert(time_dimid != NC_EBADID && level_dimid != NC_EBADID);

  int varid = NC_EBADID;
  const int m_ndims = 2;
  const int dimids[2] = {time_dimid, level_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_DOUBLE, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view);
  view_var_map.emplace(view.label(), varid);
}

template <typename RealType>
typename std::enable_if<std::is_same<RealType,float>::value,void>::type
NcWriter::define_time_var(const ekat::units::Units& units) {

  assert(time_dimid != NC_EBADID);

  int varid = NC_EBADID;
  const int m_ndims = 1;
  const int dimids[1] = {time_dimid};
  int retval = nc_def_var(ncid, "time", NC_FLOAT, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  const auto unit_str = ekat::units::to_string(units);
  retval = nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
  view_var_map.emplace("time", varid);
}

template <typename RealType>
typename std::enable_if<std::is_same<RealType,double>::value,void>::type
NcWriter::define_time_var(const ekat::units::Units& units) {

  assert(time_dimid != NC_EBADID);

  int varid = NC_EBADID;
  const int m_ndims = 1;
  const int dimids[1] = {time_dimid};
  int retval = nc_def_var(ncid, "time", NC_DOUBLE, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  const auto unit_str = ekat::units::to_string(units);
  retval = nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
  view_var_map.emplace("time", varid);
}

template <typename ViewType> VIEW_REAL_TYPE_IS_SP
NcWriter::define_interface_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {

  assert(time_dimid != NC_EBADID && interface_dimid != NC_EBADID);

  int varid = NC_EBADID;
  const int m_ndims = 2;
  const int dimids[2] = {time_dimid, interface_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_FLOAT, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view);
  view_var_map.emplace(view.label(), varid);
}

template <typename ViewType> VIEW_REAL_TYPE_IS_DP
NcWriter::define_interface_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {

  assert(time_dimid != NC_EBADID && interface_dimid != NC_EBADID);

  int varid = NC_EBADID;
  const int m_ndims = 2;
  const int dimids[2] = {time_dimid, interface_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_DOUBLE, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view);
  view_var_map.emplace(view.label(), varid);
}

// template <typename ViewType=ColumnBase::view_2d> VIEW_REAL_TYPE_IS_SP
// NcWriter::define_modal_level_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {
//
//   assert(time_dimid != NC_EBADID && level_dimid != NC_EBADID && mode_dimid != NC_EBADID);
//
//   int varid = NC_EBADID;
//   const int m_ndims=3;
//   const int dimids[3] = {time_dimid, mode_dimid, level_dimid};
//   int retval = nc_def_var(ncid, name.c_str(), NC_FLOAT, m_ndims, dimids, &varid);
//   CHECK_NCERR(retval);
//   add_var_atts(varid, units, view);
//   view_var_map.emplace(view.label(), varid);
// }
//
// template <typename ViewType=ColumnBase::view_2d> VIEW_REAL_TYPE_IS_DP
// NcWriter::define_modal_level_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {
//
//   assert(time_dimid != NC_EBADID && level_dimid != NC_EBADID && mode_dimid != NC_EBADID);
//
//   int varid = NC_EBADID;
//   const int m_ndims=3;
//   const int dimids[3] = {time_dimid, mode_dimid, level_dimid};
//   int retval = nc_def_var(ncid, name.c_str(), NC_DOUBLE, m_ndims, dimids, &varid);
//   CHECK_NCERR(retval);
//   add_var_atts(varid, units, view);
//   view_var_map.emplace(view.label(), varid);
// }

}
#endif
