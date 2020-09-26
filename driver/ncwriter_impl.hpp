#ifndef HAERO_NC_WRITER_IMPL_HPP
#define HAERO_NC_WRITER_IMPL_HPP

#include "ncwriter.hpp"
#include "ekat/ekat_assert.hpp"

namespace haero {

#define CHECK_NCERR(ec) if (ec != NC_NOERR) handle_errcode(ec, __FILE__, __FUNCTION__, __LINE__)

template <typename ViewType>
void NcWriter::add_var_atts(const int varid, const ekat::units::Units& units, const ViewType& view) {

  EKAT_ASSERT(varid != NC_EBADID);

  const std::string label_str = view.label();
  int retval = nc_put_att_text(ncid, varid, "view_label", label_str.size(), label_str.c_str());
  CHECK_NCERR(retval);
  const std::string unit_str = ekat::units::to_string(units);
  retval = nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
}

template <typename ViewType> VIEW_REAL_TYPE_IS_SP
NcWriter::define_level_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {

  static_assert(ViewType::Rank == 2, "non-aerosol variables should have views with rank = 2");

  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_ASSERT(time_dimid != NC_EBADID && col_dimid != NC_EBADID && level_dimid != NC_EBADID);

  EKAT_REQUIRE(view.extent(0) == ncol() && view.extent(1) == nlev());

  int varid = NC_EBADID;
  const int m_ndims = 3;
  const int dimids[3] = {time_dimid, col_dimid, level_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_FLOAT, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view);
  name_varid_map.emplace(name, varid);
}


template <typename ViewType> VIEW_REAL_TYPE_IS_DP
NcWriter::define_level_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {

  static_assert(ViewType::Rank == 2, "non-aerosol variables should have views with rank = 2");

  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_ASSERT(time_dimid != NC_EBADID && col_dimid != NC_EBADID && level_dimid != NC_EBADID);

  EKAT_REQUIRE(view.extent(0) == ncol() && view.extent(1) == nlev());

  int varid = NC_EBADID;
  const int m_ndims = 3;
  const int dimids[3] = {time_dimid, col_dimid, level_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_DOUBLE, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view);
  name_varid_map.emplace(name, varid);
}

template <typename ViewType> VIEW_REAL_TYPE_IS_SP
NcWriter::define_interface_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {

  static_assert(ViewType::Rank == 2, "non-aerosol variables should have views with rank = 2");
  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_ASSERT(time_dimid != NC_EBADID && col_dimid != NC_EBADID && interface_dimid != NC_EBADID);

  EKAT_REQUIRE(view.extent(0) == ncol() && view.extent(1) == ninterfaces());

  int varid = NC_EBADID;
  const int m_ndims = 3;
  const int dimids[3] = {time_dimid, col_dimid, interface_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_FLOAT, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view);
  name_varid_map.emplace(name, varid);
}

template <typename ViewType> VIEW_REAL_TYPE_IS_DP
NcWriter::define_interface_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {

  static_assert(ViewType::Rank == 2, "non-aerosol variables should have views with rank = 2");
  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_ASSERT(time_dimid != NC_EBADID && col_dimid != NC_EBADID && interface_dimid != NC_EBADID);

  EKAT_REQUIRE(view.extent(0) == ncol() && view.extent(1) == ninterfaces());

  int varid = NC_EBADID;
  const int m_ndims = 3;
  const int dimids[3] = {time_dimid, col_dimid, interface_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_DOUBLE, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view);
  name_varid_map.emplace(name, varid);
}

template <typename ViewType> VIEW_REAL_TYPE_IS_SP
NcWriter::define_modal_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {
  static_assert(ViewType::Rank == 3, "aerosol variables should have views with rank = 3");
  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_ASSERT(time_dimid != NC_EBADID && col_dimid != NC_EBADID &&
              mode_dimid != NC_EBADID && level_dimid != NC_EBADID);

  EKAT_REQUIRE(view.extent(0) == ncol() && view.extent(1) == nmodes() && view.extent(2) == nlev());

  int varid = NC_EBADID;
  const int m_ndims = 4;
  const int dimids[4] = {time_dimid, col_dimid, mode_dimid, level_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_FLOAT, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view);
  name_varid_map.emplace(name, varid);
}

template <typename ViewType> VIEW_REAL_TYPE_IS_DP
NcWriter::define_modal_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {
  static_assert(ViewType::Rank == 3, "aerosol variables should have views with rank = 3");
  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_ASSERT(time_dimid != NC_EBADID && col_dimid != NC_EBADID &&
              mode_dimid != NC_EBADID && level_dimid != NC_EBADID);

  EKAT_REQUIRE(view.extent(0) == ncol() && view.extent(1) == nmodes() && view.extent(2) == nlev());

  int varid = NC_EBADID;
  const int m_ndims = 4;
  const int dimids[4] = {time_dimid, col_dimid, mode_dimid, level_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_DOUBLE, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view);
  name_varid_map.emplace(name, varid);
}

template <typename VT> typename std::enable_if<VT::Rank == 2, void>::type
NcWriter::add_variable_impl(const std::string& varname, const size_t& time_index,
  const size_t& col_index, const VT& v) const {
  const int varid = name_varid_map.at(varname);
  int nvardims = 0;
  int retval = nc_inq_varndims(ncid, varid, &nvardims);
  CHECK_NCERR(retval);

  EKAT_REQUIRE_MSG(nvardims == 3,
    "add_variable_data (rank 2) error: unexpected number of dimensions.");

  int dimids[3] = {NC_EBADID, NC_EBADID, NC_EBADID};
  retval = nc_inq_vardimid(ncid, varid, &dimids[0]);
  CHECK_NCERR(retval);

  EKAT_REQUIRE(dimids[2] == level_dimid || dimids[2] == interface_dimid);

  int array_length = (dimids[2] == level_dimid ? nlev() : ninterfaces());
  const auto viewdata = view1d_to_vector(
    Kokkos::subview(v, col_index, Kokkos::ALL), array_length);
  if (array_length == nlev()) {
    add_level_variable_data(varname, time_index, col_index, viewdata);
  }
  else {
    add_interface_variable_data(varname, time_index, col_index, viewdata);
  }
}

template <typename VT> typename std::enable_if<VT::Rank == 3, void>::type
NcWriter::add_variable_impl(const std::string& varname, const size_t& time_index,
  const size_t& col_index, const VT& v) const {
  const int varid = name_varid_map.at(varname);
  int nvardims = 0;
  int retval = nc_inq_varndims(ncid, varid, &nvardims);
  CHECK_NCERR(retval);

  EKAT_REQUIRE_MSG(nvardims == 4,
    "add_variable_data (rank 3) error: unexpected number of dimensions.");

  int dimids[4] = {NC_EBADID, NC_EBADID, NC_EBADID, NC_EBADID};
  retval = nc_inq_vardimid(ncid, varid, &dimids[0]);
  CHECK_NCERR(retval);

  EKAT_REQUIRE(dimids[2] == mode_dimid && dimids[3] == level_dimid);
  const int array_length = nlev();
  for (size_t i=0; i<nmodes(); ++i) {
    const auto viewdata = view1d_to_vector(
      Kokkos::subview(v, col_index, i, Kokkos::ALL), array_length);
    const size_t start[4] = {time_index, col_index, i, 0};
    const size_t count[4] = {1, 1, 1, static_cast<size_t>(array_length)};
#ifdef HAERO_DOUBLE_PRECISION
    retval = nc_put_vara_double(ncid, varid, start, count, &viewdata[0]);
#else
    retval = nc_put_vara_float(ncid, varid, start, count, &viewdata[0]);
#endif
    CHECK_NCERR(retval);
  }
}

} // namespace haero
#endif
