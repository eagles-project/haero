#ifndef HAERO_NC_WRITER_IMPL_HPP
#define HAERO_NC_WRITER_IMPL_HPP

#include "ncwriter.hpp"
#include "ekat/ekat_assert.hpp"

namespace haero {

#define CHECK_NCERR(ec) if (ec != NC_NOERR) handle_errcode(ec, __FILE__, __FUNCTION__, __LINE__)

template <int Rank, bool IsSimd>
struct NcWriterImpl {
  template <typename VT>
  static void add_variable_data(const int ncid, const int varid, const int col_dimid,
    const int mode_dimid, const int nmodes,
    const int level_dimid, const int interface_dimid, const int nlev,
    const std::string& varname, const size_t time_idx, const size_t col_idx, const VT& view) {}
};

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

template <typename ViewType>
void NcWriter::define_level_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {

  static_assert(ViewType::Rank == 2, "non-aerosol variables should have views with rank = 2");

  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_ASSERT(time_dimid != NC_EBADID && col_dimid != NC_EBADID && level_dimid != NC_EBADID);

  if (ekat::ScalarTraits<typename ViewType::value_type>::is_simd) {
    EKAT_REQUIRE(view.extent(0) == num_columns() && view.extent(1) == pack_info::num_packs(num_levels()));
  }
  else {
    EKAT_REQUIRE(view.extent(0) == num_columns() && view.extent(1) == num_levels());
  }

  int varid = NC_EBADID;
  const int m_ndims = 3;
  const int dimids[3] = {time_dimid, col_dimid, level_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_REAL_KIND, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view);
  name_varid_map.emplace(name, varid);
}

template <typename ViewType>
void NcWriter::define_interface_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {

  static_assert(ViewType::Rank == 2, "non-aerosol variables should have views with rank = 2");
  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_ASSERT(time_dimid != NC_EBADID && col_dimid != NC_EBADID && interface_dimid != NC_EBADID);

  if (ekat::ScalarTraits<typename ViewType::value_type>::is_simd) {
    EKAT_REQUIRE(view.extent(0) == num_columns() && view.extent(1) == pack_info::num_packs(num_interfaces()));
  }
  else {
    EKAT_REQUIRE(view.extent(0) == num_columns() && view.extent(1) == num_interfaces());
  }

  int varid = NC_EBADID;
  const int m_ndims = 3;
  const int dimids[3] = {time_dimid, col_dimid, interface_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_REAL_KIND, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view);
  name_varid_map.emplace(name, varid);
}

template <typename ViewType>
void NcWriter::define_modal_var(const std::string& name, const ekat::units::Units& units, const ViewType& view) {
  static_assert(ViewType::Rank == 3, "aerosol variables should have views with rank = 3");
  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_ASSERT(time_dimid != NC_EBADID && col_dimid != NC_EBADID &&
              mode_dimid != NC_EBADID && level_dimid != NC_EBADID);

  if (ekat::ScalarTraits<typename ViewType::value_type>::is_simd) {
    EKAT_REQUIRE(view.extent(0) == num_columns() && view.extent(1) == num_modes() && view.extent(2) == pack_info::num_packs(num_levels()));
  }
  else {
    EKAT_REQUIRE(view.extent(0) == num_columns() && view.extent(1) == num_modes() && view.extent(2) == num_levels());
  }

  int varid = NC_EBADID;
  const int m_ndims = 4;
  const int dimids[4] = {time_dimid, col_dimid, mode_dimid, level_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_REAL_KIND, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view);
  name_varid_map.emplace(name, varid);
}

template <typename ViewType>
void NcWriter::add_variable_data(const std::string& varname, const size_t& time_idx,
  const size_t& col_idx, const ViewType& view) const {

  EKAT_REQUIRE(time_idx < num_timesteps());
  EKAT_REQUIRE(col_idx < num_columns());

  return NcWriterImpl<ViewType::Rank,
    ekat::ScalarTraits<typename ViewType::value_type>::is_simd>::add_variable_data(
    ncid, name_varid_map.at(varname), col_dimid, mode_dimid, num_modes(),
    level_dimid, interface_dimid, num_levels(), varname, time_idx, col_idx, view);
}

template <>
struct NcWriterImpl<2, false> {
  template <typename VT>
  static void add_variable_data(const int ncid, const int varid, const int col_dimid,
    const int mode_dimid, const int nmodes,
    const int level_dimid, const int interface_dimid, const int nlev,
    const std::string& varname, const size_t time_idx, const size_t col_idx, const VT& view) {

    static_assert(VT::Rank == 2, "non-aerosol views assumed to have rank 2");

    int nvardims = 0;
    int retval = nc_inq_varndims(ncid, varid, &nvardims);
    EKAT_REQUIRE(nvardims == 3 && retval == NC_NOERR);
    EKAT_REQUIRE(view.extent(1) == nlev || view.extent(1) == nlev+1);

    int dimids[3] = {NC_EBADID, NC_EBADID, NC_EBADID};
    retval = nc_inq_vardimid(ncid, varid, &dimids[0]);
    EKAT_REQUIRE(retval == NC_NOERR);
    EKAT_REQUIRE(dimids[2] == level_dimid || dimids[2] == interface_dimid);

    const auto sv = Kokkos::subview(view, col_idx, Kokkos::ALL);
    auto hsv = Kokkos::create_mirror_view(sv);
    Kokkos::deep_copy(hsv, sv);

    const size_t array_length = (dimids[2] == level_dimid ? nlev : nlev+1);
    for (size_t i=0; i<array_length; ++i) {
      const size_t idx[3] = {time_idx, col_idx, i};
#ifdef HAERO_DOUBLE_PRECISION
      retval = nc_put_var1_double(ncid, varid, &idx[0], &hsv(i));
#else
      retval = nc_put_var1_float (ncid, varid, &idx[0], &hsv(i));
#endif
      EKAT_ASSERT(retval == NC_NOERR);
    }
  }
};

template <>
struct NcWriterImpl<2, true> {
  template <typename VT>
  static void add_variable_data(const int ncid, const int varid, const int col_dimid,
    const int mode_dimid, const int nmodes,
    const int level_dimid, const int interface_dimid, const int nlev,
    const std::string& varname, const size_t time_idx, const size_t col_idx, const VT& view) {

    static_assert(VT::Rank == 2, "non-aerosol views assumed to have rank 2");

    int nvardims = 0;
    int retval = nc_inq_varndims(ncid, varid, &nvardims);
    EKAT_REQUIRE(nvardims == 3 && retval == NC_NOERR);
    EKAT_REQUIRE(view.extent(1) == pack_info::num_packs(nlev) ||
                 view.extent(1) == pack_info::num_packs(nlev+1));

    int dimids[3] = {NC_EBADID, NC_EBADID, NC_EBADID};
    retval = nc_inq_vardimid(ncid, varid, &dimids[0]);
    EKAT_REQUIRE(retval == NC_NOERR);
    EKAT_REQUIRE(dimids[2] == level_dimid || dimids[2] == interface_dimid);

    const auto sv = Kokkos::subview(view, col_idx, Kokkos::ALL);
    auto hsv = Kokkos::create_mirror_view(sv);
    Kokkos::deep_copy(hsv, sv);

    const size_t array_length = (dimids[2] == level_dimid ? nlev : nlev+1);
    for (size_t pack_idx = 0; pack_idx < pack_info::num_packs(array_length); ++pack_idx) {
      for (size_t vec_idx = 0; vec_idx < pack_info::vec_end(array_length, pack_idx); ++vec_idx) {
        const size_t idx[3] = {time_idx, col_idx,
          static_cast<size_t>(pack_info::array_idx(pack_idx, vec_idx))};
        const Real val = hsv(pack_idx)[vec_idx];
#ifdef HAERO_DOUBLE_PRECISION
        retval = nc_put_var1_double(ncid, varid, &idx[0], &val);
#else
        retval = nc_put_var1_float (ncid, varid, &idx[0], &val);
#endif
        EKAT_ASSERT(retval == NC_NOERR);
      }
    }
  }
};

template <>
struct NcWriterImpl<3,false> {
  template <typename VT>
  static void add_variable_data(const int ncid, const int varid, const int col_dimid,
    const int mode_dimid, const int nmodes,
    const int level_dimid, const int interface_dimid, const int nlev,
    const std::string& varname, const size_t time_idx, const size_t col_idx, const VT& view) {

    static_assert(VT::Rank == 3, "aerosol data assumed to have rank 3 views");

    int nvardims = 0;
    int retval = nc_inq_varndims(ncid, varid, &nvardims);
    EKAT_REQUIRE(nvardims == 4 && retval == NC_NOERR);

    int dimids[4] = {NC_EBADID, NC_EBADID, NC_EBADID, NC_EBADID};
    retval = nc_inq_vardimid(ncid, varid, &dimids[0]);

    EKAT_REQUIRE(dimids[2] == mode_dimid && dimids[3] == level_dimid);
    EKAT_REQUIRE(view.extent(1) == nmodes);
    EKAT_REQUIRE(view.extent(2) == nlev);

    const auto sv = Kokkos::subview(view, col_idx, Kokkos::ALL, Kokkos::ALL);
    auto hsv = Kokkos::create_mirror_view(sv);
    Kokkos::deep_copy(hsv, sv);

    for (size_t mode_idx = 0; mode_idx < nmodes; ++mode_idx) {
      for (size_t lev_idx = 0; lev_idx < nlev; ++lev_idx) {
        const size_t idx[4] = {time_idx, col_idx, mode_idx, lev_idx};
#ifdef HAERO_DOUBLE_PRECISION
        retval = nc_put_var1_double(ncid, varid, &idx[0], &hsv(mode_idx, lev_idx));
#else
        retval = nc_put_var1_float(ncid, varid, &idx[0], &hsv(mode_idx, lev_idx));
#endif
      }
      EKAT_ASSERT(retval == NC_NOERR);
    }
  }
};

template <>
struct NcWriterImpl<3,true> {
  template <typename VT>
  static void add_variable_data(const int ncid, const int varid, const int col_dimid,
    const int mode_dimid, const int nmodes,
    const int level_dimid, const int interface_dimid, const int nlev,
    const std::string& varname, const size_t time_idx, const size_t col_idx, const VT& view) {

    static_assert(VT::Rank == 3, "aerosol data assumed to have rank 3 views");
    EKAT_REQUIRE(view.extent(1) == nmodes);
    EKAT_REQUIRE(view.extent(2) == pack_info::num_packs(nlev));

    int nvardims = 0;
    int retval = nc_inq_varndims(ncid, varid, &nvardims);
    EKAT_REQUIRE(nvardims == 4);


    int dimids[4] = {NC_EBADID, NC_EBADID, NC_EBADID, NC_EBADID};
    retval = nc_inq_vardimid(ncid, varid, &dimids[0]);

    EKAT_REQUIRE(dimids[2] == mode_dimid && dimids[3] == level_dimid);

    const auto sv = Kokkos::subview(view, col_idx, Kokkos::ALL, Kokkos::ALL);
    auto hsv = Kokkos::create_mirror_view(sv);
    Kokkos::deep_copy(hsv, sv);

    const int array_length = nlev;

    for (size_t mode_idx=0; mode_idx<sv.extent(0); ++mode_idx) {
      for (int pack_idx = 0; pack_idx < sv.extent(1); ++pack_idx) {
        for (int vec_idx = 0; vec_idx < pack_info::vec_end(array_length, pack_idx); ++vec_idx) {
          const size_t idx[4] = {time_idx, col_idx,
            mode_idx, static_cast<size_t>(pack_info::array_idx(pack_idx, vec_idx))};
          const Real val = hsv(mode_idx, pack_idx)[vec_idx];
#ifdef HAERO_DOUBLE_PRECISION
          retval = nc_put_var1_double(ncid, varid, &idx[0], &val);
#else
          retval = nc_put_var1_float (ncid, varid, &idx[0], &val);
#endif
          EKAT_ASSERT(retval == NC_NOERR);
        }
      }
    }

  }
};

} // namespace haero
#endif
