#ifndef HAERO_NC_WRITER_IMPL_HPP
#define HAERO_NC_WRITER_IMPL_HPP

#include "ncwriter.hpp"
#include "ekat/ekat_assert.hpp"

namespace haero {
namespace driver {

#define CHECK_NCERR(ec) if (ec != NC_NOERR) handle_errcode(ec, __FILE__, __FUNCTION__, __LINE__)

template <int Rank, bool IsSimd>
struct NcWriterImpl {
  template <typename VT>
  static void add_variable_data(const int ncid, const int varid,
    const int mode_dimid, const int nmodes,
    const int level_dimid, const int interface_dimid, const int nlev,
    const int species_dimid, const int nspec,
    const size_t time_idx, const int mode_idx, const int spec_idx, const VT& view) {}
};

template <typename ViewType>
void NcWriter::add_var_atts(const int varid, const ekat::units::Units& units,
  const ViewType& view, const std::vector<text_att_type>& var_atts) {

  EKAT_ASSERT(varid != NC_EBADID);

  const std::string label_str = view.label();
  int retval = nc_put_att_text(ncid, varid, "view_label", label_str.size(), label_str.c_str());
  CHECK_NCERR(retval);
  const std::string unit_str = ekat::units::to_string(units);
  retval = nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
  for (int i=0; i<var_atts.size(); ++i) {
    retval = nc_put_att_text(ncid, varid, var_atts[i].first.c_str(), var_atts[i].second.size(),
      var_atts[i].second.c_str());
    CHECK_NCERR(retval);
  }
}

template <typename ViewType>
void NcWriter::define_level_var(const std::string& name, const ekat::units::Units& units,
  const ViewType& view, const std::vector<text_att_type>& var_atts) {

  static_assert(ViewType::Rank == 1, "non-aerosol variables should have views with rank = 1");

  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_ASSERT(time_dimid != NC_EBADID && level_dimid != NC_EBADID);

  if (ekat::ScalarTraits<typename ViewType::value_type>::is_simd) {
    EKAT_REQUIRE(view.extent(0) == pack_info::num_packs(num_levels()));
  }
  else {
    EKAT_REQUIRE(view.extent(0) == num_levels());
  }

  int varid = NC_EBADID;
  const int m_ndims = 2;
  const int dimids[2] = {time_dimid, level_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_REAL_KIND, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view, var_atts);
  name_varid_map.emplace(name, varid);
}

template <typename ViewType>
void NcWriter::define_interface_var(const std::string& name, const ekat::units::Units& units,
 const ViewType& view, const std::vector<text_att_type>& var_atts) {

  static_assert(ViewType::Rank == 1, "non-aerosol variables should have views with rank = 1");
  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_ASSERT(time_dimid != NC_EBADID && interface_dimid != NC_EBADID);

  if (ekat::ScalarTraits<typename ViewType::value_type>::is_simd) {
    EKAT_REQUIRE(view.extent(0) == pack_info::num_packs(num_interfaces()));
  }
  else {
    EKAT_REQUIRE(view.extent(0) == num_interfaces());
  }

  int varid = NC_EBADID;
  const int m_ndims = 2;
  const int dimids[2] = {time_dimid, interface_dimid};
  int retval = nc_def_var(ncid, name.c_str(), NC_REAL_KIND, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  add_var_atts(varid, units, view, var_atts);
  name_varid_map.emplace(name, varid);
}


template <typename ViewType>
void NcWriter::add_aerosol_variable_data(const std::string& varname, const size_t& time_idx,
  const int mode_idx, const int spec_idx, const ViewType& view) const {

  EKAT_REQUIRE(time_idx < num_timesteps());

  return NcWriterImpl<ViewType::Rank,
    ekat::ScalarTraits<typename ViewType::value_type>::is_simd>::add_variable_data(
    ncid, name_varid_map.at(varname),
    mode_dimid, num_modes(),
    level_dimid, interface_dimid, num_levels(),
    aerosol_dimid, num_aerosols(),
    time_idx, mode_idx, spec_idx, view);
}

template <typename ViewType>
void NcWriter::add_variable_data(const std::string& varname, const size_t& time_idx,
  const int mode_idx, const int spec_idx, const ViewType& view) const {

  EKAT_REQUIRE(time_idx < num_timesteps());

  return NcWriterImpl<ViewType::Rank,
    ekat::ScalarTraits<typename ViewType::value_type>::is_simd>::add_variable_data(
    ncid, name_varid_map.at(varname),
    mode_dimid, num_modes(),
    level_dimid, interface_dimid, num_levels(),
    aerosol_dimid, num_aerosols(),
    time_idx, mode_idx, spec_idx, view);
}

template <typename ViewType>
void NcWriter:: add_gas_variable_data(const std::string& varname, const size_t& time_idx,
  const int mode_idx, const int spec_idx, const ViewType& view) const {

  EKAT_REQUIRE(time_idx < num_timesteps());

  return NcWriterImpl<ViewType::Rank,
    ekat::ScalarTraits<typename ViewType::value_type>::is_simd>::add_variable_data(
      ncid, name_varid_map.at(varname),
      mode_dimid, num_modes(), level_dimid, interface_dimid, num_levels(),
      gas_dimid, num_gases(),
      time_idx, mode_idx, spec_idx, view);
}

template <> struct NcWriterImpl<1,false> {
  template <typename VT>
  static void add_variable_data(const int ncid, const int varid,
    const int mode_dimid, const int nmodes,
    const int level_dimid, const int interface_dimid, const int nlev,
    const int species_dimid, const int nspec,
    const size_t time_idx, const int mode_idx, const int spec_idx, const VT& view) {

      static_assert(VT::Rank == 1, "non-aerosol views must be rank 1.");

      int nvardims = 0;
      int retval = nc_inq_varndims(ncid, varid, &nvardims);
      EKAT_REQUIRE(nvardims == 2 && retval == NC_NOERR);
      EKAT_REQUIRE(view.extent(0) == nlev || view.extent(0) == nlev+1);

      int dimids[2] = {NC_EBADID, NC_EBADID};
      retval = nc_inq_vardimid(ncid, varid, &dimids[0]);
      EKAT_REQUIRE(retval == NC_NOERR);
      EKAT_REQUIRE(dimids[1] == level_dimid or dimids[1] == interface_dimid);

      auto hv = Kokkos::create_mirror_view(view);
      Kokkos::deep_copy(hv, view);

      const size_t array_length = (dimids[1] == level_dimid ? nlev : nlev+1);
      for (size_t i=0; i<array_length; ++i) {
        const size_t idx[2] = {time_idx, i};
#ifdef HAERO_DOUBLE_PRECISION
        retval = nc_put_var1_double(ncid, varid, &idx[0], &hv(i));
#else
        retval = nc_put_var1_float(ncid, varid, &idx[0], &hv(i));
#endif
        EKAT_ASSERT(retval==NC_NOERR);
      }
    }
};

template <> struct NcWriterImpl<1,true> {
  template <typename VT>
  static void add_variable_data(const int ncid, const int varid,
    const int mode_dimid, const int nmodes,
    const int level_dimid, const int interface_dimid, const int nlev,
    const int species_dimid, const int nspec,
    const size_t time_idx, const int mode_idx, const int spec_idx, const VT& view) {

      static_assert(VT::Rank == 1, "non-aerosol views must be rank 1.");
      int nvardims = 0;
      int retval = nc_inq_varndims(ncid, varid, &nvardims);
      EKAT_REQUIRE(nvardims == 2 && retval == NC_NOERR);

      int dimids[2] = {NC_EBADID, NC_EBADID};
      retval = nc_inq_vardimid(ncid, varid, &dimids[0]);
      EKAT_REQUIRE(retval == NC_NOERR);
      EKAT_REQUIRE(dimids[1] == level_dimid or dimids[1] == interface_dimid);

      auto hv = Kokkos::create_mirror_view(view);
      Kokkos::deep_copy(hv, view);

      const size_t array_length = (dimids[1] == level_dimid ? nlev : nlev+1);
      for (size_t i=0; i<array_length; ++i) {
        const size_t idx[2] = {time_idx, i};
        const int pack_idx = PackInfo::pack_idx(i);
        const int vec_idx = PackInfo::vec_idx(i);
#ifdef HAERO_DOUBLE_PRECISION
        nc_put_var1_double(ncid, varid, &idx[0], &hv(pack_idx)[vec_idx]);
#else
        nc_put_var1_float(ncid, varid, &idx[0], &hv(pack_idx)[vec_idx]);
#endif
        EKAT_ASSERT(retval == NC_NOERR);
      }
    }
};

template <>
struct NcWriterImpl<2, false> {
  template <typename VT>
  static void add_variable_data(const int ncid, const int varid,
    const int mode_dimid, const int nmodes,
    const int level_dimid, const int interface_dimid, const int nlev,
    const int species_dimid, const int nspec,
    const size_t time_idx, const int mode_idx, const int spec_idx, const VT& view) {

    static_assert(VT::Rank == 2, "aerosol views assumed to have rank 2");

    int nvardims = 0;
    int retval = nc_inq_varndims(ncid, varid, &nvardims);
    EKAT_REQUIRE(nvardims == 2 && retval == NC_NOERR);
    EKAT_REQUIRE((view.extent(0) == nmodes or view.extent(0) == nspec) and
      (view.extent(1) == nlev or view.extent(1) == nlev+1));

    int dimids[2] = {NC_EBADID, NC_EBADID};
    retval = nc_inq_vardimid(ncid, varid, &dimids[0]);
    EKAT_REQUIRE(retval == NC_NOERR);
    EKAT_REQUIRE(dimids[0] == mode_dimid or dimids[0] == species_dimid);
    EKAT_REQUIRE(dimids[1] == level_dimid);

    const size_t idx0 = (dimids[0] == mode_dimid ? mode_idx : spec_idx);

    const auto sv = Kokkos::subview(view, idx0, Kokkos::ALL);
    auto hsv = Kokkos::create_mirror_view(sv);
    Kokkos::deep_copy(hsv, sv);

    for (size_t i=0; i<nlev; ++i) {
      const size_t idx[3] = {time_idx, idx0, i};
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
  static void add_variable_data(const int ncid, const int varid,
    const int mode_dimid, const int nmodes,
    const int level_dimid, const int interface_dimid, const int nlev,
    const int species_dimid, const int nspec,
    const size_t time_idx, const int mode_idx, const int spec_idx, const VT& view) {

    static_assert(VT::Rank == 2, "aerosol views assumed to have rank 2");

    int nvardims = 0;
    int retval = nc_inq_varndims(ncid, varid, &nvardims);
    EKAT_REQUIRE(nvardims == 2 && retval == NC_NOERR);
    EKAT_REQUIRE((view.extent(0) == nmodes or view.extent(0) == nspec) and
      (view.extent(1) == nlev or view.extent(1) == nlev+1));

    int dimids[2] = {NC_EBADID, NC_EBADID};
    retval = nc_inq_vardimid(ncid, varid, &dimids[0]);
    EKAT_REQUIRE(retval == NC_NOERR);
    EKAT_REQUIRE(dimids[0] == mode_dimid or dimids[0] == species_dimid);
    EKAT_REQUIRE(dimids[1] == level_dimid);

    const size_t idx0 = (dimids[0] == mode_dimid ? mode_idx : spec_idx);

    const auto sv = Kokkos::subview(view, idx0, Kokkos::ALL);
    auto hsv = Kokkos::create_mirror_view(sv);
    Kokkos::deep_copy(hsv, sv);

    for (size_t i=0; i<nlev; ++i) {
      const size_t idx[3] = {time_idx, idx0, i};
      const int pack_idx = PackInfo::pack_idx(i);
      const int vec_idx = PackInfo::vec_idx(i);
#ifdef HAERO_DOUBLE_PRECISION
      retval = nc_put_var1_double(ncid, varid, &idx[0], &hsv(pack_idx)[vec_idx]);
#else
      retval = nc_put_var1_float (ncid, varid, &idx[0], &hsv(pack_idx)[vec_idx]);
#endif
      EKAT_ASSERT(retval == NC_NOERR);
    }
  }
};



} // namespace driver
} // namespace haero
#endif
