#include "ncwriter.hpp"

#include <exception>
#include <iostream>
#include <sstream>

#include "haero/haero.hpp"
#include "ncwriter_impl.hpp"

namespace haero {
namespace driver {

using att_type = NcWriter::text_att_type;
using var_atts = std::vector<att_type>;

void NcWriter::open() {
  int retval = nc_create(fname.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid);
  CHECK_NCERR(retval);

  text_att_type att =
      std::make_pair("HAERO", "High-performance AEROsols standalone driver");
  add_file_attribute(att);

  att = std::make_pair("HAERO_version", std::string(version()));
  add_file_attribute(att);

  att = std::make_pair("HAERO_revision", std::string(revision()));
  add_file_attribute(att);
}

void NcWriter::add_file_attribute(const text_att_type& att_pair) const {
  const int alen = att_pair.second.size();
  int retval = nc_put_att_text(ncid, NC_GLOBAL, att_pair.first.c_str(), alen,
                               att_pair.second.c_str());
  CHECK_NCERR(retval);
}

void NcWriter::add_level_dims(const int& nlev) {
  EKAT_ASSERT(ncid != NC_EBADID);  // file is open
  EKAT_REQUIRE_MSG(level_dimid == NC_EBADID,
                   "level dimension already defined.");
  EKAT_REQUIRE_MSG(interface_dimid == NC_EBADID,
                   "interface dimension already defined.");

  int retval = nc_def_dim(ncid, "level_midpts", nlev, &level_dimid);
  CHECK_NCERR(retval);
  retval = nc_def_dim(ncid, "level_interfaces", nlev + 1, &interface_dimid);
  CHECK_NCERR(retval);
  ndims += 2;
}

void NcWriter::add_aerosol_dim(const std::vector<AerosolSpecies>& species) {
  EKAT_ASSERT(ncid != NC_EBADID);  // file is open
  EKAT_REQUIRE_MSG(aerosol_dimid == NC_EBADID,
                   "aerosol dimension already defined.");
  /// Define species dimension
  const int nspec = species.size();
  int retval = nc_def_dim(ncid, "aerosol_species", nspec, &aerosol_dimid);
  CHECK_NCERR(retval);
  ++ndims;

  /// Define coordinate variables for AerosolSpecies dimension
  int species_name_var_id = NC_EBADID;
  retval = nc_def_var(ncid, "aerosol_name", NC_STRING, 1, &aerosol_dimid,
                      &species_name_var_id);
  CHECK_NCERR(retval);
  name_varid_map.emplace("aerosol_name", species_name_var_id);

  int species_symbol_var_id = NC_EBADID;
  retval = nc_def_var(ncid, "aerosol_symbol", NC_STRING, 1, &aerosol_dimid,
                      &species_symbol_var_id);
  CHECK_NCERR(retval);
  name_varid_map.emplace("aerosol_symbol", species_symbol_var_id);

  int species_molec_weight_var_id = NC_EBADID;
  retval = nc_def_var(ncid, "aerosol_molecular_weight", NC_REAL_KIND, 1,
                      &aerosol_dimid, &species_molec_weight_var_id);
  CHECK_NCERR(retval);
  name_varid_map.emplace("aerosol_molecular_weight",
                         species_molec_weight_var_id);

  int species_dryrad_varid = NC_EBADID;
  retval = nc_def_var(ncid, "aerosol_dry_radius", NC_REAL_KIND, 1,
                      &aerosol_dimid, &species_dryrad_varid);
  CHECK_NCERR(retval);
  name_varid_map.emplace("aerosol_dry_radius", species_dryrad_varid);

  int species_dens_varid = NC_EBADID;
  retval = nc_def_var(ncid, "aerosol_density", NC_REAL_KIND, 1, &aerosol_dimid,
                      &species_dens_varid);
  CHECK_NCERR(retval);
  name_varid_map.emplace("aerosol_density", species_dens_varid);

  int species_hygro_varid = NC_EBADID;
  retval = nc_def_var(ncid, "aerosol_hygroscopicity", NC_REAL_KIND, 1,
                      &aerosol_dimid, &species_hygro_varid);
  CHECK_NCERR(retval);
  name_varid_map.emplace("aerosol_hygroscopicity", species_hygro_varid);

  const auto mwstr = ekat::units::to_string(ekat::units::kg / ekat::units::mol);
  retval = nc_put_att_text(ncid, species_molec_weight_var_id, "units",
                           mwstr.size(), mwstr.c_str());
  CHECK_NCERR(retval);
  const auto radstr = ekat::units::to_string(ekat::units::m);
  retval = nc_put_att_text(ncid, species_dryrad_varid, "units", radstr.size(),
                           radstr.c_str());
  CHECK_NCERR(retval);
  const auto densstr =
      ekat::units::to_string(ekat::units::kg / pow(ekat::units::m, 3));
  retval = nc_put_att_text(ncid, species_dens_varid, "units", densstr.size(),
                           densstr.c_str());
  CHECK_NCERR(retval);
  const auto hygstr =
      ekat::units::to_string(ekat::units::Units::nondimensional());
  retval = nc_put_att_text(ncid, species_hygro_varid, "units", hygstr.size(),
                           hygstr.c_str());

  for (int i = 0; i < nspec; ++i) {
    const size_t idx = i;
    auto name = species[i].name().c_str();
    auto symb = species[i].symbol().c_str();
    retval = nc_put_var1_string(ncid, species_name_var_id, &idx, &name);
    CHECK_NCERR(retval);
    retval = nc_put_var1_string(ncid, species_symbol_var_id, &idx, &symb);
    CHECK_NCERR(retval);
#ifdef HAERO_DOUBLE_PRECISION
    retval = nc_put_var1_double(ncid, species_molec_weight_var_id, &idx,
                                &species[i].molecular_weight);
    CHECK_NCERR(retval);
/*    retval = nc_put_var1_double(ncid, species_dryrad_varid, &idx,
                                &species[i].dry_radius);
    CHECK_NCERR(retval); */
    retval =
        nc_put_var1_double(ncid, species_dens_varid, &idx, &species[i].density);
    CHECK_NCERR(retval);
    retval = nc_put_var1_double(ncid, species_hygro_varid, &idx,
                                &species[i].hygroscopicity);
    CHECK_NCERR(retval);
#else
    retval = nc_put_var1_float(ncid, species_molec_weight_var_id, &idx,
                               &species[i].molecular_weight);
    CHECK_NCERR(retval);
    retval = nc_put_var1_float(ncid, species_dryrad_varid, &idx,
                               &species[i].dry_radius);
    CHECK_NCERR(retval);
    retval =
        nc_put_var1_float(ncid, species_dens_varid, &idx, &species[i].density);
    CHECK_NCERR(retval);
    retval = nc_put_var1_float(ncid, species_hygro_varid, &idx,
                               &species[i].hygroscopicity);
    CHECK_NCERR(retval);
#endif
  }
}

void NcWriter::add_gas_dim(const std::vector<GasSpecies>& gases) {
  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_REQUIRE_MSG(gas_dimid == NC_EBADID, "gas dimension already defined.");
  const int nspec = gases.size();
  int retval = nc_def_dim(ncid, "gas_species", nspec, &gas_dimid);
  CHECK_NCERR(retval);
  ++ndims;

  int gas_name_var_id = NC_EBADID;
  retval =
      nc_def_var(ncid, "gas_name", NC_STRING, 1, &gas_dimid, &gas_name_var_id);
  CHECK_NCERR(retval);
  name_varid_map.emplace("gas_name", gas_name_var_id);

  int gas_symb_var_id = NC_EBADID;
  retval = nc_def_var(ncid, "gas_symbol", NC_STRING, 1, &gas_dimid,
                      &gas_symb_var_id);
  CHECK_NCERR(retval);
  name_varid_map.emplace("gas_symbol", gas_symb_var_id);

  int gas_mw_var_id = NC_EBADID;
  retval = nc_def_var(ncid, "gas_molecular_weight", NC_REAL_KIND, 1, &gas_dimid,
                      &gas_mw_var_id);
  CHECK_NCERR(retval);
  name_varid_map.emplace("gas_molecular_weight", gas_mw_var_id);

  const auto mwstr = ekat::units::to_string(ekat::units::kg / ekat::units::mol);
  retval = nc_put_att_text(ncid, gas_mw_var_id, "units", mwstr.size(),
                           mwstr.c_str());
  CHECK_NCERR(retval);

  for (int i = 0; i < nspec; ++i) {
    const size_t idx = i;
    auto name = gases[i].name().c_str();
    auto symb = gases[i].symbol().c_str();
    retval = nc_put_var1_string(ncid, gas_name_var_id, &idx, &name);
    CHECK_NCERR(retval);
    retval = nc_put_var1_string(ncid, gas_symb_var_id, &idx, &symb);
    CHECK_NCERR(retval);
#ifdef HAERO_DOUBLE_PRECISION
    retval = nc_put_var1_double(ncid, gas_mw_var_id, &idx,
                                &gases[i].molecular_weight);
    CHECK_NCERR(retval);
#else
    retval = nc_put_var1_float(ncid, gas_mw_var_id, &idx,
                               &gases[i].molecular_weight);
    CHECK_NCERR(retval);
#endif
  }
}

void NcWriter::add_mode_dim(const std::vector<Mode>& modes) {
  EKAT_ASSERT(ncid != NC_EBADID);  // file is open
  EKAT_REQUIRE_MSG(mode_dimid == NC_EBADID, "mode dimension already defined.");
  /// Define mode dimension
  const int nmodes = modes.size();
  int retval = nc_def_dim(ncid, "modes", nmodes, &mode_dimid);
  CHECK_NCERR(retval);
  ++ndims;

  /// define coordinate variables for Mode dimension
  int mode_name_var_id = NC_EBADID;
  retval = nc_def_var(ncid, "mode_names", NC_STRING, 1, &mode_dimid,
                      &mode_name_var_id);
  CHECK_NCERR(retval);
  name_varid_map.emplace("mode_names", mode_name_var_id);

  int mode_min_varid = NC_EBADID;
  int mode_max_varid = NC_EBADID;
  int mode_std_dev_varid = NC_EBADID;
  retval = nc_def_var(ncid, "mode_minimum_diameter", NC_REAL_KIND, 1,
                      &mode_dimid, &mode_min_varid);
  CHECK_NCERR(retval);
  name_varid_map.emplace("mode_minimum_diameter", mode_min_varid);

  retval = nc_def_var(ncid, "mode_maximum_diameter", NC_REAL_KIND, 1,
                      &mode_dimid, &mode_max_varid);
  CHECK_NCERR(retval);
  name_varid_map.emplace("mode_maximum_diameter", mode_max_varid);

  retval = nc_def_var(ncid, "mode_mean_std_deviation", NC_REAL_KIND, 1,
                      &mode_dimid, &mode_std_dev_varid);
  CHECK_NCERR(retval);
  name_varid_map.emplace("mode_mean_std_deviation", mode_std_dev_varid);

  /// define units for each variable
  const std::string meter_str = ekat::units::to_string(ekat::units::m);
  retval = nc_put_att_text(ncid, mode_min_varid, "units", meter_str.size(),
                           meter_str.c_str());
  CHECK_NCERR(retval);
  retval = nc_put_att_text(ncid, mode_max_varid, "units", meter_str.size(),
                           meter_str.c_str());
  CHECK_NCERR(retval);
  const std::string nondim_str =
      ekat::units::to_string(ekat::units::Units::nondimensional());
  retval = nc_put_att_text(ncid, mode_std_dev_varid, "units", nondim_str.size(),
                           nondim_str.c_str());
  CHECK_NCERR(retval);

  /// add variable data
  for (int i = 0; i < nmodes; ++i) {
    const size_t idx = i;
    auto name = modes[i].name().c_str();
    retval = nc_put_var1_string(ncid, mode_name_var_id, &idx, &name);
    CHECK_NCERR(retval);
#ifdef HAERO_DOUBLE_PRECISION
    retval =
        nc_put_var1_double(ncid, mode_min_varid, &idx, &modes[i].min_diameter);
    CHECK_NCERR(retval);
    retval =
        nc_put_var1_double(ncid, mode_max_varid, &idx, &modes[i].max_diameter);
    CHECK_NCERR(retval);
    retval = nc_put_var1_double(ncid, mode_std_dev_varid, &idx,
                                &modes[i].mean_std_dev);
#else
    retval =
        nc_put_var1_float(ncid, mode_min_varid, &idx, &modes[i].min_diameter);
    CHECK_NCERR(retval);
    retval =
        nc_put_var1_float(ncid, mode_max_varid, &idx, &modes[i].max_diameter);
    CHECK_NCERR(retval);
    retval = nc_put_var1_float(ncid, mode_std_dev_varid, &idx,
                               &modes[i].mean_std_dev);
#endif
  }
}

void NcWriter::define_atm_state_vars(const Atmosphere& atm) {
  const var_atts Tatts = {std::make_pair("cf_long_name", "air_temperature"),
                          std::make_pair("amip_short_name", "ta"),
                          std::make_pair("short_name", "T")};
  define_level_var("temperature", ekat::units::K, atm.temperature, Tatts);

  const var_atts rhatts = {std::make_pair("cf_long_name", "relative_humidity"),
                           std::make_pair("amip_short_name", "hur"),
                           std::make_pair("haero_short_name", "rel_humidity")};
  define_level_var("relative_humidity", ekat::units::Units::nondimensional(),
                   atm.relative_humidity, rhatts);

  const var_atts hatts = {std::make_pair("cf_long_name", "geopotential_height"),
                          std::make_pair("haero_short_name", "z")};
  define_interface_var("geopotential_height", ekat::units::m, atm.height,
                       hatts);
}

void NcWriter::add_atm_state_data(const Atmosphere& atm,
                                  const size_t time_idx) {
  const int null_idx = -1;
  add_variable_data("temperature", time_idx, null_idx, null_idx,
                    atm.temperature);
  add_variable_data("relative_humidity", time_idx, null_idx, null_idx,
                    atm.relative_humidity);
  add_variable_data("geopotential_height", time_idx, null_idx, null_idx,
                    atm.height);
}

void NcWriter::add_time_value(const Real& t) const {
  EKAT_ASSERT(time_dimid != NC_EBADID);
  const int timeid = name_varid_map.at("time");
  size_t next_time_index = 0;
  int retval = nc_inq_dimlen(ncid, time_dimid, &next_time_index);
  CHECK_NCERR(retval);
#ifdef HAERO_DOUBLE_PRECISION
  nc_put_var1_double(ncid, timeid, &next_time_index, &t);
#else
  nc_put_var1_float(ncid, timeid, &next_time_index, &t);
#endif
  CHECK_NCERR(retval);
}

void NcWriter::add_level_variable_data(const std::string& varname,
                                       const size_t& time_index,
                                       const std::vector<Real>& data) const {
  const int varid = name_varid_map.at(varname);
  size_t nlev = 0;
  int retval = nc_inq_dimlen(ncid, level_dimid, &nlev);
  CHECK_NCERR(retval);
  EKAT_REQUIRE_MSG(data.size() == nlev,
                   "add_level_variable_data called with data.size() != nlev");
  size_t nsteps = 0;
  retval = nc_inq_dimlen(ncid, time_dimid, &nsteps);
  CHECK_NCERR(retval);
  EKAT_REQUIRE_MSG(
      time_index < nsteps,
      "add_level_variable_data called for out-of-bounds time index");
  size_t start[2] = {time_index, 0};
  size_t count[2] = {1, nlev};
#ifdef HAERO_DOUBLE_PRECISION
  retval = nc_put_vara_double(ncid, varid, start, count, &data[0]);
#else
  retval = nc_put_vara_float(ncid, varid, start, count, &data[0]);
#endif
  CHECK_NCERR(retval);
}

void NcWriter::add_interface_variable_data(
    const std::string& varname, const size_t& time_index,
    const std::vector<Real>& data) const {
  const int varid = name_varid_map.at(varname);
  size_t ninterfaces = 0;
  int retval = nc_inq_dimlen(ncid, interface_dimid, &ninterfaces);
  CHECK_NCERR(retval);
  EKAT_REQUIRE_MSG(
      data.size() == ninterfaces,
      "add_interface_variable_data called with data.size() != ninterfaces");
  size_t nsteps = 0;
  retval = nc_inq_dimlen(ncid, time_dimid, &nsteps);
  EKAT_REQUIRE_MSG(
      time_index < nsteps,
      "add_interface_variable_data called for out-of-bounds time index");
  size_t start[2] = {time_index, 0};
  size_t count[2] = {1, ninterfaces};
#ifdef HAERO_DOUBLE_PRECISION
  retval = nc_put_vara_double(ncid, varid, start, count, &data[0]);
#else
  retval = nc_put_vara_float(ncid, varid, start, count, &data[0]);
#endif
  CHECK_NCERR(retval);
}

void NcWriter::add_time_dependent_scalar_value(const std::string& name,
                                               const size_t time_idx,
                                               const Real val) const {
  EKAT_ASSERT(time_idx < num_timesteps());
  const int varid = name_varid_map.at(name);

#ifdef HAERO_DOUBLE_PRECISION
  int retval = nc_put_var1_double(ncid, varid, &time_idx, &val);
#else
  int retval = nc_put_var1_float(ncid, varid, &time_idx, &val);
#endif
  CHECK_NCERR(retval);
}

void NcWriter::add_time_dim() {
  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_REQUIRE_MSG(time_dimid == NC_EBADID, "time dimension already defined.");
  int retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dimid);
  CHECK_NCERR(retval);
  ndims++;
}

void NcWriter::close() {
  int retval = nc_close(ncid);
  CHECK_NCERR(retval);
}

void NcWriter::define_scalar_var(const std::string& name,
                                 const ekat::units::Units& units,
                                 const std::vector<text_att_type>& var_atts,
                                 const Real val) {
  EKAT_ASSERT(ncid != NC_EBADID);
  int* not_used;
  int varid = NC_EBADID;
  int retval =
      nc_def_var(ncid, name.c_str(), NC_REAL_KIND, 0, not_used, &varid);
  CHECK_NCERR(retval);
  name_varid_map.emplace(name, varid);

  for (int i = 0; i < var_atts.size(); ++i) {
    retval =
        nc_put_att_text(ncid, varid, var_atts[i].first.c_str(),
                        var_atts[i].second.size(), var_atts[i].second.c_str());
    CHECK_NCERR(retval);
  }
  retval = nc_put_var(ncid, varid, &val);
  CHECK_NCERR(retval);
  const auto unitstr = ekat::units::to_string(units);
  retval =
      nc_put_att_text(ncid, varid, "units", unitstr.size(), unitstr.c_str());
  CHECK_NCERR(retval);
}

void NcWriter::define_time_dependent_scalar_var(
    const std::string& name, const ekat::units::Units& units,
    const std::vector<text_att_type>& var_atts) {
  EKAT_ASSERT(ncid != NC_EBADID);
  EKAT_REQUIRE(time_dimid != NC_EBADID);
  int varid = NC_EBADID;
  int retval =
      nc_def_var(ncid, name.c_str(), NC_REAL_KIND, 1, &time_dimid, &varid);
  CHECK_NCERR(retval);
  name_varid_map.emplace(name, varid);

  for (int i = 0; i < var_atts.size(); ++i) {
    retval =
        nc_put_att_text(ncid, varid, var_atts[i].first.c_str(),
                        var_atts[i].second.size(), var_atts[i].second.c_str());
    CHECK_NCERR(retval);
  }
  const auto unitstr = ekat::units::to_string(units);
  retval =
      nc_put_att_text(ncid, varid, "units", unitstr.size(), unitstr.c_str());
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

void NcWriter::define_level_var(const std::string& name,
                                const ekat::units::Units& units,
                                const std::vector<text_att_type>& atts) {
  EKAT_ASSERT(time_dimid != NC_EBADID && level_dimid != NC_EBADID);

  int varid = NC_EBADID;
  const int m_ndims = 2;
  const int dimids[2] = {time_dimid, level_dimid};
  int retval =
      nc_def_var(ncid, name.c_str(), NC_REAL_KIND, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  const std::string unit_str = ekat::units::to_string(units);
  retval =
      nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
  if (!atts.empty()) {
    for (int i = 0; i < atts.size(); ++i) {
      retval = nc_put_att_text(ncid, varid, atts[i].first.c_str(),
                               atts[i].second.size(), atts[i].second.c_str());
      CHECK_NCERR(retval);
    }
  }
  name_varid_map.emplace(name, varid);
}

void NcWriter::define_time_var(const ekat::units::Units& units, const Real t0) {
  EKAT_ASSERT(time_dimid != NC_EBADID);

  int varid = NC_EBADID;
  int retval = nc_def_var(ncid, "time", NC_REAL_KIND, 1, &time_dimid, &varid);
  CHECK_NCERR(retval);
  const std::string unit_str = ekat::units::to_string(units);
  retval =
      nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
  name_varid_map.emplace("time", varid);
  add_time_value(t0);
}

void NcWriter::define_interface_var(const std::string& name,
                                    const ekat::units::Units& units,
                                    const std::vector<text_att_type>& atts) {
  EKAT_ASSERT(time_dimid != NC_EBADID and interface_dimid != NC_EBADID);

  int varid = NC_EBADID;
  const int m_ndims = 2;
  const int dimids[2] = {time_dimid, interface_dimid};
  int retval =
      nc_def_var(ncid, name.c_str(), NC_REAL_KIND, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  const std::string unit_str = ekat::units::to_string(units);
  retval =
      nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
  if (!atts.empty()) {
    for (int i = 0; i < atts.size(); ++i) {
      retval = nc_put_att_text(ncid, varid, atts[i].first.c_str(),
                               atts[i].second.size(), atts[i].second.c_str());
      CHECK_NCERR(retval);
    }
  }
  name_varid_map.emplace(name, varid);
}

void NcWriter::define_const_1dvar(const std::string& name,
                                  const ekat::units::Units& units,
                                  const std::vector<Real>& vals,
                                  const std::vector<text_att_type>& atts) {
  EKAT_REQUIRE_MSG(
      vals.size() == num_levels() or vals.size() == num_levels() + 1,
      "coordinate variables must be either levels or interfaces.");

  int varid = NC_EBADID;
  int retval = nc_def_var(
      ncid, name.c_str(), NC_REAL_KIND, 1,
      (vals.size() == num_levels() ? &level_dimid : &interface_dimid), &varid);
  CHECK_NCERR(retval);
  const std::string unit_str = ekat::units::to_string(units);
  retval =
      nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
  name_varid_map.emplace(name, varid);
  for (int i = 0; i < atts.size(); ++i) {
    retval = nc_put_att_text(ncid, varid, atts[i].first.c_str(),
                             atts[i].second.size(), atts[i].second.c_str());
    CHECK_NCERR(retval);
  }
  const size_t start = 0;
  const size_t count = vals.size();
#ifdef HAERO_DOUBLE_PRECISION
  retval = nc_put_vara_double(ncid, varid, &start, &count, &vals[0]);
#else
  retval = nc_put_vara_float(ncid, varid, &start, &count, &vals[0]);
#endif
  CHECK_NCERR(retval);
}

int NcWriter::num_levels() const {
  size_t nl;
  int retval = nc_inq_dimlen(ncid, level_dimid, &nl);
  CHECK_NCERR(retval);
  return int(nl);
}

int NcWriter::num_modes() const {
  size_t nm;
  int retval = nc_inq_dimlen(ncid, mode_dimid, &nm);
  if (retval == NC_EBADDIM) {
    nm = 0;
  } else {
    CHECK_NCERR(retval);
  }
  return int(nm);
}

int NcWriter::num_aerosols() const {
  size_t ns;
  int retval = nc_inq_dimlen(ncid, aerosol_dimid, &ns);
  if (retval == NC_EBADDIM) {
    ns = 0;
  } else {
    CHECK_NCERR(retval);
  }
  return int(ns);
}

int NcWriter::num_gases() const {
  size_t ns;
  int retval = nc_inq_dimlen(ncid, gas_dimid, &ns);
  if (retval == NC_EBADDIM) {
    ns = 0;
  } else {
    CHECK_NCERR(retval);
  }
  return int(ns);
}

int NcWriter::num_timesteps() const {
  size_t ns;
  int retval = nc_inq_dimlen(ncid, time_dimid, &ns);
  CHECK_NCERR(retval);
  return int(ns);
}

void NcWriter::define_modal_var(const std::string& name,
                                const ekat::units::Units& units,
                                const std::vector<text_att_type>& atts) {
  EKAT_ASSERT(time_dimid != NC_EBADID && mode_dimid != NC_EBADID &&
              level_dimid != NC_EBADID);

  const int m_ndims = 3;
  const int dimids[3] = {time_dimid, mode_dimid, level_dimid};
  int varid = NC_EBADID;
  int retval =
      nc_def_var(ncid, name.c_str(), NC_REAL_KIND, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  const std::string unit_str = ekat::units::to_string(units);
  retval =
      nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
  if (!atts.empty()) {
    for (int i = 0; i < atts.size(); ++i) {
      retval = nc_put_att_text(ncid, varid, atts[i].first.c_str(),
                               atts[i].second.size(), atts[i].second.c_str());
      CHECK_NCERR(retval);
    }
  }
  name_varid_map.emplace(name, varid);
}

void NcWriter::define_aerosol_var(const std::string& name,
                                  const ekat::units::Units& units,
                                  const std::vector<text_att_type>& atts) {
  EKAT_ASSERT(time_dimid != NC_EBADID && aerosol_dimid != NC_EBADID &&
              level_dimid != NC_EBADID);

  const int m_ndims = 3;
  const int dimids[3] = {time_dimid, aerosol_dimid, level_dimid};
  int varid = NC_EBADID;
  int retval =
      nc_def_var(ncid, name.c_str(), NC_REAL_KIND, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  const auto unit_str = ekat::units::to_string(units);
  retval =
      nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
  if (!atts.empty()) {
    for (int i = 0; i < atts.size(); ++i) {
      retval = nc_put_att_text(ncid, varid, atts[i].first.c_str(),
                               atts[i].second.size(), atts[i].second.c_str());
      CHECK_NCERR(retval);
    }
  }
  name_varid_map.emplace(name, varid);
}

void NcWriter::define_gas_var(const std::string& name,
                              const ekat::units::Units& units,
                              const std::vector<text_att_type>& atts) {
  EKAT_ASSERT(time_dimid != NC_EBADID && gas_dimid != NC_EBADID &&
              level_dimid != NC_EBADID);

  const int m_ndims = 3;
  const int dimids[3] = {time_dimid, gas_dimid, level_dimid};
  int varid = NC_EBADID;
  int retval =
      nc_def_var(ncid, name.c_str(), NC_REAL_KIND, m_ndims, dimids, &varid);
  CHECK_NCERR(retval);
  const auto unit_str = ekat::units::to_string(units);
  retval =
      nc_put_att_text(ncid, varid, "units", unit_str.size(), unit_str.c_str());
  CHECK_NCERR(retval);
  if (!atts.empty()) {
    for (int i = 0; i < atts.size(); ++i) {
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

void NcWriter::handle_errcode(const int& ec, const std::string& file,
                              const std::string& fn, const int& line) const {
  std::ostringstream ss;
  ss << "NcWriter error in file: " << file << ", function: " << fn
     << ", at line: " << line << '\n';
  ss << "\terror " << ec << " decodes to: ";
  switch (ec) {
    case (NC_NOERR): {
      // no error: should not have called this routine
      EKAT_ASSERT(ec != NC_NOERR);
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
    case (NC_EBADDIM): {
      ss << "invalid dimension id or name";
      break;
    }
    default: {
      ss << "unknown netcdf error";
    }
  }
  throw std::runtime_error(ss.str());
}

}  // namespace driver
}  // namespace haero
