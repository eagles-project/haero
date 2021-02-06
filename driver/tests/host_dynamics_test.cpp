#include "driver/host_state.hpp"
#include "driver/host_params.hpp"
#include "driver/ncwriter.hpp"
#include "driver/ncwriter_impl.hpp"
#include "catch2/catch.hpp"
#include "ekat/util/ekat_units.hpp"
#include <iostream>

using namespace haero;
using namespace driver;

using att_type = NcWriter::text_att_type;
using var_atts = std::vector<att_type>;

TEST_CASE("driver dynamics", "") {
  const AtmosphericConditions conds;
  SECTION("height init") {
    /// In actual examples, these would come from an input yaml file
    const std::vector<Real> z_vals = {10000, 9000, 8000, 7000,
                                        6000, 5000, 4000, 3000,
                                        2000, 1500, 1000,  750,
                                         500,  400,  300,  200,
                                         100,   80,   60,   40,
                                          20,   10, 0};
    /// Assume input is interfaces
    const int nlev = z_vals.size()-1;
    
    HostDynamics zdyn(nlev);
    zdyn.init_from_interface_heights(z_vals, conds);
    std::cout << zdyn.info_string();
    
    /// Create a new netcdf file
    const std::string fname = "host_dynamics_test_zinit.nc";
    NcWriter writer(fname);
    writer.add_level_dims(nlev);
    
    // Interface variables
    const var_atts w_atts = {std::make_pair("cf_long_name", "upward_air_velocity"),
      std::make_pair("short_name", "w")};
    const auto w_units = ekat::units::m/ekat::units::s;

    const var_atts phi_atts = {std::make_pair("cf_long_name", "geopotential"),
      std::make_pair("short_name","phi")};
    const auto phi_units = ekat::units::pow(ekat::units::m,2)*ekat::units::pow(ekat::units::s,-2);
    
    writer.define_interface_var("vertical_velocity", w_units, zdyn.w, w_atts);
    writer.define_interface_var("geopotential", phi_units, zdyn.phi, phi_atts);
    
    // level variables
    const var_atts thetav_atts = {std::make_pair("cf_long_name", "null"),
      std::make_pair("haero_long_name", "virtual_potential_temperature"),
      std::make_pair("short_name", "theta_v")};
    const auto thetav_units = ekat::units::K;
    
    const var_atts p_atts = {std::make_pair("cf_long_name", "air_pressure"),
      std::make_pair("amip_short_name", "plev"),
      std::make_pair("short_name", "p")};
    const auto p_units = ekat::units::Pa;
    
    const var_atts qv_atts = {std::make_pair("cf_long_name", "humidity_mixing_ratio"),
      std::make_pair("haero_long_name", "water_vapor_mass_mixing_ratio"),
      std::make_pair("short_name", "q_v")};
    const auto qv_units = ekat::units::kg / ekat::units::kg;
    
    const var_atts rho_atts = {std::make_pair("cf_long_name", "air_density"), 
      std::make_pair("short_name", "rho")};
    const auto rho_units = ekat::units::kg * ekat::units::pow(ekat::units::m,-3);
    
    writer.define_level_var("density", rho_units, zdyn.rho, rho_atts);
    writer.define_level_var("thetav", thetav_units, zdyn.thetav, thetav_atts);
    writer.define_level_var("qv", qv_units, zdyn.qv, qv_atts);
    writer.define_level_var("p", p_units, zdyn.p, p_atts);
    
    // surface variables
    const var_atts ps_atts = {std::make_pair("cf_long_name", "surface_air_pressure"),
      std::make_pair("short_name", "psurf"), std::make_pair("amip_short_name", "ps")};
    const auto ps_units = ekat::units::Pa;
    writer.define_time_dependent_scalar_var("surface_pressure", ps_units, ps_atts);
    
    const var_atts p0_atts = {
    std::make_pair("cf_long_name", "reference_air_pressure_for_atmospheric_vertical_coordinate"),
    std::make_pair("short_name", "p0")};
    writer.define_scalar_var("p0", ekat::units::Pa, p0_atts, AtmosphericConditions::pref);

    const var_atts Tv0_atts = {std::make_pair("cf_long_name", "null"),
      std::make_pair("haero_long_name", "reference_virtual_temperature")};
    const auto Tv0_units = ekat::units::K;
    writer.define_scalar_var("Tv0", Tv0_units, Tv0_atts, conds.Tv0);

  const var_atts Gammav_atts = {std::make_pair("cf_long_name", "null"),
    std::make_pair("haero_long_name", "virtual_temperature_lapse_rate"),
    std::make_pair("short_name", "Gamma_v")};
  const auto Gammav_units = ekat::units::K / ekat::units::m;
  writer.define_scalar_var("Tv_lapse_rate", Gammav_units, Gammav_atts,
    conds.Gammav);

  const var_atts qv0_atts = {std::make_pair("cf_long_name", "null"),
    std::make_pair("haero_long_name", "surface_water_vapor_mixing_ratio"),
    std::make_pair("short_name", "qv0")};
  writer.define_scalar_var("qv0", ekat::units::Units::nondimensional(), qv0_atts,
    conds.qv0);

  const var_atts qv1_atts = {std::make_pair("cf_long_name", "null"),
    std::make_pair("haero_long_name", "water_vapor_mixing_ratio_decay_rate"),
    std::make_pair("short_name", "qv1")};
  writer.define_scalar_var("qv1", ekat::units::pow(ekat::units::m, -1), qv1_atts,
    conds.qv1);
    
    const var_atts ztop_atts = {std::make_pair("cf_long_name", "altitude_at_top_of_atmosphere_model"),
      std::make_pair("short_name", "z_top")};
    writer.define_scalar_var("ztop", ekat::units::m, ztop_atts, conds.ztop);
    
    const var_atts w0_atts = {std::make_pair("cf_long_name", "null"), 
      std::make_pair("haero_long_name", "host_dynamics_maximum_velocity"),
      std::make_pair("short_name", "w0")};
    writer.define_scalar_var("w0", ekat::units::m/ekat::units::s, w0_atts, conds.w0);
    
    const var_atts tp_atts = {std::make_pair("cf_long_name", "null"),
      std::make_pair("haero_long_name", "host_dynamics_oscillation_period"),
      std::make_pair("short_name", "t_period")};
    writer.define_scalar_var("tperiod", ekat::units::s, tp_atts, conds.tperiod);
    
    writer.close();
  }

  SECTION("pressure_init") {
    /// In actual examples, these would come from an input yaml file
    const std::vector<Real> p_vals = {20000, 30000, 40000, 50000, 60000, 70000, 85000, 92500, 100000};
    /// Assume input is interfaces
    const int nlev = p_vals.size()-1;
    
    HostDynamics pdyn(nlev);
    pdyn.init_from_interface_pressures(p_vals, conds);
    std::cout << pdyn.info_string();
  }
}