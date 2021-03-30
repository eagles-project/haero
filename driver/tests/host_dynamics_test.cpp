#include "driver/host_dynamics.hpp"
#include "driver/host_params.hpp"
#include "driver/ncwriter.hpp"
#include "driver/ncwriter_impl.hpp"
#include "haero/atmosphere.hpp"
#include "catch2/catch.hpp"
#include "ekat/util/ekat_units.hpp"
#include "haero/math_helpers.hpp"
#include <iostream>

using namespace haero;
using namespace driver;

TEST_CASE("driver dynamics", "") {

  const Real Tv0 = 300;
  const Real Gammav = 0.01;
  const Real w0 = 1;
  const Real ztop = 20E3;
  const Real tperiod = 900;
  const Real qv0 = 0.015;
  const Real qv1 = 2.5E-3;
  const AtmosphericConditions conds(Tv0, Gammav, w0, ztop, tperiod, qv0, qv1);

  const int nlev = 20;
  SECTION("height init") {
    /// In actual examples, these would come from an input yaml file
    std::vector<Real> z_vals(nlev+1);
    for (int i=0; i<nlev+1; ++i) {
      z_vals[i] = square(Real(i)/Real(nlev))*ztop;
    }

    HostDynamics zdyn(nlev);
    zdyn.init_from_interface_heights(z_vals, conds);
    std::cout << zdyn.info_string();

    /// Create a new netcdf file
    const std::string fname = "host_dynamics_test_zinit.nc";
    NcWriter writer(fname);
    writer.define_time_var();
    zdyn.nc_init_dynamics_variables(writer, conds);

    size_t time_idx = 0;
    zdyn.nc_write_data(writer, time_idx);
    Kokkos::View<PackType*> temperature("temperature", PackInfo::num_packs(nlev));
    Kokkos::View<PackType*> rel_humidity("relative_humidity", PackInfo::num_packs(nlev));
    Kokkos::View<PackType*> level_heights("level_heights", PackInfo::num_packs(nlev+1));
    auto atm = zdyn.create_atmospheric_state(temperature, rel_humidity, level_heights);
    writer.define_atm_state_vars(atm);
    writer.add_atm_state_data(atm, time_idx);

    Real t = 0.5*conds.tperiod;
    ++time_idx;
    zdyn.update(t,conds);
    zdyn.update_atmospheric_state(atm);

    writer.add_time_value(t);
    zdyn.nc_write_data(writer, time_idx);
    writer.add_atm_state_data(atm, time_idx);

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

    /// Create a new netcdf file
    const std::string fname = "host_dynamics_test_pinit.nc";
    NcWriter writer(fname);
    writer.define_time_var();
    pdyn.nc_init_dynamics_variables(writer, conds);

    size_t time_idx = 0;
    pdyn.nc_write_data(writer, time_idx);
    Kokkos::View<PackType*> temperature("temperature", PackInfo::num_packs(nlev));
    Kokkos::View<PackType*> rel_humidity("relative_humidity", PackInfo::num_packs(nlev));
    Kokkos::View<PackType*> level_heights("level_heights", PackInfo::num_packs(nlev+1));
    auto atm = pdyn.create_atmospheric_state(temperature, rel_humidity, level_heights);
    writer.define_atm_state_vars(atm);
    writer.add_atm_state_data(atm, time_idx);

    Real t = 0.5*conds.tperiod;
    ++time_idx;
    pdyn.update(t,conds);
    pdyn.update_atmospheric_state(atm);

    writer.add_time_value(t);
    pdyn.nc_write_data(writer, time_idx);
    writer.add_atm_state_data(atm, time_idx);

    writer.close();
  }
}
