#include "driver/host_dynamics.hpp"
#include "driver/host_params.hpp"
#include "driver/ncwriter.hpp"
#include "driver/ncwriter_impl.hpp"
#include "haero/atmosphere.hpp"
#include "catch2/catch.hpp"
#include "ekat/util/ekat_units.hpp"
#include "haero/math_helpers.hpp"
#include <iostream>
#include <iomanip>

using namespace haero;
using namespace driver;

struct HydrostaticBalance {

  int nerr;

  HydrostaticBalance() : nerr(0) {}

  void run_test(const HostDynamics& dyn, const AtmosphericConditions& ac);
};

struct UniformThicknessHeight {
  int nerr;
  Real dzval;
  UniformThicknessHeight(const Real layer_depth) : nerr(0), dzval(layer_depth) {}
  void run_test(const HostDynamics& dyn);
};

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

  HydrostaticBalance hbtest;
  UniformThicknessHeight onedz(ztop/nlev);

  SECTION("height init -- uniform heights") {
    HostDynamics zdyn(nlev);
    zdyn.init_from_uniform_heights(nlev, conds);
    std::cout << zdyn.info_string();
    hbtest.run_test(zdyn, conds);
    onedz.run_test(zdyn);

    REQUIRE(hbtest.nerr == 0);
    REQUIRE(onedz.nerr == 0);

    /// Create a new netcdf file
    const std::string fname = "host_dynamics_test_zinit_unif.nc";
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

  SECTION("height init -- specified heights") {
    /// In actual examples, these would come from an input yaml file
    std::vector<Real> z_vals(nlev+1);
    for (int i=0; i<nlev+1; ++i) {
      z_vals[i] = square(Real(i)/Real(nlev))*ztop;
    }

    HostDynamics zdyn(nlev);
    zdyn.init_from_interface_heights(z_vals, conds);
    std::cout << zdyn.info_string();
    hbtest.run_test(zdyn, conds);

    REQUIRE(hbtest.nerr == 0);

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

  SECTION("pressure_init -- uniform pressures") {
    HostDynamics pdyn(nlev);
    pdyn.init_from_uniform_pressures(nlev, conds);
    std::cout << pdyn.info_string();
    hbtest.run_test(pdyn, conds);

    REQUIRE(hbtest.nerr == 0);

    /// Create a new netcdf file
    const std::string fname = "host_dynamics_test_pinit_unif.nc";

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

  SECTION("pressure_init -- specified pressures") {
    /// In actual examples, these would come from an input yaml file
    const std::vector<Real> p_vals = {20000, 30000, 40000, 50000, 60000, 70000, 85000, 92500, 100000};
    /// Assume input is interfaces
    const int nlev = p_vals.size()-1;

    HostDynamics pdyn(nlev);
    pdyn.init_from_interface_pressures(p_vals, conds);
    std::cout << pdyn.info_string();
    hbtest.run_test(pdyn, conds);

    REQUIRE(hbtest.nerr == 0);

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

void HydrostaticBalance::run_test(const HostDynamics& dyn, const AtmosphericConditions& ac) {
  using namespace constants;
  nerr = 0;

  const auto phi = dyn.phi;
  const auto p = dyn.p;

  Kokkos::parallel_reduce("HydrostaticBalance::run_test", dyn.nlev(),
    KOKKOS_LAMBDA (const int k, int& errct) {
      // Taylor et al. 2020 fig. 1 level idx = k+1
      const auto pack_idx = PackInfo::pack_idx(k);
      const auto vec_idx = PackInfo::vec_idx(k);

      const int kphalf_idx = k+1; // array idx of interface k + 1/2
      const int kmhalf_idx = k; // array idx of interface k - 1/2
      const auto kphalf_pack_idx = PackInfo::pack_idx(kphalf_idx);
      const auto kphalf_vec_idx = PackInfo::vec_idx(kphalf_idx);
      const auto kmhalf_pack_idx = PackInfo::pack_idx(kmhalf_idx);
      const auto kmhalf_vec_idx = PackInfo::vec_idx(kmhalf_idx);

      const Real phimid = 0.5*(phi(kmhalf_pack_idx)[kmhalf_vec_idx] +
                                        phi(kphalf_pack_idx)[kphalf_vec_idx]);
      const Real zmid = phimid/gravity;
      const Real pres = hydrostatic_pressure_at_height(zmid, ac);

      if (!FloatingPoint<Real>::zero(pres - p(pack_idx)[vec_idx])) {
        ++errct;
        printf("at level %d: pres = %f, p = %f\n", k, pres, p(pack_idx)[vec_idx]);
      }
    }, nerr);
}

void UniformThicknessHeight::run_test(const HostDynamics& dyn) {
  const auto dz = dyn.dz;
  nerr = 0;
  Kokkos::parallel_reduce("UniformThicknessHeight::run_test", dyn.nlev(),
    KOKKOS_LAMBDA (const int k, int& errct) {
      const int pack_idx_k = PackInfo::pack_idx(k);
      const int vec_idx_k = PackInfo::vec_idx(k);
      if (!FloatingPoint<Real>::equiv(dz(pack_idx_k)[vec_idx_k], dzval, 2E-12)) {
//         std::cout << "at level " << k << ": abs(dz-dzval) = " << std::setprecision(16)
//                   << std::abs(dz(pack_idx_k)[vec_idx_k]-dzval) << "\n";
        ++errct;
      }
    }, nerr);
}
