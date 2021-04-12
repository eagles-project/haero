#include "driver/host_dynamics.hpp"
#include "driver/host_params.hpp"
#include "driver/ncwriter.hpp"
#include "driver/ncwriter_impl.hpp"
#include "haero/atmosphere.hpp"
#include "catch2/catch.hpp"
#include "haero/math_helpers.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <sstream>

using namespace haero;
using namespace driver;

struct HydrostaticBalanceTest {
  int nerr;
  HydrostaticBalanceTest() : nerr(0) {}
  void run_test(const HostDynamics& dyn, const AtmosphericConditions& ac,
    const Real tol = FloatingPoint<Real>::zero_tol);
};

struct UniformThicknessHeightTest {
  int nerr;
  Real dzval;
  UniformThicknessHeightTest(const Real layer_depth) : nerr(0), dzval(layer_depth) {}
  void run_test(const HostDynamics& dyn, const Real tol= FloatingPoint<Real>::zero_tol);
};

struct UniformThicknessPressureTest {
  int nerr;
  Real dpval;
  UniformThicknessPressureTest(const Real layer_depth) : nerr(0), dpval(layer_depth) {}
  void run_test(const HostDynamics& dyn, const Real tol= FloatingPoint<Real>::zero_tol);
};

struct HypsometricLevelsTest {
  int nerr;
  HypsometricLevelsTest() : nerr(0) {}
  void run_test(const HostDynamics& dyn, const AtmosphericConditions& ac,
  const Real tol= FloatingPoint<Real>::zero_tol);
};

struct VerticalConvergenceTests {
  int nlevstart;
  int ntests;
  std::vector<Real> hypso_maxres;
  std::vector<Real> hypso_rate;
  std::vector<Real> hydro_maxres;
  std::vector<Real> hydro_rate;
  std::vector<Real> ps_res;
  std::vector<Real> ps_rate;
  std::vector<Real> ztop_res;
  std::vector<Real> ztop_rate;
  std::vector<int>  nlevs;

  VerticalConvergenceTests(const int nlev0 = 10, const int nt = 5) :
    nlevstart(nlev0),
    ntests(nt),
    hypso_maxres(nt,0),
    hypso_rate(nt,0),
    hydro_maxres(nt,0),
    hydro_rate(nt,0),
    ps_res(nt,0),
    ps_rate(nt,0),
    ztop_res(nt,0),
    ztop_rate(nt,0),
    nlevs(nt) {
      EKAT_ASSERT(nlev0 > 0);
      EKAT_ASSERT(nt > 0);
      for (int i=0; i<nt; ++i) {
        nlevs[i] = std::pow(2,i)*nlevstart;
      }
    }

  Real run_hypsometric_test(const HostDynamics& dyn, const AtmosphericConditions& ac);
  Real run_hydrostatic_test(const HostDynamics& dyn, const AtmosphericConditions& ac);
  Real run_ps_test(const HostDynamics& dyn, const AtmosphericConditions& ac);
  Real run_ztop_test(const HostDynamics& dyn, const AtmosphericConditions& ac);
  void compute_appx_conv_rates();
  std::string info_string() const;
};

TEST_CASE("driver dynamics", "") {

  const Real Tv0 = 300;
  const Real Gammav = 0.01;
  const Real w0 = 1;
  const Real tperiod = 900;
  const Real qv0 = 0.015;
  const Real qv1 = 2.5E-3;

  HydrostaticBalanceTest hbtest;
  HypsometricLevelsTest hypsotest;

  SECTION("height init -- uniform heights") {
    const int nlev = 20;
    const Real ztop = 20E3;
    const AtmosphericConditions conds(Tv0, Gammav, w0, ztop, tperiod, qv0, qv1);
    std::cout << conds.info_string();

    UniformThicknessHeightTest onedz(ztop/nlev);

    HostDynamics zdyn(nlev);
    zdyn.init_from_uniform_heights(conds);
    std::cout << zdyn.info_string();
    hbtest.run_test(zdyn, conds);
    onedz.run_test(zdyn, 12*FloatingPoint<Real>::zero_tol);
    hypsotest.run_test(zdyn, conds, 1.5e-2);

    REQUIRE(hbtest.nerr == 0);
    REQUIRE(onedz.nerr == 0);
    REQUIRE(hypsotest.nerr == 0);

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
    const int nlev = 20;
    /// In actual examples, these would come from an input yaml file
    std::vector<Real> z_vals(nlev+1);
    const Real ztop = 20E3;
    for (int i=0; i<nlev+1; ++i) {
      z_vals[i] = square(Real(i)/Real(nlev))*ztop;
    }
    AtmosphericConditions conds(Tv0, Gammav, w0, ztop, tperiod, qv0, qv1);
    std::cout << conds.info_string();

    HostDynamics zdyn(nlev);
    zdyn.init_from_interface_heights(z_vals, conds);
    std::cout << zdyn.info_string();
    hbtest.run_test(zdyn, conds);
    hypsotest.run_test(zdyn, conds, 0.024);

    REQUIRE(hbtest.nerr == 0);
    REQUIRE(hypsotest.nerr == 0);

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
    const int nlev = 20;
    const Real ztop = 20E3;
    const AtmosphericConditions conds(Tv0, Gammav, w0, ztop, tperiod, qv0, qv1);
    std::cout << conds.info_string();

    HostDynamics pdyn(nlev);
    pdyn.init_from_uniform_pressures(conds);
    std::cout << pdyn.info_string();

    UniformThicknessPressureTest onedp(-(conds.pref-conds.ptop)/nlev);

    hbtest.run_test(pdyn, conds);
    onedp.run_test(pdyn);
    hypsotest.run_test(pdyn,conds, 0.07);

    REQUIRE(hbtest.nerr == 0);
    REQUIRE(onedp.nerr == 0);
    REQUIRE(hypsotest.nerr == 0);

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
    const Real ztop = height_at_pressure(p_vals[0], AtmosphericConditions::pref, Tv0, Gammav);
    AtmosphericConditions conds(Tv0, Gammav, w0, ztop, tperiod, qv0, qv1);
    std::cout << conds.info_string();
    HostDynamics pdyn(nlev);
    pdyn.init_from_interface_pressures(p_vals, conds);
    std::cout << pdyn.info_string();
    hbtest.run_test(pdyn, conds);
    hypsotest.run_test(pdyn,conds,0.27);

    REQUIRE(hbtest.nerr == 0);
    REQUIRE(hypsotest.nerr == 0);

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

TEST_CASE("vertical_convergence_dynamics_init", "[convergence]") {
  VerticalConvergenceTests convtests(10,8);
  const Real Tv0 = 300;
  const Real Gammav = 0.01;
  const Real w0 = 1;
  const Real tperiod = 900;
  const Real qv0 = 0.015;
  const Real qv1 = 2.5E-3;
  const Real ztop = 20E3;
  const AtmosphericConditions conds(Tv0, Gammav, w0, ztop, tperiod, qv0, qv1);

  SECTION("uniform dz tests") {
    for (int i=0; i<convtests.ntests; ++i) {
      const int nlev = convtests.nlevs[i];
      HostDynamics zdyn(nlev);
      zdyn.init_from_uniform_heights(conds);
      convtests.hydro_maxres[i] = convtests.run_hydrostatic_test(zdyn, conds);
      convtests.hypso_maxres[i] = convtests.run_hypsometric_test(zdyn, conds);
      convtests.ps_res[i] = convtests.run_ps_test(zdyn, conds);
      convtests.ztop_res[i] = convtests.run_ztop_test(zdyn, conds);
    }

    convtests.compute_appx_conv_rates();
    std::cout << convtests.info_string();
  }

  SECTION("uniform dp tests") {

  }
}

void HydrostaticBalanceTest::run_test(const HostDynamics& dyn, const AtmosphericConditions& ac,
  const Real tol) {
  using namespace constants;
  nerr = 0;

  const auto phi = ekat::scalarize(dyn.phi);
  const auto p = ekat::scalarize(dyn.p);

  Kokkos::parallel_reduce("HydrostaticBalanceTest::run_test", dyn.nlev(),
    KOKKOS_LAMBDA (const int k, int& errct) {
      // Taylor et al. 2020 fig. 1 level idx = k+1

      const int kphalf_idx = k+1; // array idx of interface k + 1/2
      const int kmhalf_idx = k; // array idx of interface k - 1/2

      const Real phimid = 0.5*(phi(kmhalf_idx) + phi(kphalf_idx));
      const Real zmid = phimid/gravity;
      const Real pres = hydrostatic_pressure_at_height(zmid, ac);

      if (!FloatingPoint<Real>::zero(pres - p(k), tol)) {
        ++errct;
        printf("hydrostatic test level %d: p = %f, expected = %f; |diff| = %18.15g\n", k, p(k), pres, std::abs(pres-p(k))/p(k));
      }
    }, nerr);

    if (nerr == 0) {
      std::cout << "Hydrostatic balance test passed with tolerance = " << tol << "\n";
    }
    else {
      std::cout << "Hydrostatic balance test failed\n";
    }

}

void UniformThicknessHeightTest::run_test(const HostDynamics& dyn, const Real tol) {
  const auto dz = ekat::scalarize(dyn.dz);
  nerr = 0;
  Kokkos::parallel_reduce("UniformThicknessHeightTest::run_test", dyn.nlev(),
    KOKKOS_LAMBDA (const int k, int& errct) {
      if (!FloatingPoint<Real>::equiv(dz(k), dzval, tol)) {
        printf("unif. dz level %d, dz = %f; expected %f; |diff| = %18.15g\n", k, dz(k), dzval, std::abs(dz(k)-dzval));
        ++errct;
      }
    }, nerr);
  if (nerr == 0) {
    std::cout << "Uniform height thickness test passed with tolerance " << tol << "\n";
  }
  else {
    std::cout << "Uniform height thickness test failed.\n";
  }
}

void UniformThicknessPressureTest::run_test(const HostDynamics& dyn, const Real tol) {
  const auto dp = ekat::scalarize(dyn.dp);
  nerr = 0;
  const int nlev = dyn.nlev();
  Kokkos::parallel_reduce("UniformThicknessPressureTest::run_test", nlev+1,
    KOKKOS_LAMBDA (const int k, int& errct) {
      if (!FloatingPoint<Real>::zero((dp(k)-dpval)/AtmosphericConditions::pref, tol)) {
        printf("unif. dp level %d, dp = %f; expected %f; rel. |diff| = %18.15g\n", k, dp(k), dpval, std::abs(dp(k)-dpval)/dp(k));
        ++errct;
      }
    },
  nerr);
  if (nerr == 0) {
    std::cout << "Uniform pressure thickness test passed with tolerance " << tol << "\n";
  }
  else {
    std::cout << "Uniform pressure thickness test failed.\n";
  }
}

void HypsometricLevelsTest::run_test(const HostDynamics& dyn, const AtmosphericConditions& ac, const Real tol) {
  using namespace constants;
  nerr = 0;
  const int nlev = dyn.nlev();
  const auto dz = ekat::scalarize(dyn.dz);
  const auto thetav = ekat::scalarize(dyn.thetav);
  const auto p = ekat::scalarize(dyn.p);
  const auto dp = ekat::scalarize(dyn.dp);
  Kokkos::parallel_reduce("HypsometricLevelsTest::run_test", nlev,
    KOKKOS_LAMBDA (const int k, int& errct) {
      if (k > 0 && k < nlev-1) {
      const Real Tv = thetav(k)*exner_function(p(k));
      const Real dpavg = 0.5*(dp(k) + dp(k+1));
      const Real p1 = p(k) - 0.5*dpavg; // (k>0 ? p(k) - 0.5*dpavg : ac.ptop);
      const Real p2 = p(k) + 0.5*dpavg; // (k<nlev-1 ? p(k) + 0.5*dpavg : AtmosphericConditions::pref);
      const Real dzhypso = r_gas_dry_air * Tv * std::log(p1/p2) / gravity;
      if (!FloatingPoint<Real>::zero( (dzhypso-dz(k))/dz(k),tol)) {
        printf("at level %d: dz = %f, dzhypso = %f; reldiff = %f\n",k,dz(k), dzhypso, std::abs(dz(k)-dzhypso)/dz(k));
        ++errct;
      }
      }
    }, nerr);

    if (nerr == 0) {
      std::cout << "Hypsometric layer thickness test passed with tolerance " << tol << "\n";
    }
    else {
      std::cout << "Hysometric layer thickness test failed.\n";
    }
}

Real VerticalConvergenceTests::run_hypsometric_test(
  const HostDynamics& dyn, const AtmosphericConditions& ac) {
  using namespace constants;

  const auto dz = ekat::scalarize(dyn.dz);
  const auto thetav = ekat::scalarize(dyn.thetav);
  const auto p = ekat::scalarize(dyn.p);
  const auto dp = ekat::scalarize(dyn.dp);
  Real maxres = 0;
  Kokkos::parallel_reduce("VerticalConvergenceTests::run_hypsometric_test", dyn.nlev(),
    KOKKOS_LAMBDA (const int k, Real& res) {
      const Real Tv = thetav(k)*exner_function(p(k));
      const Real dpavg = 0.5*(dp(k) + dp(k+1));
      const Real p1 = p(k) - 0.5*dpavg; // (k>0 ? p(k) - 0.5*dpavg : ac.ptop);
      const Real p2 = p(k) + 0.5*dpavg; // (k<nlev-1 ? p(k) + 0.5*dpavg : AtmosphericConditions::pref);
      const Real dzhypso = r_gas_dry_air * Tv * std::log(p1/p2) / gravity;
      if ( std::abs(dzhypso - dz(k)) > res)
        res = std::abs(dzhypso - dz(k));
    }, Kokkos::Max<Real>(maxres));
  return maxres;
}

Real VerticalConvergenceTests::run_hydrostatic_test(
  const HostDynamics& dyn, const AtmosphericConditions& ac) {
  using namespace constants;
  const auto dp = ekat::scalarize(dyn.dp);
  const auto dz = ekat::scalarize(dyn.dz);
  const auto rho = ekat::scalarize(dyn.rho);
  Real maxres;
  Kokkos::parallel_reduce("VerticalConvergenceTests::run_hydrostatic_test", dyn.nlev(),
  KOKKOS_LAMBDA (const int k, Real& res) {
    const Real dpavg = 0.5*(dp(k) + dp(k+1));
    const Real dpdz = dpavg/dz(k);
    const Real rhog = rho(k)*gravity;
    if (std::abs(dpdz + rhog) > res)
      res = std::abs(dpdz + rhog);
  }, Kokkos::Max<Real>(maxres));
  return maxres;
}

Real VerticalConvergenceTests::run_ps_test(
  const HostDynamics& dyn, const AtmosphericConditions& ac) {

  Real ps = 0;
  const auto dp = ekat::scalarize(dyn.dp);
  Kokkos::parallel_reduce("VerticalConvergenceTests::run_ps_test", dyn.nlev(),
    KOKKOS_LAMBDA (const int k, Real& psum) {
      psum -= 0.5*(dp(k) + dp(k+1));
    }, ps);
  return std::abs(ps + ac.ptop - AtmosphericConditions::pref);
}

void VerticalConvergenceTests::compute_appx_conv_rates() {
  for (int i=1; i<nlevs.size(); ++i) {
    const Real run = std::log(nlevs[i-1]) - std::log(nlevs[i]);
    hypso_rate[i] = (std::log(hypso_maxres[i]) - std::log(hypso_maxres[i-1])) / run;
    hydro_rate[i] = (std::log(hydro_maxres[i]) - std::log(hydro_maxres[i-1])) / run;
    ps_rate[i] = (std::log(ps_res[i]) - std::log(ps_res[i-1])) / run;
    ztop_rate[i] = (std::log(ztop_res[i]) - std::log(ztop_res[i-1])) / run;
  }
}

Real VerticalConvergenceTests::run_ztop_test(const HostDynamics& dyn, const AtmosphericConditions& ac) {
  Real zt = 0;
  const auto dz = ekat::scalarize(dyn.dz);
  Kokkos::parallel_reduce("VerticalConvergenceTests::run_ztop_test", dyn.nlev(),
    KOKKOS_LAMBDA (const int k, Real& zsum) {
      zsum += dz(k);
    }, zt);
  return std::abs(ac.ztop - zt);
}

std::string VerticalConvergenceTests::info_string() const {
  std::ostringstream ss;
  ss << "Convergence test info:\n";
  ss << "\t" << std::setw(30) << "test name" << std::setw(8) << "nlev" << std::setw(20) << "residue"
             << std::setw(20) << "appx. rate\n";
  for (int i=0; i<ntests; ++i) {
    ss << "\t" << std::setw(30) << (i==0 ? "dp/dz + \\rho g = 0" : " ") << std::setw(8) << nlevs[i] << std:: setw(20)
       << std::scientific << hydro_maxres[i] << std::setw(20) << std::scientific << hydro_rate[i] << "\n";
  }
  for (int i=0; i<ntests; ++i) {
  ss << "\t" << std::setw(30) << (i==0 ? "z2-z1 = (R T_v)/g log(p2/p1)" : " ") << std::setw(8) << nlevs[i] << std::setw(20)
       << std::scientific << hypso_maxres[i] << std::setw(20) << std::scientific << hypso_rate[i] << "\n";
  }
  for (int i=0; i<ntests; ++i) {
    ss << "\t" << std::setw(30) << (i==0 ? "ps = \\sum \\delta p + ptop" : " ") << std::setw(8) << nlevs[i] << std::setw(20)
       << std::scientific << ps_res[i] << std::setw(20) << std::scientific << ps_rate[i] << "\n";
  }
  for (int i=0; i<ntests; ++i) {
    ss << "\t" << std::setw(30) << (i==0 ? "ztop = \\sum \\delta z" : " ") << std::setw(8) << nlevs[i] << std::setw(20)
       << std::scientific << ztop_res[i] << std::setw(20) << std::scientific << ztop_rate[i] << "\n";
  }
  return ss.str();
}

