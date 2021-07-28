#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <vector>

#include "catch2/catch.hpp"
#include "driver/host_params.hpp"
#include "driver/ncwriter.hpp"
#include "driver/ncwriter_impl.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "haero/atmosphere.hpp"
#include "haero/math.hpp"

#define private \
 public  // unit tests need to access HostDynamics private variables
#include "driver/host_dynamics.hpp"
#undef private  // back to normal

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
  UniformThicknessHeightTest(const Real layer_depth)
      : nerr(0), dzval(layer_depth) {}
  void run_test(const HostDynamics& dyn,
                const Real tol = FloatingPoint<Real>::zero_tol);
};

struct HypsometricLevelsTest {
  int nerr;
  HypsometricLevelsTest() : nerr(0) {}
  void run_test(const HostDynamics& dyn, const AtmosphericConditions& ac,
                const Real tol = FloatingPoint<Real>::zero_tol);
};

struct VerticalConvergenceTests {
  int nlevstart;
  int ntests;
  std::vector<Real> hypso_maxres;
  std::vector<Real> hypso_avgres;
  std::vector<Real> hypso_avgrate;
  std::vector<Real> hypso_rate;
  std::vector<Real> hydro_maxres;
  std::vector<Real> hydro_avgres;
  std::vector<Real> hydro_rate;
  std::vector<Real> hydro_avgrate;
  std::vector<Real> ps_res;
  std::vector<Real> ps_rate;
  std::vector<Real> ztop_res;
  std::vector<Real> ztop_rate;
  std::vector<int> nlevs;

  Real avg_rate_hydro_max;
  Real avg_rate_hydro_avg;
  Real avg_rate_hypso_max;
  Real avg_rate_hypso_avg;
  Real max_ps_err;
  Real max_ztop_err;

  VerticalConvergenceTests(const int nlev0 = 10, const int nt = 5)
      : nlevstart(nlev0),
        ntests(nt),
        hypso_maxres(nt, 0),
        hypso_avgres(nt, 0),
        hypso_avgrate(nt, 0),
        hypso_rate(nt, 0),
        hydro_maxres(nt, 0),
        hydro_avgres(nt, 0),
        hydro_rate(nt, 0),
        hydro_avgrate(nt, 0),
        ps_res(nt, 0),
        ps_rate(nt, 0),
        ztop_res(nt, 0),
        ztop_rate(nt, 0),
        nlevs(nt) {
    EKAT_ASSERT(nlev0 > 0);
    EKAT_ASSERT(nt > 0);
    for (int i = 0; i < nt; ++i) {
      nlevs[i] = std::pow(2, i) * nlevstart;
    }
  }

  VerticalConvergenceTests(const std::vector<int>& levels)
      : nlevstart(levels[0]),
        ntests(levels.size()),
        hypso_maxres(levels.size(), 0),
        hypso_avgres(levels.size(), 0),
        hypso_avgrate(levels.size(), 0),
        hypso_rate(levels.size(), 0),
        hydro_maxres(levels.size(), 0),
        hydro_avgres(levels.size(), 0),
        hydro_rate(levels.size(), 0),
        hydro_avgrate(levels.size(), 0),
        ps_res(levels.size(), 0),
        ps_rate(levels.size(), 0),
        ztop_res(levels.size(), 0),
        ztop_rate(levels.size(), 0),
        nlevs(levels) {}

  void run_hypsometric_test(const int test_idx, const HostDynamics& dyn,
                            const AtmosphericConditions& ac);
  void run_hydrostatic_test(const int test_idx, const HostDynamics& dyn,
                            const AtmosphericConditions& ac);
  void run_ps_test(const int test_idx, const HostDynamics& dyn,
                   const AtmosphericConditions& ac);
  void run_ztop_test(const int test_idx, const HostDynamics& dyn,
                     const AtmosphericConditions& ac);
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

  SECTION("z_unif_init") {
    std::cout << "Uniform height levels\n";
    const int nlev = 320;
    const Real ztop = 20E3;
    const AtmosphericConditions conds(Tv0, Gammav, w0, ztop, tperiod, qv0, qv1);
    std::cout << conds.info_string();

    UniformThicknessHeightTest onedz(ztop / nlev);

    HostDynamics zdyn(nlev);
    zdyn.init_from_uniform_heights(conds);
    std::cout << zdyn.info_string();
    hbtest.run_test(zdyn, conds, 4.9 * FloatingPoint<Real>::zero_tol);
    onedz.run_test(zdyn, (std::is_same<float, Real>::value ? 1.1e-3 : 2.1e-12));
    hypsotest.run_test(zdyn, conds, 1.5e-2);

    REQUIRE(hbtest.nerr == 0);
    REQUIRE(onedz.nerr == 0);
    REQUIRE(hypsotest.nerr == 0);

    std::cout << "\tinitialization unit tests complete\n";

    /// Create a new netcdf file
    const std::string fname = "host_dynamics_test_zinit_unif.nc";
    NcWriter writer(fname);
    writer.define_time_var();
    zdyn.nc_init_dynamics_variables(writer, conds);

    std::cout << "ncwriter initialized\n";

    size_t time_idx = 0;
    zdyn.nc_write_data(writer, time_idx);
    Kokkos::View<PackType*> temperature("temperature",
                                        PackInfo::num_packs(nlev));
    Kokkos::View<PackType*> rel_humidity("relative_humidity",
                                         PackInfo::num_packs(nlev));
    Kokkos::View<PackType*> level_heights("level_heights",
                                          PackInfo::num_packs(nlev + 1));
    auto atm =
        zdyn.create_atmospheric_state(temperature, rel_humidity, level_heights);
    writer.define_atm_state_vars(atm);
    writer.add_atm_state_data(atm, time_idx);

    Real t = 0.5 * conds.tperiod;
    ++time_idx;

    zdyn.update(t, conds);
    zdyn.update_atmospheric_state(atm);
    writer.add_time_value(t);
    zdyn.nc_write_data(writer, time_idx);
    writer.add_atm_state_data(atm, time_idx);

    writer.close();
  }

  SECTION("height init -- specified heights") {
    std::cout << "User-specified height levels\n";

    const int nlev = 20;
    /// In actual examples, these would come from an input yaml file
    std::vector<Real> z_vals(nlev + 1);
    const Real ztop = 20E3;
    for (int i = 0; i < nlev + 1; ++i) {
      z_vals[i] = square(Real(i) / Real(nlev)) * ztop;
    }
    AtmosphericConditions conds(Tv0, Gammav, w0, ztop, tperiod, qv0, qv1);
    std::cout << conds.info_string();

    HostDynamics zdyn(nlev);
    zdyn.init_from_interface_heights(z_vals, conds);
    std::cout << zdyn.info_string();
    hbtest.run_test(zdyn, conds, 3.4 * FloatingPoint<Real>::zero_tol);
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
    Kokkos::View<PackType*> temperature("temperature",
                                        PackInfo::num_packs(nlev));
    Kokkos::View<PackType*> rel_humidity("relative_humidity",
                                         PackInfo::num_packs(nlev));
    Kokkos::View<PackType*> level_heights("level_heights",
                                          PackInfo::num_packs(nlev + 1));
    auto atm =
        zdyn.create_atmospheric_state(temperature, rel_humidity, level_heights);
    writer.define_atm_state_vars(atm);
    writer.add_atm_state_data(atm, time_idx);

    Real t = 0.5 * conds.tperiod;
    ++time_idx;
    zdyn.update(t, conds);
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

    hbtest.run_test(pdyn, conds, 2.7 * FloatingPoint<Real>::zero_tol);
    hypsotest.run_test(pdyn, conds, 0.07);

    REQUIRE(hbtest.nerr == 0);
    REQUIRE(hypsotest.nerr == 0);

    /// Create a new netcdf file
    const std::string fname = "host_dynamics_test_pinit_unif.nc";

    NcWriter writer(fname);
    writer.define_time_var();
    pdyn.nc_init_dynamics_variables(writer, conds);

    size_t time_idx = 0;
    pdyn.nc_write_data(writer, time_idx);
    Kokkos::View<PackType*> temperature("temperature",
                                        PackInfo::num_packs(nlev));
    Kokkos::View<PackType*> rel_humidity("relative_humidity",
                                         PackInfo::num_packs(nlev));
    Kokkos::View<PackType*> level_heights("level_heights",
                                          PackInfo::num_packs(nlev + 1));
    auto atm =
        pdyn.create_atmospheric_state(temperature, rel_humidity, level_heights);
    writer.define_atm_state_vars(atm);
    writer.add_atm_state_data(atm, time_idx);

    Real t = 0.5 * conds.tperiod;
    ++time_idx;
    pdyn.update(t, conds);
    pdyn.update_atmospheric_state(atm);

    writer.add_time_value(t);
    pdyn.nc_write_data(writer, time_idx);
    writer.add_atm_state_data(atm, time_idx);

    writer.close();
  }

  SECTION("pressure_init -- specified pressures") {
    /// In actual examples, these would come from an input yaml file
    const std::vector<Real> p_vals = {20000, 30000, 40000, 50000, 60000,
                                      70000, 85000, 92500, 100000};
    /// Assume input is interfaces
    const int nlev = p_vals.size() - 1;
    const Real ztop =
        height_at_pressure(p_vals[0], AtmosphericConditions::pref, Tv0, Gammav);
    AtmosphericConditions conds(Tv0, Gammav, w0, ztop, tperiod, qv0, qv1);
    std::cout << conds.info_string();
    HostDynamics pdyn(nlev);
    pdyn.init_from_interface_pressures(p_vals, conds);
    std::cout << pdyn.info_string();
    hbtest.run_test(pdyn, conds, 3.5 * FloatingPoint<Real>::zero_tol);
    hypsotest.run_test(pdyn, conds, 0.27);

    REQUIRE(hbtest.nerr == 0);
    REQUIRE(hypsotest.nerr == 0);

    /// Create a new netcdf file
    const std::string fname = "host_dynamics_test_pinit.nc";
    NcWriter writer(fname);
    writer.define_time_var();
    pdyn.nc_init_dynamics_variables(writer, conds);

    size_t time_idx = 0;
    pdyn.nc_write_data(writer, time_idx);
    Kokkos::View<PackType*> temperature("temperature",
                                        PackInfo::num_packs(nlev));
    Kokkos::View<PackType*> rel_humidity("relative_humidity",
                                         PackInfo::num_packs(nlev));
    Kokkos::View<PackType*> level_heights("level_heights",
                                          PackInfo::num_packs(nlev + 1));
    auto atm =
        pdyn.create_atmospheric_state(temperature, rel_humidity, level_heights);
    writer.define_atm_state_vars(atm);
    writer.add_atm_state_data(atm, time_idx);

    Real t = 0.5 * conds.tperiod;
    ++time_idx;
    pdyn.update(t, conds);
    pdyn.update_atmospheric_state(atm);

    writer.add_time_value(t);
    pdyn.nc_write_data(writer, time_idx);
    writer.add_atm_state_data(atm, time_idx);

    writer.close();
  }
}

TEST_CASE("vertical_convergence_dynamics_init", "[convergence]") {
  const std::vector<int> nlevs = {10,  20,  40,  80,  120, 160, 180,
                                  200, 220, 240, 280, 320, 400, 640};
  //   VerticalConvergenceTests convtests(10,8);
  VerticalConvergenceTests convtests(nlevs);
  const Real Tv0 = 300;
  const Real Gammav = 0.01;
  const Real w0 = 1;
  const Real tperiod = 900;
  const Real qv0 = 0.015;
  const Real qv1 = 2.5E-3;
  const Real ztop = 20E3;
  const AtmosphericConditions conds(Tv0, Gammav, w0, ztop, tperiod, qv0, qv1);

  SECTION("uniform dz tests") {
    for (int i = 0; i < convtests.ntests; ++i) {
      const int nlev = convtests.nlevs[i];
      HostDynamics zdyn(nlev);
      zdyn.init_from_uniform_heights(conds);
      convtests.run_hydrostatic_test(i, zdyn, conds);
      convtests.run_hypsometric_test(i, zdyn, conds);
      convtests.run_ps_test(i, zdyn, conds);
      convtests.run_ztop_test(i, zdyn, conds);
    }

    convtests.compute_appx_conv_rates();
    std::cout << convtests.info_string();

    // ptop + sum of all levels' pressure must equal surface pressure
    REQUIRE(FloatingPoint<Real>::zero(
        convtests.max_ps_err,
        (std::is_same<float, Real>::value ? 1.6e-2
                                          : FloatingPoint<Real>::zero_tol)));
    // sum of level thicknesses must equal ztop
    REQUIRE(FloatingPoint<Real>::zero(
        convtests.max_ztop_err,
        (std::is_same<float, Real>::value ? 5.3e-2 : 5.1e-11)));

#if HAERO_DOUBLE_PRECISION
    // rate of average error in hydrostatic equation should converge at 2nd
    // order
    REQUIRE(FloatingPoint<Real>::equiv(convtests.avg_rate_hydro_max, 2, 0.01));
    // rate of max error in hydrostatic equation should converge at 2nd order
    REQUIRE(FloatingPoint<Real>::equiv(convtests.avg_rate_hydro_avg, 2, 0.01));
    // rate of average error in hypsometric equation should converge at 3rd
    // order
    REQUIRE(FloatingPoint<Real>::equiv(convtests.avg_rate_hypso_max, 3, 0.05));
    // rate of max error in hypsometric equation should converge at 3rd order
    REQUIRE(FloatingPoint<Real>::equiv(convtests.avg_rate_hypso_avg, 3, 0.01));
#endif
  }

  //   SECTION("uniform dp tests") {
  //
  //   }
}

void HydrostaticBalanceTest::run_test(const HostDynamics& dyn,
                                      const AtmosphericConditions& ac,
                                      const Real tol) {
  nerr = 0;

  const auto phi = ekat::scalarize(dyn.phi);
  const auto p = ekat::scalarize(dyn.p);

  Kokkos::parallel_reduce(
      "HydrostaticBalanceTest::run_test", dyn.nlev(),
      KOKKOS_LAMBDA(const int k, int& errct) {
        // Taylor et al. 2020 fig. 1 level idx = k+1

        const int kphalf_idx = k + 1;  // array idx of interface k + 1/2
        const int kmhalf_idx = k;      // array idx of interface k - 1/2

        const Real phimid = 0.5 * (phi(kmhalf_idx) + phi(kphalf_idx));
        const Real zmid = phimid / Constants::gravity;
        const Real pres = hydrostatic_pressure_at_height(zmid, ac);

        if (!FloatingPoint<Real>::zero((pres - p(k)) / p(k), tol)) {
          ++errct;
          printf(
              "hydrostatic test level %d: p = %f, expected = %f; |diff| = "
              "%18.15g, tol = %g\n",
              k, p(k), pres, std::abs(pres - p(k)) / p(k), tol);
        }
      },
      nerr);

  if (nerr == 0) {
    std::cout << "Hydrostatic balance test passed with tolerance = " << tol
              << "\n";
  } else {
    std::cout << "Hydrostatic balance test failed\n";
  }
}

void UniformThicknessHeightTest::run_test(const HostDynamics& dyn,
                                          const Real tol) {
  const auto dz = ekat::scalarize(dyn.dz);
  nerr = 0;
  const Real dzval_ = dzval;
  Kokkos::parallel_reduce(
      "UniformThicknessHeightTest::run_test", dyn.nlev(),
      KOKKOS_LAMBDA(const int k, int& errct) {
        if (!FloatingPoint<Real>::equiv(dz(k), dzval_, tol)) {
          printf("unif. dz level %d, dz = %f; expected %f; |diff| = %18.15g\n",
                 k, dz(k), dzval_, std::abs(dz(k) - dzval_));
          ++errct;
        }
      },
      nerr);
  if (nerr == 0) {
    std::cout << "Uniform height thickness test passed with tolerance " << tol
              << "\n";
  } else {
    std::cout << "Uniform height thickness test failed with tolerance " << tol
              << "\n";
  }
}

void HypsometricLevelsTest::run_test(const HostDynamics& dyn,
                                     const AtmosphericConditions& ac,
                                     const Real tol) {
  nerr = 0;
  const int nlev = dyn.nlev();
  const auto dz = ekat::scalarize(dyn.dz);
  const auto thetav = ekat::scalarize(dyn.thetav);
  const auto p = ekat::scalarize(dyn.p);
  const auto pint = ekat::scalarize(dyn.phydro_int);
  Kokkos::parallel_reduce(
      "HypsometricLevelsTest::run_test", nlev,
      KOKKOS_LAMBDA(const int k, int& errct) {
        if (k > 0 && k < nlev - 1) {
          const Real Tv = thetav(k) * exner_function(p(k));
          const Real p1 = pint(k);
          const Real p2 = pint(k + 1);
          const Real dzhypso =
              -Constants::r_gas_dry_air * Tv * std::log(p1 / p2) / Constants::gravity;
          if (!FloatingPoint<Real>::zero((dzhypso - dz(k)) / dz(k), tol)) {
            printf("at level %d: dz = %f, dzhypso = %f; reldiff = %f\n", k,
                   dz(k), dzhypso, std::abs(dz(k) - dzhypso) / dz(k));
            ++errct;
          }
        }
      },
      nerr);

  if (nerr == 0) {
    std::cout << "Hypsometric layer thickness test passed with tolerance "
              << tol << "\n";
  } else {
    std::cout << "Hysometric layer thickness test failed.\n";
  }
}

void VerticalConvergenceTests::run_hypsometric_test(
    const int test_idx, const HostDynamics& dyn,
    const AtmosphericConditions& ac) {

  const auto dz = ekat::scalarize(dyn.dz);
  const auto thetav = ekat::scalarize(dyn.thetav);
  const auto p = ekat::scalarize(dyn.p);
  const auto pint = ekat::scalarize(dyn.phydro_int);

  Real maxres = 0;
  Kokkos::parallel_reduce(
      "VerticalConvergenceTests::run_hypsometric_test", dyn.nlev(),
      KOKKOS_LAMBDA(const int k, Real& res) {
        const Real Tv = thetav(k) * exner_function(p(k));
        const Real p1 = pint(k);
        const Real p2 = pint(k + 1);
        const Real dzhypso = -Constants::r_gas_dry_air * Tv * std::log(p1 / p2) / Constants::gravity;
        if (std::abs(dzhypso - dz(k)) > res) res = std::abs(dzhypso - dz(k));
      },
      Kokkos::Max<Real>(maxres));
  hypso_maxres[test_idx] = maxres;

  Real rsum = 0;
  Kokkos::parallel_reduce(
      "VerticalConvergenceTests::run_hypsometric_avg_test", dyn.nlev(),
      KOKKOS_LAMBDA(const int k, Real& ressum) {
        const Real Tv = thetav(k) * exner_function(p(k));
        const Real p1 = pint(k);
        const Real p2 = pint(k + 1);
        const Real dzhypso = -Constants::r_gas_dry_air * Tv * std::log(p1 / p2) / Constants::gravity;
        ressum += std::abs(dzhypso - dz(k));
      },
      rsum);
  hypso_avgres[test_idx] = rsum / dyn.nlev();
}

void VerticalConvergenceTests::run_hydrostatic_test(
    const int test_idx, const HostDynamics& dyn,
    const AtmosphericConditions& ac) {

  const auto dp = ekat::scalarize(dyn.hydrostatic_dp);
  const auto dz = ekat::scalarize(dyn.dz);
  const auto rho = ekat::scalarize(dyn.rho);

  Real maxres = 0;
  Kokkos::parallel_reduce(
      "VerticalConvergenceTests::run_hydrostatic_test", dyn.nlev(),
      KOKKOS_LAMBDA(const int k, Real& res) {
        const Real dpdz = dp(k) / dz(k);
        const Real rhog = rho(k) * Constants::gravity;
        if (std::abs(dpdz + rhog) > res) res = std::abs(dpdz + rhog);
      },
      Kokkos::Max<Real>(maxres));
  hydro_maxres[test_idx] = maxres;

  Real rsum = 0;
  Kokkos::parallel_reduce(
      "VerticalConvergenceTests::run_hydrostatic_test_avg", dyn.nlev(),
      KOKKOS_LAMBDA(const int k, Real& ressum) {
        const Real dpdz = dp(k) / dz(k);
        const Real rhog = rho(k) * Constants::gravity;
        ressum += std::abs(dpdz + rhog);
      },
      rsum);
  hydro_avgres[test_idx] = rsum / dyn.nlev();
}

void VerticalConvergenceTests::run_ps_test(const int test_idx,
                                           const HostDynamics& dyn,
                                           const AtmosphericConditions& ac) {
  Real ps = 0;
  const auto dp = ekat::scalarize(dyn.hydrostatic_dp);
  Kokkos::parallel_reduce(
      "VerticalConvergenceTests::run_ps_test", dyn.nlev(),
      KOKKOS_LAMBDA(const int k, Real& psum) { psum -= dp(k); }, ps);
  ps_res[test_idx] = std::abs(ps + ac.ptop - AtmosphericConditions::pref);
}

void VerticalConvergenceTests::compute_appx_conv_rates() {
  for (int i = 1; i < nlevs.size(); ++i) {
    const Real run = std::log(nlevs[i - 1]) - std::log(nlevs[i]);
    hypso_rate[i] =
        (std::log(hypso_maxres[i]) - std::log(hypso_maxres[i - 1])) / run;
    hypso_avgrate[i] =
        (std::log(hypso_avgres[i]) - std::log(hypso_avgres[i - 1])) / run;
    hydro_rate[i] =
        (std::log(hydro_maxres[i]) - std::log(hydro_maxres[i - 1])) / run;
    hydro_avgrate[i] =
        (std::log(hydro_avgres[i]) - std::log(hydro_avgres[i - 1])) / run;
    ps_rate[i] = (std::log(ps_res[i]) - std::log(ps_res[i - 1])) / run;
    ztop_rate[i] = (std::log(ztop_res[i]) - std::log(ztop_res[i - 1])) / run;
  }

  max_ps_err = *std::max_element(ps_res.begin(), ps_res.end());
  max_ztop_err = *std::max_element(ztop_res.begin(), ztop_res.end());
  avg_rate_hydro_max =
      std::accumulate(hydro_rate.begin() + 1, hydro_rate.end(), 0.0) /
      (ntests - 1);
  avg_rate_hydro_avg =
      std::accumulate(hydro_avgrate.begin() + 1, hydro_avgrate.end(), 0.0) /
      (ntests - 1);
  avg_rate_hypso_max =
      std::accumulate(hypso_rate.begin() + 1, hypso_rate.end(), 0.0) /
      (ntests - 1);
  avg_rate_hypso_avg =
      std::accumulate(hypso_avgrate.begin() + 1, hypso_avgrate.end(), 0.0) /
      (ntests - 1);
}

void VerticalConvergenceTests::run_ztop_test(const int test_idx,
                                             const HostDynamics& dyn,
                                             const AtmosphericConditions& ac) {
  Real zt = 0;
  const auto dz = ekat::scalarize(dyn.dz);
  Kokkos::parallel_reduce(
      "VerticalConvergenceTests::run_ztop_test", dyn.nlev(),
      KOKKOS_LAMBDA(const int k, Real& zsum) { zsum += dz(k); }, zt);
  ztop_res[test_idx] = std::abs(ac.ztop - zt);
}

std::string VerticalConvergenceTests::info_string() const {
  std::ostringstream ss;
  ss << "Convergence test info:\n";
  ss << "\t" << std::setw(35) << "test name" << std::setw(8) << "nlev"
     << std::setw(20) << "residue" << std::setw(20) << "appx. rate\n";
  for (int i = 0; i < ntests; ++i) {
    ss << "\t" << std::setw(35) << (i == 0 ? "dp/dz + \\rho g = 0 [MAX]" : " ")
       << std::setw(8) << nlevs[i] << std::setw(20) << std::scientific
       << hydro_maxres[i] << std::setw(20) << std::scientific << hydro_rate[i]
       << "\n";
  }
  for (int i = 0; i < ntests; ++i) {
    ss << "\t" << std::setw(35) << (i == 0 ? "dp/dz + \\rho g = 0 [AVG]" : " ")
       << std::setw(8) << nlevs[i] << std::setw(20) << std::scientific
       << hydro_avgres[i] << std::setw(20) << std::scientific
       << hydro_avgrate[i] << "\n";
  }
  for (int i = 0; i < ntests; ++i) {
    ss << "\t" << std::setw(35)
       << (i == 0 ? "z2-z1 = (R T_v)/g log(p2/p1) [MAX]" : " ") << std::setw(8)
       << nlevs[i] << std::setw(20) << std::scientific << hypso_maxres[i]
       << std::setw(20) << std::scientific << hypso_rate[i] << "\n";
  }
  for (int i = 0; i < ntests; ++i) {
    ss << "\t" << std::setw(35)
       << (i == 0 ? "z2-z1 = (R T_v)/g log(p2/p1) [AVG]" : " ") << std::setw(8)
       << nlevs[i] << std::setw(20) << std::scientific << hypso_avgres[i]
       << std::setw(20) << std::scientific << hypso_avgrate[i] << "\n";
  }
  for (int i = 0; i < ntests; ++i) {
    ss << "\t" << std::setw(35)
       << (i == 0 ? "ps = \\sum \\delta p + ptop" : " ") << std::setw(8)
       << nlevs[i] << std::setw(20) << std::scientific << ps_res[i]
       << std::setw(20) << std::scientific << ps_rate[i] << "\n";
  }
  for (int i = 0; i < ntests; ++i) {
    ss << "\t" << std::setw(35) << (i == 0 ? "ztop = \\sum \\delta z" : " ")
       << std::setw(8) << nlevs[i] << std::setw(20) << std::scientific
       << ztop_res[i] << std::setw(20) << std::scientific << ztop_rate[i]
       << "\n";
  }
  ss << "-------------------\n";
  ss << "max. PS err = " << max_ps_err << "\n";
  ss << "max. ztop err = " << max_ztop_err << "\n";
  ss << "conv. rate for hydro. max. : " << avg_rate_hydro_max << "\n";
  ss << "conv. rate for hydro. avg. : " << avg_rate_hydro_avg << "\n";
  ss << "conv. rate for hypso. max. : " << avg_rate_hypso_max << "\n";
  ss << "conv. rate for hypso. avg. : " << avg_rate_hypso_avg << "\n";
  return ss.str();
}
