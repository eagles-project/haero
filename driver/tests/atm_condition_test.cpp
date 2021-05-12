#include <iostream>

#include "catch2/catch.hpp"
#include "driver/host_params.hpp"
#include "haero/physical_constants.hpp"

using namespace haero;
using namespace haero::driver;

struct InitialThicknessTest {
  int nerrz;
  int nerrp;
  InitialThicknessTest() : nerrz(0), nerrp(0) {}
  void run_test(const AtmosphericConditions& ac,
                const Real tol = FloatingPoint<Real>::zero_tol);
};

TEST_CASE("atmosphere_conditions", "") {
  InitialThicknessTest ittest;
  SECTION("dynamics") {
    const Real T0 = 300;
    const Real Gamma = 0.005;
    const Real w0 = 1;
    const Real ztop = 10E3;
    const Real tperiod = 900;
    const Real qv0 = 0.015;
    const Real qv1 = 1E-3;

    const auto conds =
        AtmosphericConditions(T0, Gamma, w0, ztop, tperiod, qv0, qv1);

    std::cout << conds.info_string();

    const Real z1000 = 1000;
    const Real qv1000 = water_vapor_mixing_ratio(z1000, conds);
    const Real Tv1000 = virtual_temperature(z1000, conds);
    const Real T1000 = temperature_from_virtual_temperature(Tv1000, qv1000);

    REQUIRE(T1000 < Tv1000);

    const Real p1000 = hydrostatic_pressure_at_height(z1000, conds);
    const Real zp1000 = height_at_pressure(p1000, conds);

    std::cout << "qv1000 = " << qv1000 << "\n";
    std::cout << "Tv1000 = " << Tv1000 << "\n";
    std::cout << "T1000  = " << T1000 << "\n";
    std::cout << "p1000  = " << p1000 << "\n";
    std::cout << "zp1000 = " << zp1000 << "\n";

    std::cout << "abs(zp1000 - z1000) = " << std::abs(z1000 - zp1000) << "\n";

    REQUIRE(FloatingPoint<Real>::equiv(z1000, zp1000,
                                       z1000 * FloatingPoint<Real>::zero_tol));

    ittest.run_test(conds);
    REQUIRE(ittest.nerrz == 0);
    REQUIRE(ittest.nerrp == 0);
  }
}

void InitialThicknessTest::run_test(const AtmosphericConditions& ac,
                                    const Real tol) {
  using namespace constants;
  nerrz = 0;
  nerrp = 0;
  const int nlev = 100;
  Kokkos::View<Real*> zvals("zvals", nlev + 1);
  Kokkos::View<Real*> pvals("pvals", nlev + 1);
  const Real dz = ac.ztop / nlev;
  Kokkos::parallel_for(
      "initialize uniform in z", nlev + 1, KOKKOS_LAMBDA(const int k) {
        zvals(k) = ac.ztop - k * dz;
        pvals(k) = hydrostatic_pressure_at_height(zvals(k), ac);
      });
  Kokkos::parallel_reduce(
      "thickness test 1", nlev,
      KOKKOS_LAMBDA(const int k, int& errct) {
        const Real pratio = pvals(k + 1) / pvals(k);
        const Real tv2 = virtual_temperature(zvals(k), ac);
        const Real tv1 = virtual_temperature(zvals(k + 1), ac);
        const Real rhs =
            std::pow(tv2 / tv1, -gravity / (r_gas_dry_air * ac.Gammav));
        if (!FloatingPoint<Real>::equiv(pratio, rhs, tol)) {
          ++errct;
          printf("ztest: at index k = %d: pratio = %f rhs = %f; |diff| = %f\n",
                 k, pratio, rhs, std::abs(pratio - rhs));
        }
      },
      nerrz);

  if (nerrz == 0) {
    std::cout << "uniform z thickness test passed with tolerance = " << tol
              << "\n";
  } else {
    std::cout << "uniform z thickness test failed.\n";
  }

  const Real dp = (AtmosphericConditions::pref - ac.ptop) / nlev;
  Kokkos::parallel_for(
      "initialize uniform in p", nlev + 1, KOKKOS_LAMBDA(const int k) {
        pvals(k) = ac.ptop + k * dp;
        zvals(k) = height_at_pressure(pvals(k), ac);
      });
  Kokkos::parallel_reduce(
      "thickness test 2", nlev,
      KOKKOS_LAMBDA(const int k, int& errct) {
        const Real pratio = pvals(k + 1) / pvals(k);
        const Real tv2 = virtual_temperature(zvals(k), ac);
        const Real tv1 = virtual_temperature(zvals(k + 1), ac);
        const Real rhs =
            std::pow(tv2 / tv1, -gravity / (r_gas_dry_air * ac.Gammav));
        if (!FloatingPoint<Real>::equiv(pratio, rhs, tol)) {
          ++errct;
          printf("ptest: at index k = %d: pratio = %f rhs = %f; |diff| = %f\n",
                 k, pratio, rhs, std::abs(pratio - rhs));
        }
      },
      nerrp);

  if (nerrz == 0) {
    std::cout << "uniform p thickness test passed with tolerance = " << tol
              << "\n";
  } else {
    std::cout << "uniform p thickness test failed.\n";
  }
}
