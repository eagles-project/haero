#include "driver/atmosphere.hpp"
#include "catch2/catch.hpp"
#include "haero/physical_constants.hpp"
#include <iostream>

using namespace haero;

TEST_CASE("atmosphere_conditions", "") {
  SECTION("uniform") {
    std::cout << "uniform\n";
    AtmosphericConditions conds;
    conds.params.uniform.mu = molec_weight_dry_air_g_per_mole;
    conds.params.uniform.H = 20.0;
    conds.params.uniform.p0 = p_std_pa;
    conds.params.uniform.T0 = t_freeze_h2o_k;
    conds.params.uniform.phi0 = 0.25;
    conds.params.uniform.N0 = 0.1;

    std::cout << "\tmu   = " << conds.params.uniform.mu << '\n'
              << "\tH    = " << conds.params.uniform.H << '\n'
              << "\tp0   = " << conds.params.uniform.p0 << '\n'
              << "\tT0   = " << conds.params.uniform.T0 << '\n'
              << "\tphi0 = " << conds.params.uniform.phi0 << '\n'
              << "\tN0   = " << conds.params.uniform.N0 << '\n';
  }

  SECTION("hydrostatic") {
    const Real p0 = 100000;
    const Real T0 = 300;
    const Real Gamma = 0.005;
    const Real qv0 = 0.015;
    const Real qv1 = 1E-3;

    REQUIRE(FloatingPoint<Real>::in_bounds(p0, 100000, p_std_pa));

    const auto conds = hydrostatic_conditions(p0, T0, Gamma, qv0, qv1);

    const Real z1000 = 1000;
    const Real qv1000 = water_vapor_mixing_ratio(z1000, qv0, qv1);
    const Real Tv1000 = virtual_temperature(z1000, T0, Gamma);
    const Real T1000 = temperature_from_virtual_temperature(Tv1000, qv1000);

    REQUIRE(T1000 < Tv1000);

    const Real p1000 = hydrostatic_pressure_at_height(z1000, p0, T0, Gamma);
    const Real zp1000 = height_at_pressure(p1000, p0, T0, Gamma);

    std::cout << "qv1000 = " << qv1000 << "\n";
    std::cout << "Tv1000 = " << Tv1000 << "\n";
    std::cout << "T1000  = " << T1000 << "\n";
    std::cout << "p1000  = " << p1000 << "\n";
    std::cout << "zp1000 = " << zp1000 << "\n";

    std::cout << "abs(zp1000 - z1000) = " << std::abs(z1000 - zp1000) << "\n";

    REQUIRE(FloatingPoint<Real>::equiv(z1000, zp1000, z1000*FloatingPoint<Real>::zero_tol));


  }
}
