#include "driver/host_params.hpp"
#include "catch2/catch.hpp"
#include "haero/physical_constants.hpp"
#include <iostream>

using namespace haero;
using namespace haero::driver;

TEST_CASE("atmosphere_conditions", "") {

  SECTION("dynamics") {
    const Real T0 = 300;
    const Real Gamma = 0.005;
    const Real w0 = 1;
    const Real ztop = 10E3;
    const Real tperiod = 900;
    const Real qv0 = 0.015;
    const Real qv1 = 1E-3;

    const auto conds = AtmosphericConditions(T0,Gamma,w0,ztop,tperiod,qv0,qv1);
    
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

    REQUIRE(FloatingPoint<Real>::equiv(z1000, zp1000, z1000*FloatingPoint<Real>::zero_tol));


  }
}
