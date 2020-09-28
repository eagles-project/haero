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
}
