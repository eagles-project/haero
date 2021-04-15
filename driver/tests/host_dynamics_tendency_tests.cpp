#include "driver/host_dynamics_tends.hpp"
#include "driver/host_dynamics.hpp"
#include "driver/host_params.hpp"
#include "haero/atmosphere.hpp"
#include "catch2/catch.hpp"
#include <iostream>

using namespace haero;
using namespace driver;

TEST_CASE("basic tendency unit tests", "") {
  const int nlev = 10;
  const Real t = 0;
  const auto conds = AtmosphericConditions();

  HostDynamics dyn(nlev);
  dyn.init_from_uniform_pressures(conds);

  DynamicsTendencies dtends(nlev);
  dtends.compute(t, dyn.phi, dyn.rho, conds);
}
