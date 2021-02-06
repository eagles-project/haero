#include "driver/host_state.hpp"
#include "driver/host_params.hpp"
#include "catch2/catch.hpp"
#include <iostream>

using namespace haero;
using namespace driver;

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