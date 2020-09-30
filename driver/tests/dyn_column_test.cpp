#include "driver/dyn_column.hpp"
#include <iostream>

#include "catch2/catch.hpp"

using namespace haero;

TEST_CASE("dynamics_column", "") {
  const Real p0 = 100000;
  const Real T0 = 300;
  const Real Gamma = 0.001;
  const Real qv0 = 0.015;
  const Real qv1 = 5E-4;
  const auto conds = hydrostatic_conditions(p0, T0, Gamma, qv0, qv1);

  const int ncol = 1;

  SECTION("height_init") {
    const std::vector<Real> z_vals = {10000, 9000, 8000, 7000,
                                        6000, 5000, 4000, 3000,
                                        2000, 1500, 1000,  750,
                                         500,  400,  300,  200,
                                         100,   80,   60,   40,
                                          20,   10, 0};
    const int nlev = z_vals.size()-1;

    DynColumn mycol(ncol, nlev);
    mycol.init_from_interface_heights(z_vals, conds);

    std::cout << mycol.info_string();
    auto writer = mycol.write_new_ncdata("dyn_column_test.nc");
    mycol.update_ncdata(writer, 0);
    writer.close();
  }

  SECTION("pressure_init") {
    const std::vector<Real> p_vals = {20000, 30000, 40000, 50000, 60000, 70000, 85000, 92500, 100000};
    const int nlev = p_vals.size()-1;

    DynColumn mycol(ncol, nlev);


  }
}
