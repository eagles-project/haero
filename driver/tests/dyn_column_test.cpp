#include "driver/dyn_column.hpp"
#include "haero/floating_point.hpp"
#include <iostream>

#include "catch2/catch.hpp"

using namespace haero;
using namespace haero::driver;

TEST_CASE("dynamics_column", "") {
  const Real p0 = 100000;
  const Real T0 = 300;
  const Real Gamma = 0.001;
  const Real qv0 = 0.015;
  const Real qv1 = 5E-4;
  const auto conds = hydrostatic_conditions(p0, T0, Gamma, qv0, qv1);

  const int ncol = 1;

  SECTION("height_init") {
    /// In actual examples, these would come from an input yaml file
    const std::vector<Real> z_vals = {10000, 9000, 8000, 7000,
                                        6000, 5000, 4000, 3000,
                                        2000, 1500, 1000,  750,
                                         500,  400,  300,  200,
                                         100,   80,   60,   40,
                                          20,   10, 0};
    /// Assume input is interfaces
    const int nlev = z_vals.size()-1;

    DynColumn mycol(ncol, nlev);
    mycol.init_from_interface_heights(z_vals, conds);
    mycol.sum_psurf();

    auto hps = Kokkos::create_mirror_view(mycol.psurf);
    Kokkos::deep_copy(hps, mycol.psurf);
    int nerr = 0;
    for (int i=0; i<ncol; ++i) {
      if (!FloatingPoint<Real>::equiv(hps(i), p0)) {
        ++nerr;
        std::cout << "surface pressure error: ps " << hps(i) << " != p0 " << p0 << '\n';
      }
    }
    CHECK (nerr == 0);

    std::cout << mycol.info_string();
    auto writer = mycol.write_new_ncdata("dyn_column_test_zinit.nc", conds);
    mycol.update_ncdata(writer, 0);
    writer.close();
  }

  SECTION("pressure_init") {
    /// In actual examples, these would come from an input yaml file
    const std::vector<Real> p_vals = {20000, 30000, 40000, 50000, 60000, 70000, 85000, 92500, 100000};
    /// Assume input is interfaces
    const int nlev = p_vals.size()-1;

    DynColumn mycol(ncol, nlev);
    mycol.init_from_interface_pressures(p_vals, conds);
    mycol.sum_psurf();

    auto hps = Kokkos::create_mirror_view(mycol.psurf);
    Kokkos::deep_copy(hps, mycol.psurf);
    int nerr = 0;
    for (int i=0; i<ncol; ++i) {
      if (!FloatingPoint<Real>::equiv(hps(i), p0)) {
        ++nerr;
        std::cout << "surface pressure error: ps " << hps(i) << " != p0 " << p0 << '\n';
      }
    }
    CHECK (nerr == 0);

    Real pi_surf = 0;
    const auto dpids_loc = mycol.dpids;
    const auto ds_lev_loc = mycol.level_ds;
    Kokkos::parallel_reduce(nlev, KOKKOS_LAMBDA (const int i, Real& ps) {
      const int pack_idx = pack_info::pack_idx(i);
      const int vec_idx = pack_info::vec_idx(i);
      ps += dpids_loc(0,pack_idx)[vec_idx] * ds_lev_loc(pack_idx)[vec_idx];
    }, pi_surf);
    pi_surf += p_vals[0];
    CHECK(FloatingPoint<Real>::equiv(pi_surf, p_vals[nlev]));

    std::cout << mycol.info_string();
    auto writer = mycol.write_new_ncdata("dyn_column_test_pinit.nc", conds);

    mycol.update_ncdata(writer, 0);
    writer.close();
  }
}
