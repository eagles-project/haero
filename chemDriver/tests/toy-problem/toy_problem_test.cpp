#include "catch2/catch.hpp"
#include "chemDriver/chemDriver.hpp"
#include "toy_problem.hpp"
#include "chemDriver/read_chem_input.hpp"
#include "haero/floating_point.hpp"

using namespace haero;
using namespace haero::chemDriver;

TEST_CASE("TChem tendency computation tests", "haero_unit_tests"){

  SimulationInput sim_input = read_chem_input("toy_problem_input.yml");
  // this creates the toy-problem-specific input file for TChem
      // FIXME: learn more about this and whether they've switched to the yaml input spec
  create_chem_files();

  SECTION("light side of terminator"){
    {
      ChemSolver chem_solver(sim_input);

      // run the problem on device
      real_type_2d_view results = chem_solver.get_results();

      // create mirror view and deep copy to host
      auto results_host = Kokkos::create_mirror_view(results);
      Kokkos::deep_copy(results_host, results);

      // eq (4) tendency should be positive, eq (5) negative
      const Real val00 = results_host(0, 0);
      const Real val01 = results_host(0, 1);
      REQUIRE(FloatingPoint<Real>::in_bounds(val00, 0.0, 1.0e14,
                                             FloatingPoint<Real>::zero_tol));
      REQUIRE(FloatingPoint<Real>::in_bounds(val01, -1.0e14, 0.0,
                                             FloatingPoint<Real>::zero_tol));
      // we ran for two batches, so be sure that the same inputs give the same
      // outputs for different batches
      const Real val10 = results_host(1, 0);
      const Real val11 = results_host(1, 1);
      REQUIRE(FloatingPoint<Real>::zero(val00 - val10,
                                        FloatingPoint<Real>::zero_tol));
      REQUIRE(FloatingPoint<Real>::zero(val01 - val11,
                                        FloatingPoint<Real>::zero_tol));
    }
  }

  SECTION("zero tendencies"){

    // calculate initial X2 concentration such that tendencies will be zero
    Real initX = sim_input.species[0].initial_value;
    Real k1 = sim_input.reactions[0].rate_coefficients["A"];
    Real k2 = sim_input.reactions[1].rate_coefficients["A"];
    sim_input.species[1].initial_value = (k2 / k1) * pow(initX, 2);

    ChemSolver chem_solver(sim_input);

    // run the problem on device
    real_type_2d_view results = chem_solver.get_results();

    // create mirror view and deep copy to host
    auto results_host = Kokkos::create_mirror_view(results);
    Kokkos::deep_copy(results_host, results);

    // eq (4) tendency should be positive, eq (5) negative
    const Real val00 = results_host(0, 0);
    const Real val01 = results_host(0, 1);
    REQUIRE(FloatingPoint<Real>::zero(val00, FloatingPoint<Real>::zero_tol));
    REQUIRE(FloatingPoint<Real>::zero(val01, FloatingPoint<Real>::zero_tol));

    // we ran for two batches, so be sure that the same inputs give the same
    // outputs for different batches
    const Real val10 = results_host(1, 0);
    const Real val11 = results_host(1, 1);
    REQUIRE(FloatingPoint<Real>::zero(val00 - val10,
                                      FloatingPoint<Real>::zero_tol));
    REQUIRE(FloatingPoint<Real>::zero(val01 - val11,
                                      FloatingPoint<Real>::zero_tol));

  }

  SECTION("dark side of terminator"){

    // change the k1 reaction rate to zero, corresponding to column on dark side
    sim_input.reactions[0].rate_coefficients["A"] = 0.0;
    {

      ChemSolver chem_solver(sim_input);

      // run the problem on device
      real_type_2d_view results = chem_solver.get_results();

      // create mirror view and deep copy to host
      auto results_host = Kokkos::create_mirror_view(results);
      Kokkos::deep_copy(results_host, results);

      // eq (4) tendency should be positive, eq (5) negative
      const Real val00 = results_host(0, 0);
      const Real val01 = results_host(0, 1);
      REQUIRE(FloatingPoint<Real>::in_bounds(val00, -1.0e14, 0.0,
                                             FloatingPoint<Real>::zero_tol));
      REQUIRE(FloatingPoint<Real>::in_bounds(val01, 0.0, 1.0e14,
                                             FloatingPoint<Real>::zero_tol));
      // we ran for two batches, so be sure that the same inputs give the same
      // outputs for different batches
      const Real val10 = results_host(1, 0);
      const Real val11 = results_host(1, 1);
      REQUIRE(FloatingPoint<Real>::zero(val00 - val10,
                                        FloatingPoint<Real>::zero_tol));
      REQUIRE(FloatingPoint<Real>::zero(val01 - val11,
                                        FloatingPoint<Real>::zero_tol));
    }
  }

}

