#include "haero/haero.hpp"
#include "haero/math_helpers.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("powers", "") {
  const int argint = 5;
  const Real argreal = 5;

  REQUIRE(square(argint) == 25);
  REQUIRE(FloatingPoint<Real>::equiv(square(argreal), 25.0));

  REQUIRE(cube(argint) == 125);
  REQUIRE(FloatingPoint<Real>::equiv(cube(argreal), 125.0));
}

TEST_CASE("rootfinding-Real","") {

  using namespace math;

  const Real cubic_root = sqrt(3.0/5.0);
  const Real quartic_root = sqrt((15.0+2.0*sqrt(30.0))/35.0);

  const Real x0 = 1.0;

  const LegendreCubic<Real> p3;
  const LegendreQuartic<Real> p4;

  const Real conv_tol = 100*FloatingPoint<Real>::zero_tol;

  SECTION("newton_solve") {
    auto cubic_solver = ScalarNewtonSolver<LegendreCubic<Real>>(x0, conv_tol, p3);
    const Real cubic_sol = cubic_solver.solve();
    std::cout << "newton cubic_sol rel. error = " << std::abs(cubic_sol - cubic_root)/cubic_root
       << " n_iter = " << cubic_solver.counter << "\n";

    auto quartic_solver = ScalarNewtonSolver<LegendreQuartic<Real>>(x0, conv_tol, p4);
    const Real quartic_sol = quartic_solver.solve();
    std::cout << "newton quartic_sol rel. error = " << std::abs(quartic_sol - quartic_root)/quartic_root
      << " n_iter = " << quartic_solver.counter << "\n";

    REQUIRE(FloatingPoint<Real>::equiv(cubic_sol, cubic_root, conv_tol));
    REQUIRE(FloatingPoint<Real>::equiv(quartic_sol, quartic_root, conv_tol));
  }

  SECTION("bisection_solve") {

    const Real a0 = 0.5;
    const Real b0 = 1.0;

    auto cubic_solver = BisectionSolver<LegendreCubic<Real>>(a0, b0, conv_tol, p3);
    const Real cubic_sol = cubic_solver.solve();
    std::cout << "bisection cubic_sol rel. error = " << std::abs(cubic_sol - cubic_root)/cubic_root
      << " n_iter = " << cubic_solver.counter << "\n";
    REQUIRE(FloatingPoint<Real>::equiv(cubic_sol, cubic_root, conv_tol));

    auto quartic_solver = BisectionSolver<LegendreQuartic<Real>>(a0, b0, conv_tol, p4);
    const Real quartic_sol = quartic_solver.solve();
    std::cout << "bisection quartic_sol rel. error = " << std::abs(quartic_sol - quartic_root)/quartic_root
      << " n_iter = " << quartic_solver.counter << "\n";
    REQUIRE(FloatingPoint<Real>::equiv(quartic_sol, quartic_root, conv_tol));
  }
}

TEST_CASE("rootfinding-PackType", "") {
  using namespace math;

  const auto cubic_root = PackType(sqrt(3.0/5.0));
  const auto quartic_root = PackType(sqrt((15.0+2.0*sqrt(30.0))/35.0));

  const Real conv_tol = 100*FloatingPoint<Real>::zero_tol;

  const Real x0 = 1.0;

  Kokkos::View<PackType[2]> exact_roots("exact_roots");
  auto h_exact_roots = Kokkos::create_mirror_view(exact_roots);
  h_exact_roots(0) = cubic_root;
  h_exact_roots(1) = quartic_root;
  Kokkos::deep_copy(exact_roots, h_exact_roots);

  Kokkos::View<PackType[2]> num_sol("num_sol");
  auto h_num_sol = Kokkos::create_mirror_view(num_sol);
  h_num_sol(0) = PackType(x0);
  h_num_sol(1) = PackType(x0);
  Kokkos::deep_copy(num_sol, h_num_sol);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA (const int i) {
      const LegendreCubic<PackType> p3;
      const LegendreQuartic<PackType> p4;
      auto cubic_solver = ScalarNewtonSolver<LegendreCubic<PackType>>(num_sol(0), conv_tol, p3);
      auto quartic_solver = ScalarNewtonSolver<LegendreQuartic<PackType>>(num_sol(1), conv_tol, p4);
      num_sol(0) = cubic_solver.solve();
      num_sol(1) = quartic_solver.solve();
  });

  Kokkos::deep_copy(h_num_sol, num_sol);
  REQUIRE(FloatingPoint<PackType>::equiv(h_num_sol(0), cubic_root, conv_tol));
  REQUIRE(FloatingPoint<PackType>::equiv(h_num_sol(1), quartic_root, conv_tol));
}

