#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/haero.hpp"
#include "haero/math_helpers.hpp"

using namespace haero;

TEST_CASE("powers", "") {
  const int argint = 5;
  const Real argreal = 5;

  REQUIRE(square(argint) == 25);
  REQUIRE(FloatingPoint<Real>::equiv(square(argreal), 25.0));

  REQUIRE(cube(argint) == 125);
  REQUIRE(FloatingPoint<Real>::equiv(cube(argreal), 125.0));
}

TEST_CASE("rootfinding-Real", "") {
  using namespace math;

  const Real a0 = 0.5;
  const Real b0 = 1.0;

  const Real cubic_root = sqrt(3.0 / 5.0);
  const Real quartic_root = sqrt((15.0 + 2.0 * sqrt(30.0)) / 35.0);

  const Real x0 = 1.0;

  const LegendreCubic<Real> p3;
  const LegendreQuartic<Real> p4;

  const Real conv_tol = 100 * FloatingPoint<Real>::zero_tol;

  SECTION("newton_solve") {
    auto cubic_solver =
        ScalarNewtonSolver<LegendreCubic<Real>>(x0, conv_tol, p3);
    const Real cubic_sol = cubic_solver.solve();
    std::cout << "newton cubic_sol rel. error = "
              << std::abs(cubic_sol - cubic_root) / cubic_root
              << " n_iter = " << cubic_solver.counter << "\n";

    auto quartic_solver =
        ScalarNewtonSolver<LegendreQuartic<Real>>(x0, conv_tol, p4);
    const Real quartic_sol = quartic_solver.solve();
    std::cout << "newton quartic_sol rel. error = "
              << std::abs(quartic_sol - quartic_root) / quartic_root
              << " n_iter = " << quartic_solver.counter << "\n";

    REQUIRE(FloatingPoint<Real>::equiv(cubic_sol, cubic_root, conv_tol));
    REQUIRE(FloatingPoint<Real>::equiv(quartic_sol, quartic_root, conv_tol));
  }

  SECTION("bisection_solve") {
    auto cubic_solver =
        BisectionSolver<LegendreCubic<Real>>(a0, b0, conv_tol, p3);
    const Real cubic_sol = cubic_solver.solve();
    std::cout << "bisection cubic_sol rel. error = "
              << std::abs(cubic_sol - cubic_root) / cubic_root
              << " n_iter = " << cubic_solver.counter << "\n";
    REQUIRE(FloatingPoint<Real>::equiv(cubic_sol, cubic_root, conv_tol));

    auto quartic_solver =
        BisectionSolver<LegendreQuartic<Real>>(a0, b0, conv_tol, p4);
    const Real quartic_sol = quartic_solver.solve();
    std::cout << "bisection quartic_sol rel. error = "
              << std::abs(quartic_sol - quartic_root) / quartic_root
              << " n_iter = " << quartic_solver.counter << "\n";
    REQUIRE(FloatingPoint<Real>::equiv(quartic_sol, quartic_root, conv_tol));
  }

  SECTION("bracketed newton solve") {
    auto cubic_solver =
        BracketedNewtonSolver<LegendreCubic<Real>>(x0, a0, b0, conv_tol, p3);
    const Real cubic_sol = cubic_solver.solve();
    std::cout << "bracketed newton cubic sol. rel. error = "
              << std::abs(cubic_sol - cubic_root) / cubic_root
              << " n_iter = " << cubic_solver.counter << "\n";

    auto quartic_solver =
        BracketedNewtonSolver<LegendreQuartic<Real>>(x0, a0, b0, conv_tol, p4);
    const Real quartic_sol = quartic_solver.solve();
    std::cout << "bracketed newton quartic sol rel. error = "
              << std::abs(quartic_sol - quartic_root) / quartic_root
              << " n_iter = " << quartic_solver.counter << "\n";

    REQUIRE(FloatingPoint<Real>::equiv(cubic_sol, cubic_root, conv_tol));
    REQUIRE(FloatingPoint<Real>::equiv(quartic_sol, quartic_root, conv_tol));
  }
}

TEST_CASE("rootfinding-PackType", "") {
  using namespace math;

  const auto cubic_root = PackType(sqrt(3.0 / 5.0));
  const auto quartic_root = PackType(sqrt((15.0 + 2.0 * sqrt(30.0)) / 35.0));

  const Real conv_tol = 100 * FloatingPoint<Real>::zero_tol;

  const Real a0 = 0.5;
  const Real b0 = 1.0;
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

  Kokkos::View<ekat::Pack<int, HAERO_PACK_SIZE>[2]> niterations(
        "niterations");
  auto h_niterations = Kokkos::create_mirror_view(niterations);
  h_niterations(0) = 0;
  h_niterations(1) = 0;
  Kokkos::deep_copy(niterations, h_niterations);

  SECTION("NewtonSolver") {
    Kokkos::parallel_for(
        1, KOKKOS_LAMBDA(const int i) {
          const LegendreCubic<PackType> p3;
          const LegendreQuartic<PackType> p4;
          auto cubic_solver = ScalarNewtonSolver<LegendreCubic<PackType>>(
              num_sol(0), conv_tol, p3);
          auto quartic_solver = ScalarNewtonSolver<LegendreQuartic<PackType>>(
              num_sol(1), conv_tol, p4);
          num_sol(0) = cubic_solver.solve();
          num_sol(1) = quartic_solver.solve();
          niterations(0) = cubic_solver.counter;
          niterations(1) = quartic_solver.counter;
        });

    Kokkos::deep_copy(h_num_sol, num_sol);
    Kokkos::deep_copy(h_niterations, niterations);

    std::cout << "packed newton rel. error (cubic) = "
              << abs(h_num_sol(0) - cubic_root) / cubic_root
              << ", niterations = " << h_niterations(0) << "\n";
    std::cout << "packed newton rel. error (quartic) = "
              << abs(h_num_sol(1) - quartic_root) / quartic_root
              << ", niterations = " << h_niterations(1) << "\n";

    REQUIRE(FloatingPoint<PackType>::equiv(h_num_sol(0), cubic_root, conv_tol));
    REQUIRE(
        FloatingPoint<PackType>::equiv(h_num_sol(1), quartic_root, conv_tol));
  }

  SECTION("Bisection Solver") {

    Kokkos::parallel_for(1, KOKKOS_LAMBDA (const int i) {
      const LegendreCubic<PackType> p3;
      const LegendreQuartic<PackType> p4;
      auto cubic_solver = BisectionSolver<LegendreCubic<PackType>>(PackType(a0), PackType(b0), conv_tol, p3);
      auto quartic_solver = BisectionSolver<LegendreQuartic<PackType>>(PackType(a0), PackType(b0), conv_tol, p4);
      num_sol(0) = cubic_solver.solve();
      num_sol(1) = quartic_solver.solve();
      niterations(0) = cubic_solver.counter;
      niterations(1) = quartic_solver.counter;
    });

    Kokkos::deep_copy(h_num_sol, num_sol);
    Kokkos::deep_copy(h_niterations, niterations);

    std::cout << "packed bisection rel. error (cubic) = "
              << abs(h_num_sol(0) - cubic_root) / cubic_root
              << ", niterations = " << h_niterations(0) << "\n";
    std::cout << "packed bisection rel. error (quartic) = "
              << abs(h_num_sol(1) - quartic_root) / quartic_root
              << ", niterations = " << h_niterations(1) << "\n";

    REQUIRE(FloatingPoint<PackType>::equiv(h_num_sol(0), cubic_root, conv_tol));
    REQUIRE(
        FloatingPoint<PackType>::equiv(h_num_sol(1), quartic_root, conv_tol));
  }

  SECTION("Bracketed Newton Solver") {

    Kokkos::parallel_for(
        1, KOKKOS_LAMBDA(const int i) {
          const LegendreCubic<PackType> p3;
          const LegendreQuartic<PackType> p4;
          auto cubic_solver =
              BracketedNewtonSolver<LegendreCubic<PackType>>(
                  num_sol(0), a0, b0, conv_tol, p3);
          auto quartic_solver =
              BracketedNewtonSolver<LegendreQuartic<PackType>>(
                  num_sol(1), a0, b0, conv_tol, p4);
          num_sol(0) = cubic_solver.solve();
          num_sol(1) = quartic_solver.solve();
          niterations(0) = cubic_solver.counter;
          niterations(1) = quartic_solver.counter;
        });

    Kokkos::deep_copy(h_num_sol, num_sol);
    Kokkos::deep_copy(h_niterations, niterations);

    std::cout << "packed bracketed newton rel. error (cubic) = "
              << abs(h_num_sol(0) - cubic_root) / cubic_root
              << ", niterations = " << h_niterations(0) << "\n";
    std::cout << "packed bracketed newton rel. error (quartic) = "
              << abs(h_num_sol(1) - quartic_root) / quartic_root
              << ", niterations = " << h_niterations(1) << "\n";

    REQUIRE(FloatingPoint<PackType>::equiv(h_num_sol(0), cubic_root, conv_tol));
    REQUIRE(
        FloatingPoint<PackType>::equiv(h_num_sol(1), quartic_root, conv_tol));
  }
}
