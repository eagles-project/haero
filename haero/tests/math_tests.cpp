// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "math_tests.hpp"

using namespace haero;

TEST_CASE("haero_math_powers", "") {
  const Real arg = 5;

  REQUIRE(square(arg) == 25);
  REQUIRE(FloatingPoint<Real>::equiv(square(arg), 25.0));

  REQUIRE(cube(arg) == 125);
  REQUIRE(FloatingPoint<Real>::equiv(cube(arg), 125.0));
}

TEST_CASE("haero_math_rootfinding", "") {
  const Real cubic_root = sqrt(3.0 / 5.0);
  const Real quartic_root = sqrt((15.0 + 2 * sqrt(30.0)) / 35.0);

  using cubic_leg_poly = math::LegendreCubic<Real>;
  using quartic_leg_poly = math::LegendreQuartic<Real>;

  cubic_leg_poly p3;
  quartic_leg_poly p4;

  const Real conv_tol = 100 * FloatingPoint<Real>::zero_tol;
  std::cout << "convergence tolerance = " << conv_tol << "\n";

  const Real a0 = 0.5;
  const Real b0 = 1.0;
  const Real x0 = 1.0;

  SECTION("Newton solve") {

    auto cubic_solver =
        math::NewtonSolver<cubic_leg_poly>(x0, a0, b0, conv_tol, p3);
    const Real cubic_sol = cubic_solver.solve();

    std::cout << "newton cubic_sol rel. error = "
              << std::abs(cubic_sol - cubic_root) / cubic_root
              << " n_iter = " << cubic_solver.counter << "\n";

    REQUIRE(cubic_sol == Approx(cubic_root));

    auto quartic_solver =
        math::NewtonSolver<quartic_leg_poly>(x0, a0, b0, conv_tol, p4);
    const Real quartic_sol = quartic_solver.solve();

    std::cout << "newton quartic_sol rel. error = "
              << std::abs(quartic_sol - quartic_root) / quartic_root
              << " n_iter = " << quartic_solver.counter << "\n";

    REQUIRE(quartic_sol == Approx(quartic_root));
  }

  SECTION("bisection_solve") {
    auto cubic_solver =
        math::BisectionSolver<cubic_leg_poly>(x0, a0, b0, conv_tol, p3);
    const Real cubic_sol = cubic_solver.solve();
    std::cout << "bisection cubic_sol rel. error = "
              << std::abs(cubic_sol - cubic_root) / cubic_root
              << " n_iter = " << cubic_solver.counter << "\n";

    REQUIRE(cubic_sol == Approx(cubic_root));

    auto quartic_solver =
        math::BisectionSolver<quartic_leg_poly>(x0, a0, b0, conv_tol, p4);
    const Real quartic_sol = quartic_solver.solve();
    std::cout << "bisection quartic_sol rel. error = "
              << std::abs(quartic_sol - quartic_root) / quartic_root
              << " n_iter = " << quartic_solver.counter << "\n";
    REQUIRE(quartic_sol == Approx(quartic_root));
  }

  SECTION("bracketed_newton_solve") {
    auto cubic_solver =
        math::BracketedNewtonSolver<cubic_leg_poly>(x0, a0, b0, conv_tol, p3);
    const Real cubic_sol = cubic_solver.solve();
    std::cout << "bracketed newton cubic_sol rel. error = "
              << std::abs(cubic_sol - cubic_root) / cubic_root
              << " n_iter = " << cubic_solver.counter << "\n";

    REQUIRE(cubic_sol == Approx(cubic_root));

    auto quartic_solver =
        math::BracketedNewtonSolver<quartic_leg_poly>(x0, a0, b0, conv_tol, p4);
    const Real quartic_sol = quartic_solver.solve();
    std::cout << "bracketed newton quartic_sol rel. error = "
              << std::abs(quartic_sol - quartic_root) / quartic_root
              << " n_iter = " << quartic_solver.counter << "\n";
    REQUIRE(quartic_sol == Approx(quartic_root));
  }
}

TEST_CASE("no_root", "") {
  using quadratic_poly = math::MonicParabola<Real>;

  quadratic_poly qpoly;

  const Real conv_tol = 100 * FloatingPoint<Real>::zero_tol;
  std::cout << "convergence tolerance = " << conv_tol << "\n";

  const Real x0 = 0;
  const Real a0 = -1;
  const Real b0 = 1;

  auto quadratic_solver =
      math::NewtonSolver<quadratic_poly>(x0, a0, b0, conv_tol, qpoly);
  const Real qsol = quadratic_solver.solve();
  std::cout << "newton quadratic_sol = " << qsol
            << " n_iter = " << quadratic_solver.counter << "\n";

  REQUIRE(quadratic_solver.fail);
  REQUIRE(isnan(qsol));
}
