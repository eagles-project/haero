#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/haero.hpp"
#include "haero/math.hpp"
#include "haero/floating_point.hpp"
#include "ekat/ekat_pack_math.hpp"

using namespace haero;

TEST_CASE("haero_math_powers", "") {
  SECTION("no packs") {
  const Real arg = 5;

  REQUIRE(square(arg) == 25);
  REQUIRE(FloatingPoint<Real>::equiv(square(arg), 25.0));

  REQUIRE(cube(arg) == 125);
  REQUIRE(FloatingPoint<Real>::equiv(cube(arg), 125.0));

  }

#ifndef __CUDACC__
  SECTION("packs") {
    typedef ekat::Pack<Real,4> real_pack;
    const auto arg = real_pack(5);

    REQUIRE(FloatingPoint<real_pack>::equiv(square(arg), 25.0));

    REQUIRE(FloatingPoint<real_pack>::equiv(cube(arg), 125.0));
  }
#endif
}

TEST_CASE("haero_math_rootfinding_no_packs", "") {
  const Real cubic_root = sqrt(3.0/5.0);
  const Real quartic_root = sqrt((15.0 + 2*sqrt(30.0))/35.0);

  using cubic_leg_poly = math::LegendreCubic<Real>;
  using quartic_leg_poly = math::LegendreQuartic<Real>;

  cubic_leg_poly p3;
  quartic_leg_poly p4;

  const Real conv_tol = 100*FloatingPoint<Real>::zero_tol;

  SECTION("Newton solve") {
    const Real x0 = 1.0;

    auto cubic_solver = math::ScalarNewtonSolver<cubic_leg_poly>(x0, conv_tol, p3);
    const Real cubic_sol = cubic_solver.solve();

    std::cout << "newton cubic_sol rel. error = "
              << std::abs(cubic_sol - cubic_root) / cubic_root
              << " n_iter = " << cubic_solver.counter << "\n";

    REQUIRE(cubic_sol == Approx(cubic_root));

    auto quartic_solver = math::ScalarNewtonSolver<quartic_leg_poly>(x0, conv_tol, p4);
    const Real quartic_sol = quartic_solver.solve();

    std::cout << "newton quartic_sol rel. error = "
              << std::abs(quartic_sol - quartic_root) / quartic_root
              << " n_iter = " << quartic_solver.counter << "\n";

    REQUIRE(quartic_sol == Approx(quartic_root));
  }

  SECTION("bisection_solve") {
    const Real a0 = 0.5;
    const Real b0 = 1.0;
    auto cubic_solver = math::BisectionSolver<cubic_leg_poly>(a0, b0, conv_tol, p3);
    const Real cubic_sol = cubic_solver.solve();
    std::cout << "bisection cubic_sol rel. error = "
              << std::abs(cubic_sol - cubic_root) / cubic_root
              << " n_iter = " << cubic_solver.counter << "\n";

    REQUIRE(cubic_sol == Approx(cubic_root));

    auto quartic_solver = math::BisectionSolver<quartic_leg_poly>(a0, b0, conv_tol, p4);
    const Real quartic_sol = quartic_solver.solve();
    std::cout << "bisection quartic_sol rel. error = "
              << std::abs(quartic_sol - quartic_root) / quartic_root
              << " n_iter = " << quartic_solver.counter << "\n";
    REQUIRE(quartic_sol == Approx(quartic_root));
  }

  SECTION("bracketed_newton_solve") {
    const Real a0 = 0.5;
    const Real b0 = 1.0;
    const Real x0 = 1.0;
    auto cubic_solver = math::BracketedNewtonSolver<cubic_leg_poly>(x0, a0, b0, conv_tol, p3);
    const Real cubic_sol = cubic_solver.solve();
    std::cout << "bracketed newton cubic_sol rel. error = "
              << std::abs(cubic_sol - cubic_root) / cubic_root
              << " n_iter = " << cubic_solver.counter << "\n";

    REQUIRE(cubic_sol == Approx(cubic_root));

    auto quartic_solver = math::BracketedNewtonSolver<quartic_leg_poly>(x0, a0, b0, conv_tol, p4);
    const Real quartic_sol = quartic_solver.solve();
    std::cout << "bracketed newton quartic_sol rel. error = "
              << std::abs(quartic_sol - quartic_root) / quartic_root
              << " n_iter = " << quartic_solver.counter << "\n";
    REQUIRE(quartic_sol == Approx(quartic_root));
  }
}
#ifndef __CUDACC__
TEST_CASE("haer_math_rootfinding_packs","") {
  typedef ekat::Pack<Real,4> pack_type;
  const Real cubic_root = sqrt(3.0/5.0);
  const Real quartic_root = sqrt((15.0 + 2*sqrt(30.0))/35.0);

  using cubic_leg_poly = math::LegendreCubic<pack_type>;
  using quartic_leg_poly = math::LegendreQuartic<pack_type>;

  cubic_leg_poly p3;
  quartic_leg_poly p4;

  const Real conv_tol = 100*FloatingPoint<Real>::zero_tol;

  SECTION("newton solve") {
    const pack_type x0(1.0);

    auto cubic_solver = math::ScalarNewtonSolver<cubic_leg_poly>(x0, conv_tol, p3);
    const pack_type cubic_sol = cubic_solver.solve();

    std::cout << "newton cubic_sol rel. error = "
              << abs(cubic_sol - cubic_root) / cubic_root
              << " n_iter = " << cubic_solver.counter << "\n";

    REQUIRE(FloatingPoint<pack_type>::rel(cubic_sol, cubic_root));

    auto quartic_solver = math::ScalarNewtonSolver<quartic_leg_poly>(x0, conv_tol, p4);
    const pack_type quartic_sol = quartic_solver.solve();

    std::cout << "newton quartic_sol rel. error = "
              << abs(quartic_sol - quartic_root) / quartic_root
              << " n_iter = " << quartic_solver.counter << "\n";

    REQUIRE(FloatingPoint<pack_type>::rel(quartic_sol, quartic_root));
  }
  SECTION("bisection solve") {
    const pack_type a0(0.5);
    const pack_type b0(1.0);

    auto cubic_solver = math::BisectionSolver<cubic_leg_poly>(a0,b0, conv_tol, p3);
    const pack_type cubic_sol = cubic_solver.solve();

    std::cout << "bisection cubic_sol rel. error = "
              << abs(cubic_sol - cubic_root) / cubic_root
              << " n_iter = " << cubic_solver.counter << "\n";

    REQUIRE(FloatingPoint<pack_type>::rel(cubic_sol, cubic_root, conv_tol));

    auto quartic_solver = math::BisectionSolver<quartic_leg_poly>(a0, b0, conv_tol, p4);
    const pack_type quartic_sol = quartic_solver.solve();

    std::cout << "bisection quartic_sol rel. error = "
              << abs(quartic_sol - quartic_root) / quartic_root
              << " n_iter = " << quartic_solver.counter << "\n";

    REQUIRE(FloatingPoint<pack_type>::rel(quartic_sol, quartic_root, conv_tol));
  }
    SECTION("bracketed newton solve") {
    const pack_type x0(1.0);
    const Real a0(0.5);
    const Real b0(1.0);

    auto cubic_solver = math::BracketedNewtonSolver<cubic_leg_poly>(x0, a0,b0, conv_tol, p3);
    const pack_type cubic_sol = cubic_solver.solve();

    std::cout << "bracketed newton cubic_sol rel. error = "
              << abs(cubic_sol - cubic_root) / cubic_root
              << " n_iter = " << cubic_solver.counter << "\n";

    REQUIRE(FloatingPoint<pack_type>::rel(cubic_sol, cubic_root));

    auto quartic_solver = math::BracketedNewtonSolver<quartic_leg_poly>(x0,a0, b0, conv_tol, p4);
    const pack_type quartic_sol = quartic_solver.solve();

    std::cout << "bracketed newton quartic_sol rel. error = "
              << abs(quartic_sol - quartic_root) / quartic_root
              << " n_iter = " << quartic_solver.counter << "\n";

    REQUIRE(FloatingPoint<pack_type>::rel(quartic_sol, quartic_root));
  }
}
#endif




