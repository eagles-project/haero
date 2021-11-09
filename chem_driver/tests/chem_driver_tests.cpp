#include "catch2/catch.hpp"
#include "chem_driver/chem_driver.hpp"
#include "util/read_solver_params.hpp"
#include "haero/floating_point.hpp"
#include "haero/haero.hpp"

using namespace haero;
using namespace haero::chem_driver;

TEST_CASE("TChem tendency computation tests", "haero_unit_tests") {

  SECTION("simple arrhenius") {
    std::string prefix = HAERO_TEST_DATA_DIR;
    std::string input_file = prefix;
    input_file += "/config_arrhenius.yaml";
    std::string refsol_file = prefix;
    refsol_file += "/camp_arrhenius.yaml";

    Real tbeg = 0.0;
    Real tend = 100.0;
    Real dt = 1.0;
    Real cur_time = tbeg;
    ordinal_type nsteps = 101;
    ordinal_type nspec = 4;

    auto refsol = Real_2d_view_host("refsol", nspec, nsteps);
    std::vector<std::string> specs = {"A", "B", "C", "D"};

    auto root = YAML::LoadFile(refsol_file);
    for (int i = 0; i < nspec; ++i) {
      auto node = root[specs[i]];
      int j = 0;
      for (auto iter : node) {
        refsol(i, j) = iter.as<Real>();
        ++j;
      }
    }

    auto solution = Real_2d_view_host("solution", nspec, nsteps);
    auto hsolution = Kokkos::create_mirror_view(solution);

    ChemSolver chem_solver(input_file);

    auto cur_sol = chem_solver.get_state();

    Kokkos::parallel_for(
        "copy_sol", nspec,
        KOKKOS_LAMBDA(const int& i) { solution(i, 0) = cur_sol(0, i + 3); });

    for (int i = 1; i < nsteps; ++i) {
      chem_solver.time_integrate(cur_time, cur_time + dt);
      cur_sol = chem_solver.get_state();
      cur_time += dt;
      Kokkos::parallel_for(
          "copy_sol", nspec,
          KOKKOS_LAMBDA(const int& j) { solution(j, i) = cur_sol(0, j + 3); });
    }
    Kokkos::deep_copy(hsolution, solution);

    Real tol = 1e-5;
    auto error = Real_1d_view_host("error", nspec);
    for (int i = 0; i < nspec; ++i) {
      error(i) = 0.0;
      for (int j = 0; j < nsteps; ++j) {
        error(i) += pow(hsolution(i, j) - refsol(i, j), 2);
      }
      error(i) = sqrt(error(i));
      REQUIRE(FloatingPoint<Real>::zero(error(i), tol));
    }

    tol = 1e-8;
    for (int i = 0; i < nspec; ++i)
    {
      REQUIRE(FloatingPoint<Real>::zero(abs(hsolution(i, nsteps - 1) - refsol(i, nsteps - 1)), tol));
    }
  }

  SECTION("simple troe") {
    std::string prefix = HAERO_TEST_DATA_DIR;
    std::string input_file = prefix;
    input_file += "/config_troe.yaml";
    std::string refsol_file = prefix;
    refsol_file += "/camp_troe.yaml";

    Real tbeg = 0.0;
    Real tend = 100.0;
    Real dt = 1.0;
    Real cur_time = tbeg;
    ordinal_type nsteps = 101;
    ordinal_type nspec = 3;

    auto refsol = Real_2d_view_host("refsol", nspec, nsteps);
    std::vector<std::string> specs = {"A", "B", "C"};

    auto root = YAML::LoadFile(refsol_file);
    for (int i = 0; i < nspec; ++i) {
      auto node = root[specs[i]];
      int j = 0;
      for (auto iter : node) {
        refsol(i, j) = iter.as<Real>();
        ++j;
      }
    }

    auto solution = Real_2d_view_host("solution", nspec, nsteps);
    auto hsolution = Kokkos::create_mirror_view(solution);

    ChemSolver chem_solver(input_file);

    auto cur_sol = chem_solver.get_state();

    Kokkos::parallel_for(
        "copy_sol", nspec,
        KOKKOS_LAMBDA(const int& i) { solution(i, 0) = cur_sol(0, i + 3); });

    for (int i = 1; i < nsteps; ++i) {
      chem_solver.time_integrate(cur_time, cur_time + dt);
      cur_sol = chem_solver.get_state();
      cur_time += dt;
      Kokkos::parallel_for(
          "copy_sol", nspec,
          KOKKOS_LAMBDA(const int& j) { solution(j, i) = cur_sol(0, j + 3); });
    }
    Kokkos::deep_copy(hsolution, solution);

    Real tol = 1e-5;
    auto error = Real_1d_view_host("error", nspec);
    for (int i = 0; i < nspec; ++i) {
      error(i) = 0.0;
      for (int j = 0; j < nsteps; ++j) {
        error(i) += pow(hsolution(i, j) - refsol(i, j), 2);
      }
      error(i) = sqrt(error(i));
      REQUIRE(FloatingPoint<Real>::zero(error(i), tol));
    }

    tol = 1e-7;
    for (int i = 0; i < nspec; ++i)
    {
      REQUIRE(FloatingPoint<Real>::zero(abs(hsolution(i, nsteps - 1) - refsol(i, nsteps - 1)), tol));
    }
  }

  SECTION("simple photolysis") {
    std::string prefix = HAERO_TEST_DATA_DIR;
    std::string input_file = prefix;
    input_file += "/config_photolysis.yaml";
    std::string refsol_file = prefix;
    refsol_file += "/camp_photolysis.yaml";

    Real tbeg = 0.0;
    Real tend = 100.0;
    Real dt = 1.0;
    Real cur_time = tbeg;
    ordinal_type nsteps = 101;
    ordinal_type nspec = 3;

    auto refsol = Real_2d_view_host("refsol", nspec, nsteps);
    std::vector<std::string> specs = {"A", "B", "C"};

    auto root = YAML::LoadFile(refsol_file);
    for (int i = 0; i < nspec; ++i) {
      auto node = root[specs[i]];
      int j = 0;
      for (auto iter : node) {
        refsol(i, j) = iter.as<Real>();
        ++j;
      }
    }

    auto solution = Real_2d_view_host("solution", nspec, nsteps);
    auto hsolution = Kokkos::create_mirror_view(solution);

    ChemSolver chem_solver(input_file);

    auto cur_sol = chem_solver.get_state();

    Kokkos::parallel_for(
        "copy_sol", nspec,
        KOKKOS_LAMBDA(const int& i) { solution(i, 0) = cur_sol(0, i + 3); });

    for (int i = 1; i < nsteps; ++i) {
      chem_solver.time_integrate(cur_time, cur_time + dt);
      cur_sol = chem_solver.get_state();
      cur_time += dt;
      Kokkos::parallel_for(
          "copy_sol", nspec,
          KOKKOS_LAMBDA(const int& j) { solution(j, i) = cur_sol(0, j + 3); });
    }
    Kokkos::deep_copy(hsolution, solution);

    Real tol = 1e-5;
    auto error = Real_1d_view_host("error", nspec);
    for (int i = 0; i < nspec; ++i) {
      error(i) = 0.0;
      for (int j = 0; j < nsteps; ++j) {
        error(i) += pow(hsolution(i, j) - refsol(i, j), 2);
      }
      error(i) = sqrt(error(i));
      REQUIRE(FloatingPoint<Real>::zero(error(i), tol));
    }

    tol = 1e-9;
    for (int i = 0; i < nspec; ++i)
    {
      REQUIRE(FloatingPoint<Real>::zero(abs(hsolution(i, nsteps - 1) - refsol(i, nsteps - 1)), tol));
    }
  }
}
