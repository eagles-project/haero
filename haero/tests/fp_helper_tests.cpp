#include "haero/floating_point_helper.hpp"
#include "catch2/catch.hpp"
#include <iostream>

using namespace haero;

TEST_CASE("FloatingPointHelper-sp", "single precision") {

  const float test_tol = 1.0E-4;

  const float lower = 1;
  const float upper = 2;

  using fp_helper = FloatingPointHelper<float>;
  typedef FloatingPointHelper<float> fp_helper;

  REQUIRE_FALSE(fp_helper::zero(test_tol));
  REQUIRE(fp_helper::zero(0.5*test_tol, test_tol));

  REQUIRE_FALSE(fp_helper::equiv(lower, 0.5*upper + test_tol));
  REQUIRE(fp_helper::equiv(lower, 0.5*upper));

  REQUIRE(fp_helper::in_bounds(1.5, lower, upper));
  REQUIRE_FALSE(fp_helper::in_bounds(1-0.5*test_tol, lower, upper));
  REQUIRE(fp_helper::in_bounds(1-0.5*test_tol, lower, upper, test_tol));

  REQUIRE(std::isinf(lower/0));
  REQUIRE(lower * fp_helper::safe_denominator(0) == 0);
}

TEST_CASE("FloatingPointHelper-dp","double precision") {
  using fp_helper = FloatingPointHelper<double>;

  const double test_tol = 1.0E-4;

  const double lower = 1;
  const double upper = 2;
  REQUIRE_FALSE(fp_helper::zero(test_tol));
  REQUIRE(fp_helper::zero(0.5*test_tol, test_tol));

  REQUIRE_FALSE(fp_helper::equiv(lower, 0.5*upper + test_tol));
  REQUIRE(fp_helper::equiv(lower, 0.5*upper));

  REQUIRE(fp_helper::in_bounds(1.5, lower, upper));
  REQUIRE_FALSE(fp_helper::in_bounds(1-0.5*test_tol, lower, upper));
  REQUIRE(fp_helper::in_bounds(1-0.5*test_tol, lower, upper, test_tol));

  REQUIRE(std::isinf(lower/0));
  REQUIRE(lower * fp_helper::safe_denominator(0) == 0);

}

