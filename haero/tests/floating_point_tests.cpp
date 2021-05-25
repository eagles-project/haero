#include <iostream>
#include <limits>

#include "catch2/catch.hpp"
#include "haero/floating_point.hpp"

using namespace haero;

TEST_CASE("FloatingPointHelper-sp", "single precision") {
  const float test_tol = 1.0E-4;

  const float lower = 1;
  const float upper = 2;

  REQUIRE_FALSE(FloatingPoint<float>::zero(test_tol));
  REQUIRE(FloatingPoint<float>::zero(0.5 * test_tol, test_tol));

  REQUIRE_FALSE(FloatingPoint<float>::equiv(lower, 0.5 * upper + test_tol));
  REQUIRE(FloatingPoint<float>::equiv(lower, 0.5 * upper));

  REQUIRE(FloatingPoint<float>::in_bounds(1.5, lower, upper));
  REQUIRE_FALSE(
      FloatingPoint<float>::in_bounds(1 - 0.5 * test_tol, lower, upper));
  REQUIRE(FloatingPoint<float>::in_bounds(1 - 0.5 * test_tol, lower, upper,
                                          test_tol));

  REQUIRE(std::isinf(lower / 0));
  REQUIRE(lower * FloatingPoint<float>::safe_denominator(0) == 0);

  std::cout << "std::numeric_limits<float>::epsilon() = "
            << std::numeric_limits<float>::epsilon() << "\n";
  std::cout << "std::numeric_limits<double>::epsilon() = "
            << std::numeric_limits<double>::epsilon() << "\n";
}

TEST_CASE("FloatingPointHelper-dp", "double precision") {
  const double test_tol = 1.0E-4;

  const double lower = 1;
  const double upper = 2;
  REQUIRE_FALSE(FloatingPoint<double>::zero(test_tol));
  REQUIRE(FloatingPoint<double>::zero(0.5 * test_tol, test_tol));

  REQUIRE_FALSE(FloatingPoint<double>::equiv(lower, 0.5 * upper + test_tol));
  REQUIRE(FloatingPoint<double>::equiv(lower, 0.5 * upper));

  REQUIRE(FloatingPoint<double>::in_bounds(1.5, lower, upper));
  REQUIRE_FALSE(
      FloatingPoint<double>::in_bounds(1 - 0.5 * test_tol, lower, upper));
  REQUIRE(FloatingPoint<double>::in_bounds(1 - 0.5 * test_tol, lower, upper,
                                           test_tol));

  REQUIRE(std::isinf(lower / 0));
  REQUIRE(lower * FloatingPoint<double>::safe_denominator(0) == 0);
}

TEST_CASE("FloatingPointHelper-packed", "packed") {
  const Real test_tol = 1.0e-4;
  const Real lower = 1;
  const Real upper = 2;

  REQUIRE_FALSE(FloatingPoint<PackType>::zero(PackType(test_tol)));
  REQUIRE(FloatingPoint<PackType>::zero(PackType(0.5 * test_tol), test_tol));

  REQUIRE_FALSE(FloatingPoint<PackType>::equiv(
      PackType(lower), 0.5 * PackType(upper) + test_tol));
  REQUIRE(
      FloatingPoint<PackType>::equiv(PackType(lower), 0.5 * PackType(upper)));

  REQUIRE(FloatingPoint<PackType>::in_bounds(PackType(1.5), lower, upper));
  REQUIRE_FALSE(FloatingPoint<PackType>::in_bounds(PackType(1) - 0.5 * test_tol,
                                                   lower, upper));
  REQUIRE(FloatingPoint<PackType>::in_bounds(PackType(1 - 0.5 * test_tol),
                                             lower, upper, test_tol));

  REQUIRE((PackType(lower) *
               FloatingPoint<PackType>::safe_denominator(PackType(0)) ==
           0)
              .all());
}
