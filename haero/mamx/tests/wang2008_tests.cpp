#include <cstdio>

#include "catch2/catch.hpp"
#include "haero/processes/wang2008.hpp"

using namespace haero;

// These tests exercise the planetary boundary layer nucleation rates published
// by Wang et al (2008).

TEST_CASE("wang2008_first_order") {
  PackType c_h2so4(1e7);
  REQUIRE(FloatingPoint<PackType>::equiv(
      1e-6 * c_h2so4, wang2008::first_order_pbl_nucleation_rate(c_h2so4)));
}

TEST_CASE("wang2008_second_order") {
  PackType c_h2so4(1e14);
  REQUIRE(FloatingPoint<PackType>::equiv(
      1e-12 * c_h2so4, wang2008::second_order_pbl_nucleation_rate(c_h2so4)));
}
