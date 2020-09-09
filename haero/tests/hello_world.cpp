#include "catch2/catch.hpp"

namespace {

void hello_world()
{
  REQUIRE(1 == 1);
  REQUIRE(2+2 != 5);
}

TEST_CASE("hello_world", "[haero_unit_tests]")
{
  hello_world();
}

} // anonymous namespace
