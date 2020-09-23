#include "driver/dyn_column.hpp"
#include <iostream>

#include "catch2/catch.hpp"

using namespace haero;

TEST_CASE("dynamics_column", "") {
  const int nlev = 10;
  DynColumn mycol(nlev);

  std::cout << mycol.info_string();
}
