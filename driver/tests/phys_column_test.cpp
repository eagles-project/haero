#include "driver/phys_column.hpp"
#include "catch2/catch.hpp"

using namespace haero;
using namespace haero::driver;

TEST_CASE("phys_column", "") {
  const int nlev = 10;
  phys_column mycol(nlev);
}
