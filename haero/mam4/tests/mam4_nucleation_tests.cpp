#include <haero/mam4/mam4.hpp>

#include <catch2/catch.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

TEST_CASE("test_constructor", "mam4_nucleation_process") {
  haero::mam4::AeroConfig mam4_config;
  haero::mam4::NucleationProcess process(mam4_config);
}

