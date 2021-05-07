#include "available_diagnostics.hpp"
#include "haero/modal_aerosol_config.hpp"
#include "haero/math_helpers.hpp"
#include "haero/mode.hpp"
#include "haero/aerosol_species.hpp"
#include "catch2/catch.hpp"
#include <iostream>

using namespace haero;

TEST_CASE ("wet_radius_diagnostic", "") {

  // setup basic level info
  const int nlev = 8;
  const int npacks = PackInfo::num_packs(nlev);

  // setup aerosol config
  const auto config = create_simple_test_config();
  std::cout << config.info_string();

}
