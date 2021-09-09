#include <iomanip>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/aerosol_species.hpp"
#include "haero/check.hpp"
#include "haero/diagnostics/kohler_solve.hpp"
#include "haero/diagnostics/mode_dry_particle_volume.hpp"
#include "haero/diagnostics/mode_hygroscopicity.hpp"
#include "haero/diagnostics/mode_wet_radius.hpp"
#include "haero/math.hpp"
#include "haero/modal_aerosol_config.hpp"
#include "haero/mode.hpp"
// #include "ekat/mpi/ekat_comm.hpp"
// #include "ekat/logging/ekat_logger.hpp"

using namespace haero;

TEST_CASE("aerosol_activation_abdul-razaak_ghan_2000", "") {
  const auto config = create_ag00_2mode_config();
  std::cout << config.info_string();
}
