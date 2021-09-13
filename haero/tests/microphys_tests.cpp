#include <vector>
#include <iostream>

#include "haero/microphysics_functions.hpp"
#include "catch2/catch.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/logging/ekat_logger.hpp"
#include "haero/constants.hpp"
#include "haero/floating_point.hpp"

using namespace haero;
using namespace ekat::logger;

TEST_CASE("microphysics functions", "") {

  ekat::Comm comm;

  Logger<> logger("microphysics_test_log", Log::level::debug, comm);

  logger.info("starting microphysics test.");

  SECTION("surface tension") {

    const int npts = 86;
    const Real Tmax = 40;
    const Real Tmin = -45;
    const Real dT = (Tmax - Tmin)/(npts-1);

    std::vector<Real> temps_c(npts);
    std::vector<Real> surf_ten(npts);

    for (int i=0; i<npts; ++i) {
      const Real Ti = Tmin + i*dT;
      temps_c[i] = Ti;
      // convert Celsius to Kelvin
      const Real T = Ti + Constants::freezing_pt_h2o;
      const Real sigma_wa = surface_tension_water_air(T);
      surf_ten[i] = sigma_wa * 1000; // convert to c.g.s. from SI
    }

    logger.info("surface tension at 273.15K: {}",
      surface_tension_water_air(Constants::freezing_pt_h2o));

    REQUIRE(FloatingPoint<Real>::equiv(
      surface_tension_water_air(Constants::freezing_pt_h2o),
      Constants::surface_tension_h2o_air_273k, 3e-4));

    std::cout << "temps_c = [";
    for (int i=0; i<npts; ++i) {
      std::cout << temps_c[i] << (i < npts-1 ? ", " : "]\n");
    }
    std::cout << "surf_ten = [";
    for (int i=0; i<npts; ++i) {
      std::cout << surf_ten[i] << (i < npts-1 ? ", " : "]\n");
    }

  }
}
