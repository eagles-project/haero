#include "haero/microphysics_functions.hpp"
#include "catch2/catch.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/logging/ekat_logger.hpp"

using namespace ekat::logger;

TEST_CASE("microphysics functions", "") {

  ekat::Comm comm;

  Logger<> logger("microphysics_test_log", Log::level::debug, comm);




}
