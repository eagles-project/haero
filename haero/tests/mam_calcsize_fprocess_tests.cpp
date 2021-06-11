#include "catch2/catch.hpp"
#include "haero/model.hpp"
#include "mam_calcsize_test_bridge.hpp"

using namespace haero;

TEST_CASE("compute_diameter", "mam_calcsize_fprocess") {

  static constexpr Real vol2num = 2;
  //compute diameter
  //MAMCalcsizeProcess::compute_diameter(vol2num);
  Real diameter = compute_diameter_bridge(vol2num);

  REQUIRE(diameter==3);
}
