#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "haero/mode.hpp"
#include "haero/physical_constants.hpp"

using namespace haero;

TEST_CASE("mode_ctor", "") {
  // Create a representation of the Aitken mode.
  Mode aitken("Aitken", 8.7e-9, 5.2e-8, 1.6, 0.35, 0.8);
  REQUIRE(aitken.name() == "Aitken");
  REQUIRE(aitken.min_diameter == Real(8.7e-9));
  REQUIRE(aitken.max_diameter == Real(5.2e-8));
  REQUIRE(aitken.mean_std_dev == Real(1.6));

  // For single precision calculations, we expect accuracy till 7 decimal digits
  // (1e-7)
  static constexpr Real exp_sngl_prc_accuracy =
      1e-7;  // expected single precision accuracy

  // compute min_vol_to_num ratio
  const Real comp_min_vol_to_num_ratio =
      1 / (constants::pi_sixth * cube(5.2e-8) * exp(4.5 * square(log(1.6))));

  // compute relative difference for min_vol_to_num_ratio
  const Real rel_diff_min_vol =
      (aitken.min_vol_to_num_ratio<Real>() - comp_min_vol_to_num_ratio) /
      aitken.min_vol_to_num_ratio<Real>();
  REQUIRE(rel_diff_min_vol <= exp_sngl_prc_accuracy);

  // compute max_vol_to_num ratio
  const Real comp_max_vol_to_num_ratio =
      1 / (constants::pi_sixth * cube(8.7e-9) * exp(4.5 * square(log(1.6))));

  // compute relative difference for max_vol_to_num_ratio
  const Real rel_diff_max_vol =
      (aitken.max_vol_to_num_ratio<Real>() - comp_max_vol_to_num_ratio) /
      aitken.max_vol_to_num_ratio<Real>();
  REQUIRE(rel_diff_max_vol <= exp_sngl_prc_accuracy);
}
