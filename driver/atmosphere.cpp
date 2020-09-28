#include "atmosphere.hpp"
#include "ekat/ekat_assert.hpp"
#include "haero/floating_point.hpp"

namespace haero {

using fp_helper = FloatingPoint<Real>;

AtmosphericConditions hydrostatic_conditions(const Real& p0, const Real& T0, const Real& Gamma,
  const Real& qv0, const Real& qv1) {

  /// check valid input
  EKAT_REQUIRE_MSG(fp_helper::in_bounds(p0, 100000, p_std_pa),
    "expected p0 appx. equal 100000 Pa; check units = Pa");
  EKAT_REQUIRE_MSG(fp_helper::in_bounds(T0, 273, 323),
    "unexpected T0, check units = K");
  EKAT_REQUIRE_MSG(fp_helper::in_bounds(Gamma, 0, 0.02),
    "unexpected lapse rate, check units = K/m");
  EKAT_REQUIRE_MSG(fp_helper::in_bounds(qv0, 0, 0.1),
    "unexpected water vapor mixing ratio; check units = kg/kg");
  EKAT_REQUIRE_MSG(qv1 >= 0, "nonnegative decay rate required.");

  AtmosphericConditions result;
  result.params.hydrostatic.p0 = p0;
  result.params.hydrostatic.T0 = T0;
  result.params.hydrostatic.lapse_rate = Gamma;
  result.params.hydrostatic.qv0 = qv0;
  result.params.hydrostatic.qv1 = qv1;
  return result;
}

} // namespace haero
