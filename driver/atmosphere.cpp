#include "atmosphere.hpp"
#include "ekat/ekat_assert.hpp"

namespace haero {

AtmosphericConditions stable_hydrostatic(const Real& refp, const Real& reft, const Real& gamma)
{
  EKAT_REQUIRE_MSG(gamma < dry_adiabatic_lapse_rate_K_per_m,
  "stable conditions require gamma < dry_adiabatic_lapse_rate");

  AtmosphericConditions result;
  result.params.hydrostatic.p0 = refp;
  result.params.hydrostatic.T0 = reft;
  result.params.hydrostatic.lapse_rate = gamma;
  return result;
}

AtmosphericConditions neutrally_stable_hydrostatic(const Real& refp, const Real& reft)
{
  AtmosphericConditions result;
  result.params.hydrostatic.p0 = refp;
  result.params.hydrostatic.T0 = reft;
  result.params.hydrostatic.lapse_rate = dry_adiabatic_lapse_rate_K_per_m;
  return result;
}

AtmosphericConditions unstable_hydrostatic(const Real& refp, const Real& reft, const Real& gamma)
{
  EKAT_REQUIRE_MSG(gamma > dry_adiabatic_lapse_rate_K_per_m,
    "unstable conditions require gamma > dry_adiabatic_lapse_rate");

  AtmosphericConditions result;
  result.params.hydrostatic.p0 = refp;
  result.params.hydrostatic.T0 = reft;
  result.params.hydrostatic.lapse_rate = gamma;
  return result;
}

} // namespace haero
