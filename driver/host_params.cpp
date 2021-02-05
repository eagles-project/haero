#include "host_params.hpp"
#include "ekat/ekat_assert.hpp"
#include "haero/floating_point.hpp"
#include "haero/utils.hpp"
#include <sstream>

namespace haero {
namespace driver {

using fp_helper = FloatingPoint<Real>;

AtmosphericConditions::AtmosphericConditions(const Real Tv0_, const Real Gammav_,
 const Real w0_, const int ztop_,  const int tperiod_, const Real qv0_, 
 const Real qv1_ ) {

  /// check valid input
  EKAT_REQUIRE_MSG(fp_helper::in_bounds(Tv0_, 273, 323),
    "unexpected T0, check units = K");
  Tv0 = Tv0_;
  
  EKAT_REQUIRE_MSG(fp_helper::in_bounds(w0_, 0, 10), "unexpected w0, check units = m/s");
  w0 = w0_;
  
  EKAT_REQUIRE_MSG(fp_helper::in_bounds(Gammav_, 0, 0.02),
    "unexpected lapse rate, check units = K/m");
  Gammav = Gammav_;
  
  EKAT_REQUIRE_MSG(fp_helper::in_bounds(ztop_, 3E3,40E3),
    "unexpected model top, check units = m");
  ztop = ztop_;
  
  EKAT_REQUIRE_MSG(tperiod_>0, "nonnegative oscillation period required.");
  tperiod = tperiod_;
  
  EKAT_REQUIRE_MSG(fp_helper::in_bounds(qv0_, 0, 0.1),
    "unexpected water vapor mixing ratio; check units = kg/kg");
  qv0 = qv0_;
  
  EKAT_REQUIRE_MSG(qv1 >= 0, "nonnegative decay rate required.");
  qv1 = qv1_;

}

std::string AtmosphericConditions::info_string(const int tab_level) const {
  std::string tabstr = indent_string(tab_level);
  std::ostringstream ss;
  ss << tabstr << "AtmosphericConditions info:\n";
  tabstr += "\t";
  ss << tabstr << "alpha_v = " << AtmosphericConditions::alpha_v << '\n';
  ss << tabstr << "kappa = " << AtmosphericConditions::kappa << '\n';
  ss << tabstr << "Tv0 = " << Tv0 << " K\n";
  ss << tabstr << "Gammav = " << Gammav << " K/m\n";
  ss << tabstr << "w0 = " << w0 << " m/s\n";
  ss << tabstr << "ztop = " << ztop << " m\n";
  ss << tabstr << "tperiod = " << tperiod << " s\n";
  ss << tabstr << "qv0 = " << qv0 << " kg/kg\n";
  ss << tabstr << "qv1 = " << qv1 << " per m\n";
  return ss.str();
}

} // namespace driver
} // namespace haero
