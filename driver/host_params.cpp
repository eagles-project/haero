#include "host_params.hpp"

#include <sstream>

#include "ekat/ekat_assert.hpp"

namespace haero {
namespace driver {

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
  ss << tabstr << "ptop = " << ptop << " Pa\n";
  ss << tabstr << "tperiod = " << tperiod << " s\n";
  ss << tabstr << "qv0 = " << qv0 << " kg/kg\n";
  ss << tabstr << "qv1 = " << qv1 << " per m\n";
  return ss.str();
}

}  // namespace driver
}  // namespace haero
