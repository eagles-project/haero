#include "column_base.hpp"
#include "dyn_column.hpp"
#include "haero/utils.hpp"
#include <sstream>

namespace haero {

void DynColumn::initialize(const AtmosphericConditions& conds) {
  const Real gamma = conds.params.hydrostatic.lapse_rate;
  const Real T0 = conds.params.hydrostatic.T0;
  const Real p0 = conds.params.hydrostatic.p0;

}

std::string DynColumn::info_string(const int& tab_lev) const {
  std::ostringstream ss;
  std::string tabstr = indent_string(tab_lev);
  ss << tabstr << "haero dynamics column info:\n";
  tabstr += "\t";
  ss << tabstr << "nlev = " << nlev << '\n';
  ss << tabstr << "npack_mid = " << npack_mid << '\n';
  ss << tabstr << "npack_interface = " << npack_interface << '\n';
  return ss.str();
}

} // namespace haero
