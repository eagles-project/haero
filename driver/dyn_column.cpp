#include "column_base.hpp"
#include "dyn_column.hpp"
#include "haero/utils.hpp"
#include "haero/physical_constants.hpp"
#include "kokkos/Kokkos_Core.hpp"
#include <sstream>

namespace haero {

std::string DynColumn::info_string(const int& tab_lev) const {
  std::ostringstream ss;
  std::string tabstr = indent_string(tab_lev);
  ss << tabstr << "haero dynamics column info:\n";
  tabstr += "\t";
  ss << tabstr << "nlev = " << num_levels() << '\n';
  ss << tabstr << "npack_mid = " << num_packs_lev() << '\n';
  ss << tabstr << "npack_interface = " << num_packs_int() << '\n';
  return ss.str();
}

void DynColumn::init_from_interface_heights(const std::vector<Real>& z_vals,
  const AtmosphericConditions& conds) {

  /// check input
  EKAT_REQUIRE_MSG(vector_is_monotone(z_vals), "input height interface values must be monotone.");
  EKAT_REQUIRE_MSG(z_vals.size() == num_levels() + 1,
    "z-values are stored on level interfaces; num_levels() + 1 values are required.");

  /// copy input in case we need to reverse the z_vals so that they go from model top to the surface
  std::vector<Real> interface_heights(z_vals);
  std::vector<Real> level_heights(num_levels());
  const bool increasing = (z_vals[1] > z_vals[0]);
  if (increasing) std::reverse(interface_heights.begin(), interface_heights.end());
  /// set model top
  m_ztop = interface_heights[0];
  m_ptop = hydrostatic_pressure_at_height(m_ztop, conds);

  for (int col_idx=0; col_idx<m_ncol; ++col_idx) {
    /// set independent interface variables
    for (int pack_idx=0; pack_idx<num_packs_int(); ++pack_idx) {
      for (int vec_idx=0; vec_idx<HAERO_PACK_SIZE; ++vec_idx) {
        const int array_idx = pack_info(pack_idx, vec_idx);
        if (array_idx < num_levels()+1) { // skip padding
          host_w(col_idx, pack_idx)[vec_idx] = 0;
          host_phi(col_idx, pack_idx)[vec_idx] = gravity_m_per_s2*interface_heights[array_idx];
          host_pi(col_idx, pack_idx)[vec_idx] =
            hydrostatic_pressure_at_height(interface_heights[array_idx], conds);
          host_mu(col_idx, pack_idx)[vec_idx] = 1;
          host_interface_scoord(pack_idx)[vec_idx] = s_coord(interface_heights[i]);
        }
      }
    }
    /// set independent level variables
    for (int pack_idx=0; pack_idx<num_packs_lev(); ++pack_idx) {
      for (int vec_idx=0; vec_idx<HAERO_PACK_SIZE; ++vec_idx) {
        const int array_idx = pack_info(pack_idx, vec_idx);
        if (array_idx < num_levels()) { // skip padding
          const Real zmid = 0.5*(interface_heights[i] + interface_heights[i+1]);
          const Real Tv = virtual_temperature(zmid, conds);
          const Real qv = water_vapor_mixing_ratio(qv, conds);
          const Real p = hydrostatic_pressure_at_height(zmid, conds);
          host_p(pack_idx)[vec_idx] = p;
          host_theta_v(pack_idx)[vec_idx] = potential_temperature(Tv, p, conds);
          host_qv(pack_idx)[vec_idx] = qv;
          host_exner(pack_idx)[vec_idx] = exner_function(p, conds);
          host_level_scoord(pack_idx)[vec_idx] = s_coord(zmid);
        }
      }
    }
  }
}


}

} // namespace haero
