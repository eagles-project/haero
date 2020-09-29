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
  /// set boundary data
  m_ztop = interface_heights[0];
  m_ptop = hydrostatic_pressure_at_height(m_ztop, conds);
  m_psurf = conds.params.hydrostatic.p0;

  for (int col_idx=0; col_idx<m_ncol; ++col_idx) {
    /// set independent interface variables
    for (int pack_idx=0; pack_idx<num_packs_int(); ++pack_idx) {
      for (int vec_idx=0; vec_idx<HAERO_PACK_SIZE; ++vec_idx) {
        const int array_idx = pack_info::array_idx(pack_idx, vec_idx);
        if (array_idx < num_levels()+1) { // skip padding
          host_w(col_idx, pack_idx)[vec_idx] = 0;
          host_phi(col_idx, pack_idx)[vec_idx] = gravity_m_per_s2*interface_heights[array_idx];
          host_pi(col_idx, pack_idx)[vec_idx] =
            hydrostatic_pressure_at_height(interface_heights[array_idx], conds);
          host_mu(col_idx, pack_idx)[vec_idx] = 1;
          if (col_idx == 0) {
            host_interface_scoord(pack_idx)[vec_idx] = s_coord(interface_heights[array_idx]);
          }
        }
      }
    }
    /// set independent level variables
    for (int pack_idx=0; pack_idx<num_packs_lev(); ++pack_idx) {
      for (int vec_idx=0; vec_idx<HAERO_PACK_SIZE; ++vec_idx) {
        const int array_idx = pack_info::array_idx(pack_idx, vec_idx);
        if (array_idx < num_levels()) { // skip padding
          const Real zmid = 0.5*(interface_heights[array_idx] + interface_heights[array_idx+1]);
          level_heights[array_idx] = zmid;
          const Real Tv = virtual_temperature(zmid, conds);
          const Real qv = water_vapor_mixing_ratio(zmid, conds);
          const Real p = hydrostatic_pressure_at_height(zmid, conds);
          host_p(col_idx, pack_idx)[vec_idx] = p;
          host_thetav(col_idx, pack_idx)[vec_idx] = potential_temperature(Tv, p, conds);
          host_qv(col_idx, pack_idx)[vec_idx] = qv;
          host_exner(col_idx, pack_idx)[vec_idx] = exner_function(p, conds);
          if (col_idx == 0) {
            host_level_scoord(pack_idx)[vec_idx] = s_coord(zmid);
          }
        }
      }
    }
    /// set level-dependent interface variables
    for (int pack_idx=0; pack_idx<num_packs_int(); ++pack_idx) {
      for (int vec_idx=0; vec_idx<HAERO_PACK_SIZE; ++vec_idx) {
        const int array_idx = pack_info::array_idx(pack_idx, vec_idx);
        if (array_idx < num_levels()+1) {//skip padding
          /** Below, array_idx maps to
              interface id and level id comments that correspond
              to indices in Taylor et. al., Figure 1.
          */
          if (array_idx == 0) { // model top
            // my interface id = 1/2 = "half"
            const Real ds_half = host_interface_scoord(pack_info::pack_idx(1))[pack_info::vec_idx(1)] -
                                 host_interface_scoord(pack_info::pack_idx(0))[pack_info::vec_idx(0)];
            const Real p_half = m_ptop;
            const Real p_1 = host_p(col_idx,pack_info::pack_idx(0))[pack_info::vec_idx(0)];

            /// Taylor et. al. eqn. (33)
            host_dpds(col_idx, pack_idx)[vec_idx] = 2*(p_1 - p_half)/ds_half;
            if (col_idx == 0) {
              host_interface_ds(pack_idx)[vec_idx] = ds_half;
            }
          }
          else if (array_idx == num_levels()) { // model surface
            // my interface id = n + 1/2, "nphalf", where "p" => "plus"
            const Real ds_nphalf = 1 -
              host_interface_scoord(pack_info::pack_idx(num_levels()))[pack_info::vec_idx(num_levels())];
            const Real p_nphalf = m_psurf;
            const Real p_n = host_p(col_idx, pack_info::pack_idx(num_levels()))[pack_info::vec_idx(num_levels())];

            /// Taylor et. al. eqn. (33)
            host_dpds(col_idx, pack_idx)[vec_idx] = 2*(p_nphalf - p_n)/ds_nphalf;
            if (col_idx == 0) {
              host_interface_ds(pack_idx)[vec_idx] = ds_nphalf;
            }
          }
          else { // general case
            // my interface id = i+1/2, "iphalf"
            const int level_i = array_idx-1; // level i
            const int level_ip1 = array_idx; // level i+1
            const int ipack = pack_info::pack_idx(level_i);
            const int ivec  = pack_info::vec_idx(level_i);
            const int ip1pack = pack_info::pack_idx(level_ip1);
            const int ip1vec = pack_info::vec_idx(level_ip1);

            const Real p_i   = host_p(col_idx, ipack)[ivec];
            const Real p_ip1 = host_p(col_idx, ip1pack)[ip1vec];
            const Real ds_iphalf = host_level_scoord(ip1pack)[ip1vec] -
                            host_level_scoord(ipack)[ivec];

            /// Taylor et. al. eqn. (32)
            host_dpds(col_idx, pack_idx)[vec_idx] = (p_ip1 - p_i) / (ds_iphalf);
            if (col_idx == 0) {
              host_interface_ds(pack_idx)[vec_idx] = ds_iphalf;
            }
          }
        }
      }
      /// set interface-dependent level variables
      for (int pack_idx=0; pack_idx<num_packs_lev(); ++pack_idx) {
        for (int vec_idx=0; vec_idx<HAERO_PACK_SIZE; ++vec_idx) {
          const int array_idx = pack_info::array_idx(pack_idx, vec_idx);
          if (array_idx < num_levels()) { // skip padding
            // my level id = i = "array_idx"
            const int interface_imhalf = array_idx;
            const int interface_iphalf = array_idx+1;
            const int imhalf_pack = pack_info::pack_idx(interface_imhalf);
            const int imhalf_vec  = pack_info::vec_idx(interface_imhalf);
            const int iphalf_pack = pack_info::pack_idx(interface_iphalf);
            const int iphalf_vec = pack_info::vec_idx(interface_iphalf);

            const Real ds_i = host_interface_scoord(iphalf_pack)[iphalf_vec] -
                              host_interface_scoord(imhalf_pack)[imhalf_vec];
            const Real pi_mhalf = host_pi(col_idx, imhalf_pack)[imhalf_vec];
            const Real pi_phalf = host_pi(col_idx, iphalf_pack)[iphalf_vec];
            const Real phi_mhalf = host_phi(col_idx, imhalf_pack)[imhalf_vec];
            const Real phi_phalf = host_phi(col_idx, iphalf_pack)[iphalf_vec];

            /// Taylor et. al. eqn. (32)
            host_dpids(col_idx, pack_idx)[vec_idx] = (pi_phalf - pi_mhalf)/ds_i;
            host_dphids(col_idx, pack_idx)[vec_idx] = (phi_phalf - phi_mhalf)/ds_i;
            if (col_idx == 0) {
              host_level_ds(pack_idx)[vec_idx] = ds_i;
            }
          }
        }
      }
    }
  }
}


} // namespace haero
