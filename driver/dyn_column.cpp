#include "column_base.hpp"
#include "dyn_column.hpp"
#include "ncwriter_impl.hpp"
#include "haero/utils.hpp"
#include "haero/physical_constants.hpp"
#include "ekat/util/ekat_units.hpp"
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

  for (int col_idx=0; col_idx<m_ncol; ++col_idx) {
    host_psurf(col_idx) = conds.params.hydrostatic.p0;

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
              host_interface_scoord(pack_info::pack_idx(num_levels()-1))[pack_info::vec_idx(num_levels()-1)];
            const Real p_nphalf = m_psurf;
            const Real p_n = host_p(col_idx, pack_info::pack_idx(num_levels()-1))[pack_info::vec_idx(num_levels()-1)];

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
  Kokkos::deep_copy(w, host_w);
  Kokkos::deep_copy(phi, host_phi);
  Kokkos::deep_copy(pi, host_pi);
  Kokkos::deep_copy(mu, host_mu);
  Kokkos::deep_copy(dpds, host_dpds);
  Kokkos::deep_copy(dpids, host_dpids);
  Kokkos::deep_copy(thetav, host_thetav);
  Kokkos::deep_copy(p, host_p);
  Kokkos::deep_copy(qv, host_qv);
  Kokkos::deep_copy(exner, host_exner);
  Kokkos::deep_copy(dphids, host_dphids);
  Kokkos::deep_copy(interface_scoord, host_interface_scoord);
  Kokkos::deep_copy(interface_ds, host_interface_ds);
  Kokkos::deep_copy(level_scoord, host_level_scoord);
  Kokkos::deep_copy(level_ds, host_level_ds);
  Kokkos::deep_copy(psurf, host_psurf);
}

NcWriter DynColumn::write_new_ncdata(const std::string& filename) const {
  NcWriter writer(filename);

  using att_type = NcWriter::text_att_type;
  using var_atts = std::vector<att_type>;

  writer.add_column_dim(w.extent(0));
  writer.add_level_dims(num_levels());
  writer.define_time_var();

  /// Interface variables
  const var_atts w_atts = {std::make_pair("cf_long_name", "upward_air_velocity"),
    std::make_pair("short_name", "w")};
  const auto w_units = ekat::units::m/ekat::units::s;

  const var_atts phi_atts = {std::make_pair("cf_long_name", "geopotential"),
    std::make_pair("short_name","phi")};
  const auto phi_units = ekat::units::pow(ekat::units::m,2)*ekat::units::pow(ekat::units::s,-2);

  const var_atts pi_atts = {std::make_pair("cf_long_name", "null"),
                            std::make_pair("haero_long_name", "hydrostatic_air_pressure"),
                            std::make_pair("short_name", "pi")};
  const auto pi_units = ekat::units::Pa;

  const var_atts mu_atts = {std::make_pair("cf_long_name", "null"),
            std::make_pair("haero_long_name", "full_pressure_to_hydrostatic_pressure_ratio"),
            std::make_pair("short_name", "mu")};
  const auto mu_units = ekat::units::Units::nondimensional();

  const var_atts dpds_atts = {std::make_pair("cf_long_name", "null"),
            std::make_pair("haero_long_name", "change_in_air_pressure_at_level_interface"),
            std::make_pair("short_name", "dp/ds")};
  const auto dpds_units = ekat::units::Pa;

  writer.define_interface_var("vertical_velocity", w_units, w, w_atts);
  writer.define_interface_var("geopotential", phi_units, phi, phi_atts);
  writer.define_interface_var("hydrostatic_pressure", pi_units, pi, pi_atts);
  writer.define_interface_var("mu", mu_units, mu, mu_atts);
  writer.define_interface_var("dpds", dpds_units, dpds, dpds_atts);

  /// Level variables
  const var_atts dpids_atts = {std::make_pair("cf_long_name", "null"),
      std::make_pair("haero_long_name", "pseudodensity_of_air"),
      std::make_pair("short_name", "dpi/ds")};
  const auto dpids_units = ekat::units::Pa;

  const var_atts thetav_atts = {std::make_pair("cf_long_name", "null"),
    std::make_pair("haero_long_name", "virtual_potential_temperature"),
    std::make_pair("short_name", "theta_v")};
  const auto thetav_units = ekat::units::K;

  const var_atts p_atts = {std::make_pair("cf_long_name", "air_pressure"),
    std::make_pair("amip_short_name", "plev"),
    std::make_pair("short_name", "p")};
  const auto p_units = ekat::units::Pa;

  const var_atts qv_atts = {std::make_pair("cf_long_name", "humidity_mixing_ratio"),
    std::make_pair("haero_long_name", "water_vapor_mass_mixing_ratio"),
    std::make_pair("short_name", "q_v")};
  const auto qv_units = ekat::units::kg / ekat::units::kg;

  const var_atts exner_atts = {std::make_pair("cf_long_name", "dimensionless_exner_function"),
    std::make_pair("short_name", "Pi")};
  const auto exner_units = ekat::units::Units::nondimensional();

  const var_atts dphids_atts = {std::make_pair("cf_long_name", "null"),
    std::make_pair("haero_long_name", "change_in_geopotential_at_level_midpoint"),
    std::make_pair("short_name", "dphi/ds")};

  writer.define_level_var("pseudodensity", dpids_units, dpids, dpids_atts);
  writer.define_level_var("virtual_potential_temperature", thetav_units, thetav,
    thetav_atts);
  writer.define_level_var("pressure", p_units, p, p_atts);
  writer.define_level_var("water_vapor_mixing_ratio", qv_units, qv, qv_atts);
  writer.define_level_var("exner", exner_units, exner, exner_atts);
  writer.define_level_var("dphids", phi_units, dphids, dphids_atts);

  /// surface variables
  const var_atts ps_atts = {std::make_pair("cf_long_name", "surface_air_pressure"),
    std::make_pair("short_name", "psurf"), std::make_pair("amip_short_name", "ps")};
  const auto ps_units = ekat::units::Pa;
  writer.define_time_dependent_scalar_var("surface_pressure", ps_units, ps_atts);

  /// Scalar variables
  const var_atts ptop_atts = {
  std::make_pair("cf_long_name","air_pressure_at_top_of_atmosphere_model"),
  std::make_pair("short_name", "p_top")};
  const auto ptop_units = ekat::units::Pa;

  writer.define_scalar_var("p_top", ptop_units, ptop_atts, m_ptop);

  /// Coordinate variables
  const var_atts scoord_atts = {
    std::make_pair("cf_long_name","atmosphere_hybrid_sigma_pressure_coordinate"),
    std::make_pair("short_name", "scoord")};
  const auto scoord_units = ekat::units::Units::nondimensional();
  const auto sinterface_vec = view1d_to_vector(interface_scoord, num_levels()+1);
  const auto slevel_vec = view1d_to_vector(level_scoord, num_levels());
  const auto ds_levels = view1d_to_vector(level_ds, num_levels());
  const auto ds_interfaces = view1d_to_vector(interface_ds, num_levels()+1);
  writer.define_const_1dvar("scoord_interfaces", scoord_units, sinterface_vec,
    scoord_atts);
  writer.define_const_1dvar("scoord_levels", scoord_units, slevel_vec, scoord_atts);
  writer.define_const_1dvar("ds_levels", scoord_units, ds_levels);
  writer.define_const_1dvar("ds_interfaces", scoord_units, ds_interfaces);

  return writer;
}

void DynColumn::update_ncdata(NcWriter& writer, const size_t time_idx) const {

  const auto ps = view1d_to_vector(psurf, m_ncol);
  writer.add_time_dependent_scalar_values("surface_pressure", time_idx, ps);

  for (size_t col_idx = 0; col_idx < m_ncol; ++col_idx) {
    writer.add_variable_data("vertical_velocity", time_idx, col_idx, w);
    writer.add_variable_data("geopotential", time_idx, col_idx, phi);
    writer.add_variable_data("hydrostatic_pressure", time_idx, col_idx, pi);
    writer.add_variable_data("mu", time_idx, col_idx, mu);
    writer.add_variable_data("dpds", time_idx, col_idx, dpds);

    writer.add_variable_data("pseudodensity", time_idx, col_idx, dpids);
    writer.add_variable_data("virtual_potential_temperature", time_idx, col_idx, thetav);
    writer.add_variable_data("pressure", time_idx, col_idx, p);
    writer.add_variable_data("water_vapor_mixing_ratio", time_idx, col_idx, qv);
    writer.add_variable_data("exner", time_idx, col_idx, exner);
    writer.add_variable_data("dphids", time_idx, col_idx, dphids);
  }

}

} // namespace haero
