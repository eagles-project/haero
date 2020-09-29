#ifndef HAERO_DYN_COLUMN_HPP
#define HAERO_DYN_COLUMN_HPP

#include "haero/haero_config.hpp"
#include "column_base.hpp"
#include "atmosphere.hpp"
#include <map>
#include <string>
#include <vector>

namespace haero {

class DynColumn : public ColumnBase {
  public:
  typedef typename view_2d::HostMirror host_view2d;
  typedef typename view_1d::HostMirror host_view1d;

  DynColumn(const int ncol, const int nl) : ColumnBase(nl),
    /// interface variables
    w("vertical_velocity", ncol, num_packs_int()),
    phi("geopotential", ncol, num_packs_int()),
    pi("hydrostatic_presssure", ncol, num_packs_int()),
    mu("mu", ncol, num_packs_int()),
    dpds("dpds", ncol, num_packs_int()),
    /// midpoint variables
    dpids("pseudodensity", ncol, num_packs_lev()),
    thetav("virtual_potential_temperature", ncol, num_packs_lev()),
    p("pressure", ncol, num_packs_lev()),
    qv("water_vapor_mixing_ratio", ncol, num_packs_lev()),
    exner("exner", ncol, num_packs_lev()),
    dphids("dphids", ncol, num_packs_lev()),
    /// coordinate variables
    interface_scoord("interface_scoord", num_packs_int()),
    interface_ds("inteface_ds", num_packs_int()),
    level_scoord("level_scoord", num_packs_lev()),
    level_ds("level_ds", num_packs_lev()),
    /// scalars
    m_ztop(0), m_ptop(0), m_ncol(ncol) {

      host_w = Kokkos::create_mirror_view(w);
      host_phi = Kokkos::create_mirror_view(phi);
      host_pi = Kokkos::create_mirror_view(pi);
      host_mu = Kokkos::create_mirror_view(mu);
      host_dpds = Kokkos::create_mirror_view(dpds);

      host_dpids = Kokkos::create_mirror_view(dpids);
      host_thetav = Kokkos::create_mirror_view(thetav);
      host_p = Kokkos::create_mirror_view(p);
      host_qv = Kokkos::create_mirror_view(qv);
      host_exner = Kokkos::create_mirror_view(exner);
      host_dphids = Kokkos::create_mirror_view(dphids);

      host_level_scoord = Kokkos::create_mirror_view(level_scoord);
      host_level_ds = Kokkos::create_mirror_view(level_ds);
      host_interface_scoord = Kokkos::create_mirror_view(interface_scoord);
      host_interface_ds = Kokkos::create_mirror_view(interface_ds);
    }

  void init_from_interface_heights(const std::vector<Real>& z_vals, const AtmosphericConditions& conds);

  void init_from_interface_pressures(const std::vector<Real>& p_vals, const AtmosphericConditions& conds);

  std::string info_string(const int& tab_lev=0) const ;

  inline int num_columns() const {return m_ncol;}

  inline Real s_coord(const Real z) const {return (m_ztop - z)/m_ztop;}

  /// level interface variables
  view_2d w;
  view_2d phi;
  view_2d pi;
  view_2d mu;
  view_2d dpds;

  /// level midpoint variables
  view_2d dpids;
  view_2d thetav;
  view_2d p;
  view_2d qv;
  view_2d exner;

  view_1d level_scoord;
  view_1d interface_scoord;
  view_1d level_ds;
  view_1d interface_ds;

  DynColumn() = delete;

  protected:
    Real m_ztop;
    Real m_ptop;
    Real m_ncol;

    host_view2d host_w;
    host_view2d host_phi;
    host_view2d host_dpids;
    host_view2d host_thetav;
    host_view2d host_p;
    host_view2d host_qv;
    host_view2d host_pi;
    host_view2d host_mu;
    host_view2d host_exner;

    host_view1d host_level_scoord;
    host_view1d host_level_ds;
    host_view1d host_interface_scoord;
    host_view1d host_interface_ds;
};


} // namespace haero
#endif
