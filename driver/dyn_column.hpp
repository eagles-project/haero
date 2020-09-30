#ifndef HAERO_DYN_COLUMN_HPP
#define HAERO_DYN_COLUMN_HPP

#include "haero/haero_config.hpp"
#include "driver/ncwriter.hpp"
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
    dpds("dp/ds", ncol, num_packs_int()),
    /// midpoint variables
    dpids("pseudodensity_dpi/ds", ncol, num_packs_lev()),
    thetav("virtual_potential_temperature", ncol, num_packs_lev()),
    p("pressure", ncol, num_packs_lev()),
    qv("water_vapor_mixing_ratio", ncol, num_packs_lev()),
    exner("exner", ncol, num_packs_lev()),
    dphids("dphi/ds", ncol, num_packs_lev()),
    /// surface variables
    psurf("surface_pressure", ncol),
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

      host_psurf = Kokkos::create_mirror_view(psurf);
  }

  /** @brief Return a netcdf writer (and open a new netcdf file)

    Defines all dimensions and variables required by the column dynamics model.
    Does not write any variable data except for constant coordinate variables.

    @param [in] filename
  */
  NcWriter write_new_ncdata(const std::string& filename) const;

  /** @brief Write current variable data to netCDF filename

    @param [in] writer
    @param [in] time_idx
  */
  void update_ncdata(NcWriter& writer, const size_t time_idx) const;

  /** @brief Initialize all column data to a hydrostatically balanced, stationary reference state
  using a list of z values to define interface heights.

    @param [in] z_vals
    @param [in] conds
  */
  void init_from_interface_heights(const std::vector<Real>& z_vals, const AtmosphericConditions& conds);

  /** @brief Initialize all column data to a hydrostatically balanced, stationary reference state
  using a list of p values to define interface hydrostatic pressures.

    @param [in] p_vals
    @param [in] conds
  */
  void init_from_interface_pressures(const std::vector<Real>& p_vals, const AtmosphericConditions& conds);

  /** @brief return a string with basic column state data for use with console output.

  */
  std::string info_string(const int& tab_lev=0) const ;

  inline int num_columns() const {return m_ncol;}

  /** @brief Returns the nondimensional s-coordinate associated with height z.

    @warning Must be called after either initialize method, since it requires z_top to be defined.

    @param [in] z height [m]
    @return s
  */
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
  view_2d dphids;

  /// surface variables
  view_1d psurf;

  /// constant views
  view_1d interface_scoord;
  view_1d interface_ds;
  view_1d level_scoord;
  view_1d level_ds;

  DynColumn() = delete;

  protected:
    Real m_ztop;
    Real m_ptop;
    Real m_ncol;
    Real m_psurf;

    host_view2d host_w;
    host_view2d host_phi;
    host_view2d host_pi;
    host_view2d host_mu;
    host_view2d host_dpds;

    host_view2d host_dpids;
    host_view2d host_thetav;
    host_view2d host_p;
    host_view2d host_qv;
    host_view2d host_exner;
    host_view2d host_dphids;

    host_view1d host_interface_scoord;
    host_view1d host_interface_ds;
    host_view1d host_level_scoord;
    host_view1d host_level_ds;

    host_view1d host_psurf;
};


} // namespace haero
#endif
