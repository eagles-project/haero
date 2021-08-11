#ifndef HAERO_mam_nucleation_TEST_BRIDGE_HPP
#define HAERO_mam_nucleation_TEST_BRIDGE_HPP

#include "haero/haero.hpp"

using haero::Real;

extern "C" {

// Fortran subroutines that implement this process.

extern void nucleation_init_bridge();

extern void compute_tendencies_bridge(
    const Real &factor_bin_tern_ratenucl,
    const Real &adjust_factor_pbl_ratenucl, const Real &deltat,
    const Real &temp, const Real &pmid, const Real &aircon, const Real &zmid,
    const Real &pblh, const Real &relhum, const Real &uptkrate_h2so4,
    const Real &del_h2so4_gasprod, const Real &del_h2so4_aeruptk,
    const Real *qgas_cur, const Real *qgas_avg, const Real *qnum_cur,
    const Real *qaer_cur, const Real *qwtr_cur, Real &dndt_ait_f,
    Real &dmdt_ait_f, Real &dso4dt_ait_f, Real &dnh4dt_ait_f,
    Real &nclusterdt_f);

extern void ternary_nuc_merik2007_bridge(const Real t, const Real rh,
                                         const Real c2, const Real c3,
                                         Real &j_log, Real &ntot, Real &nacid,
                                         Real &namm, Real &r);
extern void binary_nuc_vehk2002_bridge(const Real temp, const Real rh,
                                       const Real so4vol, Real &ratenucl,
                                       Real &rateloge, Real &cnum_h2so4,
                                       Real &cnum_tot, Real &radius_cluster);

extern void pbl_nuc_wang2008_bridge(const Real adjust_factor_pbl_ratenucl,
                                    const Real so4vol,
                                    const int newnuc_method_flagaa,
                                    int &newnuc_method_flagaa2, Real &ratenucl,
                                    Real &rateloge, Real &cnum_tot,
                                    Real &cnum_h2so4, Real &cnum_nh3,
                                    Real &radius_cluster);

extern void mer07_veh02_nuc_mosaic_1box_bridge(
    const double factor_bin_tern_ratenucl,
    const double adjust_factor_pbl_ratenucl, const int newnuc_method_flagaa,
    const Real dtnuc, const Real temp_in, const Real rh_in, const Real press_in,
    const Real zm_in, const Real pblh_in, const Real qh2so4_cur,
    const Real qh2so4_avg, const Real qnh3_cur, const Real h2so4_uptkrate,
    const Real mw_so4a_host, const int nsize, const int maxd_asize,
    const Real *dplom_sect,  // array size maxd_asize
    const Real *dphim_sect,  // array size maxd_asize
    int &isize_nuc, Real &qnuma_del, Real &qso4a_del, Real &qnh4a_del,
    Real &qh2so4_del, Real &qnh3_del, Real &dens_nh4so4a, const int ldiagaa,
    Real *dnclusterdt = nullptr);

}  // extern "C"

#endif
