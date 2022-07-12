#ifndef HAERO_MAM4_NUCLEATION_IMPL_HPP
#define HAERO_MAM4_NUCLEATION_IMPL_HPP

#include <haero/atmosphere.hpp>
#include <haero/mam4/aero_config.hpp>

namespace haero {
namespace mam4 {

/// @class NucleationImpl
/// This class implements MAM4's nucleation parameterization. Its structure is
/// defined by the usage of the impl_ member in the AeroProcess class in
/// ../aero_process.hpp.
class NucleationImpl {
  static const int nait       = static_cast<int>(ModeIndex::Aitken);
  static const int igas_h2so4 = static_cast<int>(GasId::H2SO4);
  static const int igas_nh3   = static_cast<int>(GasId::NH3);

  static constexpr Real avogadro = Constants::avogadro;
  static constexpr Real mw_h2so4 = Constants::molec_weight_h2so4;
  static constexpr Real pi       = Constants::pi;
  static constexpr Real rgas     = Constants::r_gas;

  // min h2so4 vapor for nuc calcs = 4.0e-16 mol/mol-air ~= 1.0e4 molecules/cm3,
  static constexpr Real qh2so4_cutoff = 4.0e-16;
  static constexpr Real ln_nuc_rate_cutoff = -13.82;

  // "Host parameters.
  Real dens_so4a_host, mw_nh4a_host, mw_so4a_host;
 public:

  using Pack = PackType;

  // name -- unique name of the process implemented by this class
  const char* name() const { return "MAM4 nucleation"; }

  // init -- initializes the implementation with MAM4's configuration
  void init(const AeroConfig& config) {
    // set "host" variables.
    dens_so4a_host = 0.0;
    mw_nh4a_host = 0.0;
    mw_so4a_host = 0.0;
  }

  // validate -- validates the given atmospheric state and prognostics against
  // assumptions made by this implementation, returning true if the states are
  // valid, false if not
  KOKKOS_INLINE_FUNCTION
  bool validate(const AeroConfig& config,
                const TeamType& team,
                const Atmosphere& atm,
                const Prognostics& progs) const {

    return atm.quantities_nonnegative(team) &&
           progs.quantities_nonnegative(team);
  }

  // compute_tendencies -- computes tendencies and updates diagnostics
  // NOTE: that both diags and tends are const below--this means their views
  // NOTE: are fixed, but the data in those views is allowed to vary.
  KOKKOS_INLINE_FUNCTION
  void compute_tendencies(const AeroConfig& config,
                          const TeamType& team, Real t, Real dt,
                          const Atmosphere& atm,
                          const Prognostics& progs,
                          const Diagnostics& diags,
                          const Tendencies& tends) const {
    const int nk = PackInfo::num_packs(atm.num_levels());
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nk),
      KOKKOS_LAMBDA(int k) {
        // extract column data at level k
        Pack temp = atm.temperature(k);
        Pack pmid = atm.pressure(k);
        Pack aircon = 0;
        Pack zmid = 0;
        Real pblh = atm.planetary_boundary_height;
        Pack relhum = 0;
        Pack uptkrate_so4 = 0;
        Pack del_h2so4_gasprod = 0;
        Pack del_h2so4_aeruptk = 0;
        Pack qgas_cur[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        Pack qgas_avg[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        Pack qnum_cur[4] = {0, 0, 0, 0};
        Pack qaer_cur[4][7] = {
          {0, 0, 0, 0, 0, 0, 0},
          {0, 0, 0, 0, 0, 0, 0},
          {0, 0, 0, 0, 0, 0, 0},
          {0, 0, 0, 0, 0, 0, 0},
        };
        Pack qwtr_cur[4] = {0, 0, 0, 0};
        Pack dndt_ait, dmdt_ait, dso4dt_ait, dnh4dt_ait, dnclusterdt;

        // compute tendencies at this level
        compute_tendencies_(dt, temp, pmid, aircon, zmid, pblh, relhum,
                            uptkrate_so4, del_h2so4_gasprod, del_h2so4_aeruptk,
                            qgas_cur, qgas_avg, qnum_cur, qaer_cur,
                            qwtr_cur, dndt_ait, dmdt_ait, dso4dt_ait,
                            dnh4dt_ait, dnclusterdt);
      });
  }

  // This function computes relevant tendencies at a single vertical level. It
  // was ported directly from the compute_tendencies subroutine in the
  // modal_aero_newnuc module from the MAM4 box model.
  void compute_tendencies_(Real deltat,
    const Pack& temp, const Pack& pmid, const Pack& aircon, const Pack& zmid,
    Real pblh, const Pack& relhum, const Pack& uptkrate_h2so4,
    const Pack& del_h2so4_gasprod, const Pack& del_h2so4_aeruptk,
    const Pack qgas_cur[13], const Pack qgas_avg[13],
    const Pack qnum_cur[4], const Pack qaer_cur[4][7],
    const Pack qwtr_cur[4],
    Pack& dndt_ait, Pack& dmdt_ait, Pack& dso4dt_ait, Pack& dnh4dt_ait,
    Pack& dnclusterdt) const {

    int newnuc_method_actual, pbl_nuc_wang2008_actual;

    constexpr int nsize = 1;
    Real dplom_mode[nsize], dphim_mode[nsize];
    int isize_group;

    Pack cair;  // air density
    Pack so4vol, nh3ppt;

    Pack radius_cluster;  // radius of newly formed cluster, in nm
    Pack rateloge;        // ln(J)
    Pack cnum_h2so4;
    Pack cnum_nh3;

    Pack mass1p;
    Pack mass1p_aithi, mass1p_aitlo;
    Pack qh2so4_cur, qh2so4_avg, qh2so4_del;
    Pack qnh3_cur, qnh3_del, qnh4a_del;
    Pack qnuma_del;
    Pack qso4a_del;
    Pack relhumnn;
    Pack tmpa, tmpb, tmpc;
    Pack tmp_q2, tmp_q3;
    Pack tmp_frso4, tmp_uptkrate;
    Pack dens_nh4so4a;

    dndt_ait = 0;
    dmdt_ait = 0;
    dnh4dt_ait = 0;
    dso4dt_ait = 0;
    dnclusterdt = 0;

    //==================================
    // limit RH to between 0.1% and 99%
    //==================================
    relhumnn = max( 0.01, min( 0.99, relhum ) );

    //=========================================================================
    // prepare h2so4 mixing ratio and condensation rate that will be passed to
    // the nucleation parameterization
    //=========================================================================
    qh2so4_cur = qgas_cur[igas_h2so4];

    // E3SM: use qh2so4_avg and first-order loss rate calculated in
    // mam_gasaerexch_1subarea
    qh2so4_avg = qgas_avg[igas_h2so4];
    tmp_uptkrate = uptkrate_h2so4;

    if (qh2so4_avg <= qh2so4_cutoff) { // qh2so4_avg very low. assume no
                                       // nucleation will happen
      // Diagnose so4 and nh4 tendencies and exit
      dso4dt_ait = dmdt_ait*tmp_frso4/mw_so4a_host;
      dnh4dt_ait = dmdt_ait*(1.0 - tmp_frso4)/mw_nh4a_host;
      return;
    }

    if (igas_nh3 > 0) {
      qnh3_cur = max( 0.0, qgas_cur(igas_nh3) )
    } else {
      qnh3_cur = 0.0;
    }

    // unit conversion for gas concentrations:
    // calculate h2so4 in molecules/cm3 and nh3 in ppt
    cair   = pmid/(temp*rgas);
    so4vol = qh2so4_avg * cair * avogadro * 1.0e-6;
    nh3ppt = qnh3_cur * 1.0e12;

    //=======================================================================
    // call routine to get nucleation rate in terms of new cluster formation
    // rate (#/m3/s)
    //=======================================================================
    if (newnuc_method_user_choice != 0) {

      // Hui Wan's note from code refactoring in July 2021:
      // Subroutine mer07_veh02_wang08_nuc_1box provides
      //  - dnclusterdt (unit: #/m3/s): new cluster formation rate
      //  - rateloge (unit: ln(#/cm3/s)): logarithm of new cluster formation rate
      //  - cnum_h2so4, cnum_nh3: number of of h2so4 or nh3 molecules per cluster
      //  - radius_cluster (unit: nm): radius of new cluster
      // Output variables rateloge, cnum_h2so4, cnum_nh3, and radius_cluster
      // are used below in the calculation of cluster "growth". I chose to keep
      // these variable names the same as in the old subroutine mer07_veh02_nuc_mosaic_1box
      // to facilitate comparison.

      mer07_veh02_wang08_nuc_1box(
        newnuc_method_user_choice, newnuc_method_actual,              // in, out
        pbl_nuc_wang2008_user_choice, pbl_nuc_wang2008_actual,        // in, out
        ln_nuc_rate_cutoff,                                           // in
        adjust_factor_bin_tern_ratenucl, adjust_factor_pbl_ratenucl,  // in
        pi, so4vol, nh3ppt, temp, relhumnn, zmid, pblh,               // in
        dnclusterdt, rateloge, cnum_h2so4, cnum_nh3, radius_cluster); // out

    } else {
      rateloge    = ln_nuc_rate_cutoff;
      dnclusterdt = 0.;
      newnuc_method_actual = 0;
    }

    //======================================================================
    // "Grow" the newly formed clusters to size in the smallest bin/mode of
    // the host model
    //======================================================================
    qnuma_del = 0.0;
    qso4a_del = 0.0;
    qnh4a_del = 0.0;
    qh2so4_del = 0.0;
    qnh3_del = 0.0;

    // dry-diameter limits for "grown" new particles
    dplom_mode(1) = exp( 0.67*log(dgnumlo_aer(nait))
                       + 0.33*log(dgnum_aer(nait)) );
    dphim_mode(1) = dgnumhi_aer(nait);

    //----------------------------------------------------------------
    // Only do the cluster growth calculation when nucleation rate is
    // appreciable
    //----------------------------------------------------------------
    if (rateloge > ln_nuc_rate_cutoff ) {

      // mass1p_... = mass (kg) of so4 & nh4 in a single particle of diameter ...
      // (assuming same dry density for so4 & nh4)
      // mass1p_aitlo - dp = dplom_mode(1);
      // mass1p_aithi - dp = dphim_mode(1);

      tmpa = dens_so4a_host*pi/6.0;
      mass1p_aitlo = tmpa*cube(dplom_mode(1));
      mass1p_aithi = tmpa*cube(dphim_mode(1));

      // Cluster growth
      newnuc_cluster_growth(
        dnclusterdt, cnum_h2so4, cnum_nh3, radius_cluster,
        dplom_mode, dphim_mode, nsize,
        deltat, temp, relhumnn, cair,
        accom_coef_h2so4, mw_so4a,mw_so4a_host, mw_nh4a, avogadro, pi,
        qnh3_cur, qh2so4_cur, so4vol, tmp_uptkrate,
        isize_group, dens_nh4so4a,
        qh2so4_del, qnh3_del, qso4a_del, qnh4a_del, qnuma_del);
    } // nucleation rate is appreciable

    //=====================================
    // Deriving mass mixing ratio tendency
    //=====================================

    // convert qnuma_del from (#/mol-air) to (#/kmol-air)
    qnuma_del *= 1.0e3;

    // number nuc rate (#/kmol-air/s) from number nuc amt
    dndt_ait = qnuma_del/deltat;

    // fraction of mass nuc going to so4
    tmpa = qso4a_del*mw_so4a_host;
    if (igas_nh3 > 0) {
      tmpb = tmpa + qnh4a_del*mw_nh4a_host;
      tmp_frso4 = max( tmpa, 1.0e-35)/max( tmpb, 1.0e-35);
    } else {
      tmpb = tmpa;
      tmp_frso4 = 1.0;
    }

    // mass nuc rate (kg/kmol-air/s) from mass nuc amts
    dmdt_ait = max( 0.0, (tmpb/deltat) );

    //=====================================================
    // Various adjustments to keep the solution reasonable
    //=====================================================
    if (dndt_ait < 1.0e2) {
      // ignore newnuc if number rate < 100 #/kmol-air/s ~= 0.3 #/mg-air/d
      dndt_ait = 0.0;
      dmdt_ait = 0.0;
    } else {

      // mirage2 code checked for complete h2so4 depletion here,
      // but this is now done in mer07_veh02_nuc_mosaic_1box
      mass1p = dmdt_ait/dndt_ait;

      // apply particle size constraints
      if (mass1p < mass1p_aitlo) {
        // reduce dndt to increase new particle size
        dndt_ait = dmdt_ait/mass1p_aitlo;
      } else if (mass1p > mass1p_aithi) {
        // reduce dmdt to decrease new particle size
        dmdt_ait = dndt_ait*mass1p_aithi;
      }
    }

    // *** apply adjustment factor to avoid unrealistically high
    // aitken number concentrations in mid and upper troposphere
    dndt_ait = dndt_ait * newnuc_adjust_factor_dnaitdt;
    dmdt_ait = dmdt_ait * newnuc_adjust_factor_dnaitdt;

    //=================================
    // Diagnose so4 and nh4 tendencies
    //=================================
    // dso4dt_ait, dnh4dt_ait are (kmol/kmol-air/s)
    dso4dt_ait = dmdt_ait*tmp_frso4/mw_so4a_host;
    dnh4dt_ait = dmdt_ait*(1.0 - tmp_frso4)/mw_nh4a_host;
  }
};

} // namespace mam4
} // namespace haero

#endif
