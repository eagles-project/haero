#ifndef HAERO_mam_gasaerexch_TEST_BRIDGE_HPP
#define HAERO_mam_gasaerexch_TEST_BRIDGE_HPP

#include "haero/haero.hpp"

using haero::Real;

extern "C" {

// Fortran subroutines that implement this process.

extern void gasaerexch_init_bridge();

extern void mam_gasaerexch_1subarea_1gas_nonvolatile_bridge(
    const Real &dt, const Real &qgas_netprod_otrproc, const int &n_mode,
    const Real *uptkaer, Real &qgas_cur, Real &qgas_avg, Real *qaer_cur);

extern void mam_soaexch_1subarea_bridge(
    const int max_gas, const int max_aer, const int iaer_pom, const int nsoa,
    const int npoa, const int npca, const Real pstd, const Real r_universal,
    const int lund, const Real dt, const Real temp, const Real pmid,
    const Real aiircon, const int n_mode, const int ntot_amode,
    const int max_mode, Real *qgas_cur, Real *qgas_avg, Real *qaer_cur,
    Real *qnum_cur, const Real *uptkaer, const int *mode_aging_optaa,
    const int *lptr2_soa_a_amode);

}  // extern "C"

#endif
