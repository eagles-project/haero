// This file is autogenerated. DO NOT EDIT DIRECTLY.
#ifndef HAERO_mam_nucleation_TEST_BRIDGE_HPP
#define HAERO_mam_nucleation_TEST_BRIDGE_HPP

#include "haero/haero.hpp"

using haero::Real;

extern "C" {

// Fortran subroutines that implement this process.
extern void ternary_nuc_merik2007_bridge(const Real t, const Real rh, const Real c2, const Real c3, 
 Real &j_log, Real &ntot, Real &nacid, Real &namm, Real &r);
extern void binary_nuc_vehk2002_bridge(const Real temp, const Real rh, const Real so4vol, 
  Real &ratenucl, Real &rateloge, Real &cnum_h2so4, Real &cnum_tot, Real &radius_cluster);
extern void pbl_nuc_wang2008_bridge(const Real adjust_factor_pbl_ratenucl, const Real so4vol, const int newnuc_method_flagaa, 
  int & newnuc_method_flagaa2, Real & ratenucl, Real & rateloge, 
  Real & cnum_tot, Real & cnum_h2so4, Real & cnum_nh3, Real & radius_cluster );

} // extern "C"

#endif
