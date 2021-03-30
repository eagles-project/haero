!> This module is a bridge between the MAMNucleationFProcess C++ class and the
!> mam_nucleation Fortran module.
module mam_nucleation_test_bridge

  implicit none
  private

  ! Module functions
  public :: ternary_nuc_merik2007_bridge
  public :: binary_nuc_vehk2002_bridge

contains

subroutine ternary_nuc_merik2007_bridge(t, rh, c2, c3, j_log, ntot, nacid, namm, r) bind(c)
  use haero, only: wp
  use mam_nucleation, only: ternary_nuc_merik2007
  implicit none

  ! Arguments
  real(wp), value, intent(in) :: t, rh, c2, c3
  real(wp), intent(out) :: j_log, ntot, nacid, namm, r

  ! Call the actual subroutine.
  call ternary_nuc_merik2007(t, rh, c2, c3, j_log, ntot, nacid, namm, r)
end subroutine

subroutine binary_nuc_vehk2002_bridge(temp, rh, so4vol, ratenucl, rateloge, cnum_h2so4, cnum_tot, radius_cluster) bind(c)
  use haero, only: wp
  use mam_nucleation, only: binary_nuc_vehk2002
  implicit none

  ! Arguments
  real(wp), value, intent(in) :: temp, rh, so4vol
  real(wp), intent(out) :: ratenucl, rateloge, cnum_h2so4, cnum_tot, radius_cluster

  ! Call the actual subroutine.
  call binary_nuc_vehk2002(temp, rh, so4vol, ratenucl, rateloge, cnum_h2so4, cnum_tot, radius_cluster)
end subroutine

subroutine pbl_nuc_wang2008_bridge(factor_pbl_ratenucl, so4vol, flagaa, flagaa2, ratenucl, rateloge, &
  cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster) bind(c)
  use iso_c_binding, only: c_int
  use haero, only: wp
  use mam_nucleation, only: pbl_nuc_wang2008
  use mam_nucleation, only: adjust_factor_pbl_ratenucl
  implicit none

  ! Arguments
  real(wp), value, intent(in) :: factor_pbl_ratenucl
  real(wp), value, intent(in) :: so4vol
  integer(c_int), value, intent(in) :: flagaa
  integer(c_int), intent(inout) :: flagaa2
  real(wp), intent(inout) :: ratenucl, rateloge, cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster

  ! Call the actual subroutine.
  ! But first set this public value on the module that the function will use.
  adjust_factor_pbl_ratenucl = factor_pbl_ratenucl
  call pbl_nuc_wang2008(so4vol, flagaa, flagaa2, ratenucl, rateloge, &
    cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster)
end subroutine

end module

