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
  use iso_c_binding, only: c_ptr
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
  use iso_c_binding, only: c_ptr
  use haero, only: wp
  use mam_nucleation, only: binary_nuc_vehk2002
  implicit none

  ! Arguments
  real(wp), value, intent(in) :: temp, rh, so4vol
  real(wp), intent(out) :: ratenucl, rateloge, cnum_h2so4, cnum_tot, radius_cluster

  ! Call the actual subroutine.
  call binary_nuc_vehk2002(temp, rh, so4vol, ratenucl, rateloge, cnum_h2so4, cnum_tot, radius_cluster)
end subroutine

end module

