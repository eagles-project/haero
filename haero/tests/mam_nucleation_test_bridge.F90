!> This module is a bridge between the MAMNucleationFProcess C++ class and the
!> mam_nucleation Fortran module.
module mam_nucleation_test_bridge

  implicit none
  private

  ! Module functions
  public :: ternary_nuc_merik2007_bridge

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

end module

