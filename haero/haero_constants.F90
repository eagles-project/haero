!> This module contains physical constants available for use in Fortran
!> implementations of aerosol processes.
module haero_constants

  use haero, only: wp

  implicit none

  public

  !> Pi, to however much precision we need [-]
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)

  !> Universal molar gas constant [J/K/mol]
  real(wp), parameter :: R_gas = 8.31446261815324_wp

  !> Avogadro's constant [#/mol]
  real(wp), parameter :: Avogadro = 6.02214076e23_wp

end module

