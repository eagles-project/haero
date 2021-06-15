!> This module is a bridge between the MAMCalcsizeFProcess C++ class and the
!> mam_calcsize Fortran module.
module mam_calcsize_bridge

  implicit none
  private

  ! Module functions
  public :: mam_calcsize_init, &
       run_bridge, &
            mam_calcsize_finalize, &
            compute_diameter_bridge

contains

! Initializes the prognostic process.
subroutine mam_calcsize_init() bind(c)
  use haero, only: model
  use mam_calcsize, only: init
  implicit none

  call init(model)
end subroutine

! Runs the prognostic process, computing tendencies
subroutine run_bridge(t, dt, progs, atm, diags, tends) bind(c)
  use iso_c_binding, only: c_ptr, c_f_pointer
  use haero_precision, only: wp
  use haero, only: model, &
                   prognostics_t, atmosphere_t, diagnostics_t, tendencies_t, &
                   prognostics_from_c_ptr, atmosphere_from_c_ptr, &
                   diagnostics_from_c_ptr, tendencies_from_c_ptr
  use mam_calcsize, only: run
  implicit none

  ! Arguments
  real(wp), value, intent(in) :: t     ! simulation time
  real(wp), value, intent(in) :: dt    ! simulation time step
  type(c_ptr), value, intent(in) :: progs ! prognostic variables
  type(c_ptr), value, intent(in) :: atm   ! atmospheric state
  type(c_ptr), value, intent(in) :: diags ! diagnostic variables
  type(c_ptr), value, intent(in) :: tends ! tendencies

  ! Fortran prognostics, atmosphere, diagnostics, tendencies types
  type(prognostics_t) :: prognostics
  type(atmosphere_t)  :: atmosphere
  type(diagnostics_t) :: diagnostics
  type(tendencies_t)  :: tendencies

  ! Get Fortran data types from our C pointers.
  prognostics = prognostics_from_c_ptr(progs)
  atmosphere = atmosphere_from_c_ptr(atm)
  diagnostics = diagnostics_from_c_ptr(diags)
  tendencies = tendencies_from_c_ptr(tends)

  ! Call the actual subroutine.
  call run(model, t, dt, prognostics, atmosphere, diagnostics, tendencies)
end subroutine run_bridge

pure function compute_diameter_bridge(vol2num) result(diameter) bind(c)

  use haero_precision, only: wp
  use mam_calcsize, only:compute_diameter

  implicit none

  real(wp), value, intent(in) :: vol2num

  real(wp) :: diameter

  !diameter = vol2num + 1
  diameter = compute_diameter(vol2num)

end function compute_diameter_bridge

! Finalizes the prognostic process
subroutine mam_calcsize_finalize() bind(c)
  use haero, only: model
  use mam_calcsize, only: finalize
  implicit none

  call finalize(model)
end subroutine mam_calcsize_finalize

end module
