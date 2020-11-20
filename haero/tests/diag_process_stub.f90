!> This module implements a stub for a diagnostic process written in Fortran.
!> It demonstrates how to use Haero's Fortran data types in a way that allows
!> such a process to communicate with the underlying C++ machinery.
module diag_process_stub

  use iso_c_binding, only: c_ptr
  use haero, only: wp, model, prognostics_t, diagnostics_t, &
                   prognostics_from_c_ptr, diagnostics_from_c_ptr

  implicit none
  private

  public :: diag_stub_init, &
            diag_stub_update, &
            diag_stub_finalize

contains

!> Performs initialization
subroutine diag_stub_init() bind(c)
  implicit none
end subroutine

!> Calls the update for the process.
subroutine diag_stub_update(t, progs, diags) bind(c)
  use iso_c_binding, only: c_ptr, c_f_pointer
  implicit none

  ! Arguments
  real(wp), value, intent(in) :: t     ! simulation time
  type(c_ptr), intent(in)     :: progs ! prognostic variables
  type(c_ptr), intent(inout)  :: diags ! diagnostic variables

  ! Fortran prognostics and diagnostics types
  type(prognostics_t) :: prognostics
  type(diagnostics_t) :: diagnostics

  ! Other local variables.
  integer :: num_modes, m
  real(wp), pointer, dimension(:,:,:) :: q_a ! (interstitial) aerosol mix fracs
  real(wp), pointer, dimension(:,:,:) :: q_g ! gas mole fracs
  real(wp), pointer, dimension(:,:,:) :: f_a ! (interstitial) aerosol frequencies
  real(wp), pointer, dimension(:,:,:) :: f_g ! gas frequencies

  ! Get Fortran data types from our C pointers.
  prognostics = prognostics_from_c_ptr(progs)
  diagnostics = diagnostics_from_c_ptr(diags)

  ! Iterate over modes and diagnose aerosol oscillation frequencies.
  num_modes = size(model%modes)
  f_a = diagnostics%modal_var("aerosol_frequencies")
  do m=1,num_modes
    q_a = prognostics%interstitial_aerosols(m)
  end do

  ! Diagnose gas mole fraction oscillation frequency.
  q_g = prognostics%gas_mole_fractions()
  f_g = diagnostics%var("gas_frequencies")
end subroutine

!> Disposes of the process-specific data allocated in diag_process_init.
!> You don't need this if you haven't allocated any process-specific data.
subroutine diag_stub_finalize() bind(c)
  implicit none
end subroutine

end module


