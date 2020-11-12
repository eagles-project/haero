!> This module implements a stub for a diagnostic process written in Fortran.
!> It demonstrates how to use Haero's Fortran data types in a way that allows
!> such a process to communicate with the underlying C++ machinery.
module diag_process_stub

  use iso_c_binding, only: c_ptr
  use haero, only: wp, model, prognostics_t, diagnostics_t

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
  real(wp), value, intent(in) :: t              ! simulation time
  type(c_ptr), intent(in)     :: progs          ! prognostic variables
  type(c_ptr), intent(inout)  :: diags          ! diagnostic variables

  ! Fortran prognostics and diagnostics types
  type(prognostics_t) :: prognostics
  type(diagnostics_t) :: diagnostics

  ! Emplace the C pointer into our Fortran types so we can use them.
  prognostics%ptr = progs
  diagnostics%ptr = diags

  ! Do stuff
end subroutine

!> Disposes of the process-specific data allocated in diag_process_init.
!> You don't need this if you haven't allocated any process-specific data.
subroutine diag_stub_finalize() bind(c)
  implicit none
end subroutine

end module


