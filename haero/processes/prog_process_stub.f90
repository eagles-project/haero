!> This module implements a stub for a prognostic process written in Fortran.
!> It demonstrates how to use Haero's Fortran data types in a way that allows
!> such a process to communicate with the underlying C++ machinery.
module prog_process_stub

  use iso_c_binding, only: c_ptr
  use haero, only: wp, model, prognostics_t, diagnostics_t, tendencies_t

  implicit none
  private

  public :: prog_stub_init, &
            prog_stub_run, &
            prog_stub_finalize

contains

!> Performs initialization.
subroutine prog_stub_init() bind(c)
  implicit none
end subroutine

!> Calls the update for the process.
subroutine prog_stub_run(t, dt, progs, diags, tends) bind(c)
  use iso_c_binding, only: c_ptr, c_f_pointer
  implicit none

  ! Arguments
  real(wp), value, intent(in) :: t              ! simulation time
  real(wp), value, intent(in) :: dt             ! simulation time step
  type(c_ptr), intent(in)     :: progs          ! prognostic variables
  type(c_ptr), intent(inout)  :: diags          ! diagnostic variables
  type(c_ptr), intent(inout)  :: tends          ! tendencies

  ! Fortran prognostics and diagnostics types
  type(prognostics_t) :: prognostics
  type(diagnostics_t) :: diagnostics
  type(tendencies_t)  :: tendencies

  ! Emplace the C pointer into our Fortran types so we can use them.
  prognostics%ptr = progs
  diagnostics%ptr = diags
  tendencies%ptr = tends

  ! Do stuff
end subroutine

!> Disposes of the process-specific data allocated in prog_process_init.
!> You don't need this if you haven't allocated any process-specific data.
subroutine prog_stub_finalize() bind(c)
  implicit none
end subroutine

end module


