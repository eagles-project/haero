!> This module implements a stub for a diagnostic process written in Fortran.
!> It demonstrates how to use Haero's Fortran data types in a way that allows
!> such a process to communicate with the underlying C++ machinery.
module diag_process_stub

  use iso_c_binding, only: c_ptr
  use haero, only: wp, model, prognostics_t, diagnostics_t, &
                   prognostics_from_c_ptr, diagnostics_from_c_ptr

  implicit none
  private

  ! Process parameters
  real(wp), parameter :: R = 8.314472 ! universal gas constant [J/K/mol]
  real(wp), parameter :: T0 = 283.15  ! constant temperature [K]

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
  integer :: num_modes, m, i, k
  real(wp), pointer, dimension(:,:,:) :: n ! modal number densities
  real(wp), pointer, dimension(:,:,:) :: p ! modal pressure

  ! Get Fortran data types from our C pointers.
  prognostics = prognostics_from_c_ptr(progs)
  diagnostics = diagnostics_from_c_ptr(diags)

  ! Diagnose modal pressure using ideal gas law
  n => prognostics%modal_num_densities()
  p => diagnostics%modal_var("pressure")
  num_modes = size(model%modes)
  do m=1,num_modes
    do i=1,model%num_columns
      do k=1,model%num_levels
        p(k, i, m) = n(k, i, m) * R * T
      end do
    end do
  end do

end subroutine

!> Disposes of the process-specific data allocated in diag_process_init.
!> You don't need this if you haven't allocated any process-specific data.
subroutine diag_stub_finalize() bind(c)
  implicit none
end subroutine

end module


