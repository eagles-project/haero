!> This module implements a stub for a diagnostic process written in Fortran.
!> It demonstrates how to use Haero's Fortran data types in a way that allows
!> such a process to communicate with the underlying C++ machinery.
module diag_process_stub

  use iso_c_binding, only: c_ptr
  use haero, only: wp => real_t, model_t, prognostics_t, diagnostics_t

  implicit none
  private

  public :: diag_process_init, &
            diag_process_update, &
            diag_process_destroy

  !> This type stores parameters for data specific to this diagnostic process.
  type :: diag_process_data_t
    !> How many modes are in this aerosol model?
    integer :: nmodes
    !> This parameter stores a floating point value in the correct precision.
    real(wp) :: f
    !> Here's an array that holds data for a vertical column.
    real(wp), dimension(:), allocatable :: column_data
  end type

contains

!> Creates and returns any process-specific data. You don't need this function
!> if your process doesn't need to store its own data.
type(c_ptr) function diag_process_init(model) bind(c)
  ! We use the c_loc function to return a pointer to process-specific dataa
  use iso_c_binding, only: c_loc

  implicit none

  ! Argument: the aerosol model
  type(Model), intent(in) :: model

  ! Local variables
  integer :: nlevels

  ! We store process-specific data here.
  type(diag_process_data_t), pointer, target :: process_data

  ! Extract some information from the model.
  nlevels = model%num_levels()

  ! Initialize the process-specific information.
  allocate(process_data)
  process_data%nmodes = model%num_modes()     ! number of modes in the model
  process_data%f = 3.8                        ! bogus parameter value
  allocate(process_data%column_data(nlevels)) ! allocate column data array
  do i=1,nlevels
    process_data%column_data(i) = i
  end do

  ! Return a C pointer to the process-specific data.
  diag_process_init = c_loc(process)
end function

!> Disposes of the process-specific data allocated in diag_process_init.
!> You don't need this if you haven't allocated any process-specific data.
subroutine diag_process_destroy(c_process_data) bind(c)
  use iso_c_binding, only :: c_ptr, c_f_pointer
  implicit none

  ! Argument: a C pointer to the process-specific data that was allocated
  ! and returned in diag_process_init. If you returned c_null_ptr in that
  ! subroutine, there's no need to do anything with p here.
  type(c_ptr), intent(inout) :: c_process_data

  ! Fortran pointer to process-specific data
  type(diag_process_data_t), pointer :: process_data

  ! Convert the C pointer to a Fortran pointer.
  call c_f_pointer(c_process_data, process_data)

  ! Deallocate the column data array and the process data itself.
  deallocate(process%column_data)
  deallocate(process)
end subroutine

subroutine diag_process_update(c_process_data, model, t, progs, diags) bind(c)
  use iso_c_binding, only :: c_ptr, c_f_pointer
  implicit none

  ! Arguments
  type(c_ptr),       intent(in)      :: c_process_data ! C pointer to process-ѕpecific data
  type(model_t),       intent(in)    :: model          ! aerosol model
  real(wp), value,   intent(in)      :: t              ! simulation time
  type(prognostics_t), intent(in)    :: progs          ! prognostic variables
  type(diagnostics_t), intent(inout) :: diags          ! diagnostic variables

  ! Fortran pointer to process-specific data
  type(diag_process_data_t), pointer :: process_data

  ! Convert the C pointer to a Fortran pointer to access process-specific data.
  call c_f_pointer(c_process_data, process_data)

end subroutine

end module


