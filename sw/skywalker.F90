!> This module contains data structures that allow Fortran modules to access
!> skywalker's ensemble input and output data.
module skywalker

  use iso_c_binding
  use haero_precision, only: wp

  implicit none

  private

  public :: input_data_t, output_data_t, ensemble_t, load_ensemble

  !> This Fortran type is the equivalent of the C++ InputData struct.
  type :: input_data_t
    ! C pointer
    type(c_ptr) :: ptr

   end type input_data_t

   type :: output_data_t
    ! C pointer
    type(c_ptr) :: ptr

   end type output_data_t

   ! This type represents an ensemble and its corresponding input and output
   ! data for each member.
   type :: ensemble_t
     ! C pointer
     type(c_ptr) :: ptr
     ! The number of members in the ensemble
     integer :: size
     ! An array of input data for every member of the ensemble
     type(input_data_t), dimension(:), allocatable :: inputs
     ! An array of output data for every member of the ensemble
     type(output_data_t), dimension(:), allocatable :: outputs

     contains
       ! Writes a Python module containing input/output data to a file
       procedure :: write_py_module => e_write_py_module
   end type ensemble_t

  interface

    type(c_ptr) function sw_load_ensemble(aerosol_config, filename) bind(c)
      use iso_c_binding, only: c_ptr
      type(c_ptr), value, intent(in) :: aerosol_config, filename
    end function

    integer(c_int) function sw_ensemble_size(ensemble) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: ensemble
    end function

    subroutine sw_ensemble_get_array_sizes(ensemble, num_modes, num_pops, num_gases) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: ensemble
      integer(c_int), intent(out) :: num_modes, num_pops, num_gases
    end subroutine

    subroutine sw_ensemble_get_modal_aerosol_sizes(ensemble, aerosols_per_mode) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: ensemble
      integer(c_int), dimension(:), intent(out) :: aerosols_per_mode
    end subroutine

    type(c_ptr) function sw_ensemble_input(ensemble, i) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: ensemble
      integer(c_int), intent(in) :: i
    end function

    subroutine sw_input_get_timestepping(input, dt, total_time) bind(c)
      use iso_c_binding, only: c_ptr, c_real
      type(c_ptr), value, intent(in) :: input
      real(c_real), intent(out) :: dt, total_time
    end subroutine

    subroutine sw_input_get_atmosphere(input, temperature, pressure, rh, height,&
                                       dp, pblh) bind(c)
      use iso_c_binding, only: c_ptr, c_real
      type(c_ptr), value, intent(in) :: input
      real(c_real), intent(out) :: temperature, pressure, rh, height, dp, pblh
    end subroutine

    subroutine sw_input_get_aerosols(input, int_num_concs, cld_num_concs, &
                                     int_aero_mmrs, cld_aero_mmrs) bind(c)
      use iso_c_binding, only: c_ptr, c_real
      type(c_ptr), value, intent(in) :: input
      type(c_ptr), value, intent(in) :: int_num_concs, cld_num_concs,&
                                        int_aero_mmrs, cld_aero_mmrs
    end subroutine

    subroutine sw_input_get_gases(input, gas_mmrs) bind(c)
      use iso_c_binding, only: c_ptr, c_real
      type(c_ptr), value, intent(in) :: input
      type(c_ptr), value, intent(in) :: gas_mmrs
    end subroutine

    type(c_ptr) function sw_ensemble_output(ensemble, i) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: ensemble
      integer(c_int), intent(in) :: i
    end function

    subroutine sw_output_set_aerosols(output, int_num_concs, cld_num_concs, &
                                      int_aero_mmrs, cld_aero_mmrs) bind(c)
      use iso_c_binding, only: c_ptr, c_real
      type(c_ptr), value, intent(in) :: output
      type(c_ptr), value, intent(in) :: int_num_concs, cld_num_concs,&
                                        int_aero_mmrs, cld_aero_mmrs
    end subroutine

    subroutine sw_output_set_gases(output, gas_mmrs) bind(c)
      use iso_c_binding, only: c_ptr, c_real
      type(c_ptr), value, intent(in) :: output
      type(c_ptr), value, intent(in) :: gas_mmrs
    end subroutine

    subroutine sw_ensemble_write_py_module(ensemble, filename) bind(c)
      use iso_c_binding, only: c_ptr
      type(c_ptr), value, intent(in) :: ensemble
      type(c_ptr), value, intent(in) :: filename
    end subroutine

    subroutine sw_ensemble_free(ensemble) bind(c)
      use iso_c_binding, only: c_ptr
      type(c_ptr), value, intent(in) :: ensemble
    end subroutine

  end interface

contains

  !> This helper function converts the given Fortran string to a C string.
  function f_to_c_string(f_string) result(c_string)
    use, intrinsic :: iso_c_binding
    implicit none
    character(len=*), target :: f_string
    character(len=:), pointer :: f_ptr
    type(c_ptr) :: c_string

    interface
        function new_c_string(f_str_ptr, f_str_len) bind (c) result(c_string)
        use, intrinsic :: iso_c_binding
            type(c_ptr), value :: f_str_ptr
            integer(c_int), value :: f_str_len
            type(c_ptr) :: c_string
        end function new_c_string
    end interface

    f_ptr => f_string
    c_string = new_c_string(c_loc(f_ptr), len(f_string))
  end function f_to_c_string

  ! Given an aerosol configuration string and a skywalker input file, fetch an
  ! ensemble's worth of input data.
  function load_ensemble(aerosol_config, filename) result(ensemble)
    use iso_c_binding, only: c_int, c_ptr, c_associated
    implicit none

    character(len=*), intent(in) :: aerosol_config
    character(len=*), intent(in) :: filename
    type(ensemble_t) :: ensemble

    ensemble%ptr = sw_load_ensemble(f_to_c_string(aerosol_config), f_to_c_string(filename))
    if (.not. c_associated(ensemble%ptr)) then
      print *, "Could not load a ", aerosol_config, " from ", filename
      stop
    end if

    ! Size up the ensemble and extract data.
    ensemble%size = sw_ensemble_size(ensemble%ptr)
    allocate(ensemble%inputs(ensemble%size))
    allocate(ensemble%outputs(ensemble%size))

  end function

  subroutine e_write_py_module(ensemble, filename)
    use iso_c_binding, only: c_ptr
    implicit none

    class(ensemble_t), intent(in) :: ensemble
    character(len=*), intent(in) :: filename

    call sw_ensemble_write_py_module(ensemble%ptr, f_to_c_string(filename))
  end subroutine

end module
