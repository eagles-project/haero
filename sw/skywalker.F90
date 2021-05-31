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
    ! Timestepping data
    real(c_real) :: dt, total_time
    ! atmospheric state parameters
    real(c_real) :: temperature, pressure, relative_humidity, height, &
                    hydrostatic_dp, planetary_boundary_layer_height

    ! Modal aerosol number concentrations [# aero molecules / kg air]
    real(c_real), dimension(:), pointer :: interstitial_number_concs, &
                                           cloud_number_concs
    ! Aerosol mass mixing ratios [kg aerosol / kg air], indexed by mode and species
    real(c_real), dimension(:, :), pointer :: interstitial_aero_mmrs, &
                                              cloud_aero_mmrs
    ! Gas mass mixing ratios [kg gas / kg air]
    real(c_real), dimension(:), pointer :: gas_mmrs
   end type input_data_t

   type :: output_data_t
    ! C pointer
    type(c_ptr) :: ptr
    ! Modal aerosol number concentrations [# aero molecules / kg air]
    real(c_real), dimension(:), pointer :: interstitial_number_concs, &
                                           cloud_number_concs
    ! Aerosol mass mixing ratios [kg aerosol / kg air], indexed by mode and species
    real(c_real), dimension(:, :), pointer :: interstitial_aero_mmrs, &
                                              cloud_aero_mmrs
    ! Gas mass mixing ratios [kg gas / kg air]
    real(c_real), dimension(:), pointer :: gas_mmrs
   end type output_data_t

   ! This type represents an ensemble and its corresponding input and output
   ! data for each member.
   type :: ensemble_t
     ! C pointer
     type(c_ptr) :: ptr
     ! The number of members in the ensemble
     integer :: size
     ! Number of aerosol modes and populations, and number of gases
     integer(c_int) :: num_modes, num_populations, num_gases
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
      type(c_ptr), value, intent(in) :: aerosols_per_mode
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
    use iso_c_binding, only: c_int, c_ptr, c_associated, c_loc
    implicit none

    character(len=*), intent(in) :: aerosol_config
    character(len=*), intent(in) :: filename
    type(ensemble_t) :: ensemble

    integer(c_int), dimension(:), pointer :: mode_array_sizes
    integer :: i, m, s, p, max_mode_size
    real(c_real), dimension(:), pointer :: int_aero_data, cld_aero_data

    ensemble%ptr = sw_load_ensemble(f_to_c_string(aerosol_config), f_to_c_string(filename))
    if (.not. c_associated(ensemble%ptr)) then
      print *, "Could not load a ", aerosol_config, " from ", filename
      stop
    end if

    ! Size up the ensemble.
    ensemble%size = sw_ensemble_size(ensemble%ptr)
    allocate(ensemble%inputs(ensemble%size))
    allocate(ensemble%outputs(ensemble%size))

    ! Extract metadata.
    call sw_ensemble_get_array_sizes(ensemble%ptr, ensemble%num_modes, &
                                     ensemble%num_populations, &
                                     ensemble%num_gases)
    allocate(mode_array_sizes(ensemble%num_modes))
    call sw_ensemble_get_modal_aerosol_sizes(ensemble%ptr, c_loc(mode_array_sizes))
    max_mode_size = 0
    do m = 1, ensemble%num_modes
      max_mode_size = max(max_mode_size, mode_array_sizes(m))
    end do
    allocate(int_aero_data(max_mode_size))
    allocate(cld_aero_data(max_mode_size))

    ! Allocate input and output aerosol and gas arrays for the ensemble,
    ! and extract input data
    do i = 1, ensemble%size
      allocate(ensemble%inputs(i)%interstitial_number_concs(ensemble%num_modes))
      allocate(ensemble%inputs(i)%cloud_number_concs(ensemble%num_modes))
      allocate(ensemble%inputs(i)%interstitial_aero_mmrs(ensemble%num_modes, max_mode_size))
      allocate(ensemble%inputs(i)%cloud_aero_mmrs(ensemble%num_modes, max_mode_size))
      allocate(ensemble%inputs(i)%gas_mmrs(ensemble%num_gases))
      allocate(ensemble%outputs(i)%interstitial_aero_mmrs(ensemble%num_modes, max_mode_size))
      allocate(ensemble%outputs(i)%cloud_aero_mmrs(ensemble%num_modes, max_mode_size))
      allocate(ensemble%outputs(i)%gas_mmrs(ensemble%num_gases))

      ! Aerosol data
      call sw_input_get_aerosols(ensemble%inputs(i)%ptr, &
        c_loc(ensemble%inputs(i)%interstitial_number_concs), &
        c_loc(ensemble%inputs(i)%cloud_number_concs), &
        c_loc(int_aero_data), c_loc(cld_aero_data))
      p = 1 ! population index
      do m = 1, ensemble%num_modes
        do s = 1, mode_array_sizes(m)
          ensemble%inputs(i)%interstitial_aero_mmrs(m, s) = int_aero_data(p)
          ensemble%inputs(i)%cloud_aero_mmrs(m, s) = cld_aero_data(p)
          p = p + 1
        end do
      end do

      ! Gas data
      call sw_input_get_gases(ensemble%inputs(i)%ptr, &
        c_loc(ensemble%inputs(i)%gas_mmrs))
    end do

    ! Clean up.
    deallocate(int_aero_data, cld_aero_data, mode_array_sizes)
  end function

  ! Writes a Python module containing all input and output for the given
  ! ensemble to a file with the given name.
  subroutine e_write_py_module(ensemble, filename)
    use iso_c_binding, only: c_ptr
    implicit none

    class(ensemble_t), intent(in) :: ensemble
    character(len=*), intent(in) :: filename

    integer i, m, s, p, max_mode_size
    integer(c_int), dimension(:), pointer :: mode_array_sizes
    real(c_real), dimension(:), pointer :: int_aero_data, cld_aero_data

    ! Copy output data into place.
    max_mode_size = size(ensemble%outputs(1)%interstitial_aero_mmrs, 2)
    allocate(int_aero_data(max_mode_size))
    allocate(cld_aero_data(max_mode_size))

    ! Aerosol data
    allocate(mode_array_sizes(ensemble%num_modes))
    call sw_ensemble_get_modal_aerosol_sizes(ensemble%ptr, c_loc(mode_array_sizes))
    do i = 1, ensemble%size
      p = 1 ! population index
      do m = 1, ensemble%num_modes
        do s = 1, mode_array_sizes(m)
          int_aero_data(p) = ensemble%outputs(i)%interstitial_aero_mmrs(m, s)
          cld_aero_data(p) = ensemble%outputs(i)%cloud_aero_mmrs(m, s)
          p = p + 1
        end do
      end do
      call sw_input_set_aerosols(ensemble%outputs(i)%ptr, &
        c_loc(ensemble%outputs(i)%interstitial_number_concs), &
        c_loc(ensemble%outputs(i)%cloud_number_concs), &
        c_loc(int_aero_data), c_loc(cld_aero_data))

      ! Gas data
      call sw_input_set_gases(ensemble%inputs(i)%ptr, &
        c_loc(ensemble%inputs(i)%gas_mmrs))
    end do

    ! Generate the Python module.
    call sw_ensemble_write_py_module(ensemble%ptr, f_to_c_string(filename))

    ! Clean up.
    deallocate(mode_array_sizes, int_aero_data, cld_aero_data)
  end subroutine

end module
