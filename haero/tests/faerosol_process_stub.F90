!> This module implements a stub for a prognostic process written in Fortran.
!> It demonstrates how to use Haero's Fortran data types in a way that allows
!> such a process to communicate with the underlying C++ machinery.
!>
!> This prognostic process transfers a cloudborne aerosol to an interstitial
!> aerosol using an exponential decay process. The process can work with any
!> number of modes/species. Gases are ignored. Modal number densities are
!> unaffected, as aerosol particles are simply transferred from one population
!> to the other.
module faerosol_process_stub

  use iso_c_binding, only: c_ptr
  use haero_precision, only: wp
  use haero, only: modal_aero_config, prognostics_t, atmosphere_t, diagnostics_t, &
                   tendencies_t, prognostics_from_c_ptr, atmosphere_from_c_ptr,&
                   diagnostics_from_c_ptr, tendencies_from_c_ptr

  implicit none
  private

  ! Process parameters (mark with save attribute)
  real(wp), save :: decay_rate ! Decay rate for cloudborne aerosols

  ! Module interface
  public :: process_stub_init, process_stub_run, process_stub_finalize

  ! Parameter setters
  public :: process_stub_set_integer_param, &
            process_stub_set_logical_param, &
            process_stub_set_real_param

  ! C function for obtaining decay rate.
  interface
    real(c_real) function process_stub_decay_rate() bind(c)
      use iso_c_binding, only: c_real
    end function
  end interface

contains

!> Performs initialization.
subroutine process_stub_init() bind(c)
  implicit none
end subroutine

!> Calls the update for the process, computing tendencies for each affected
!> prognostic variable.
subroutine process_stub_run(t, dt, progs, atm, diags, tends) bind(c)
  use iso_c_binding, only: c_ptr, c_f_pointer
  implicit none

  ! Arguments
  real(wp), value, intent(in) :: t     ! simulation time
  real(wp), value, intent(in) :: dt    ! simulation time step
  type(c_ptr), value, intent(in) :: progs ! prognostic variables
  type(c_ptr), value, intent(in) :: atm   ! atmosphere state variables
  type(c_ptr), value, intent(in) :: diags ! diagnostic variables
  type(c_ptr), value, intent(in) :: tends ! tendencies

  ! Fortran prognostics, diagnostics, tendencies types
  type(prognostics_t) :: prognostics
  type(atmosphere_t)  :: atmosphere
  type(diagnostics_t) :: diagnostics
  type(tendencies_t)  :: tendencies

  ! Other local variables.
  integer :: num_modes, p, k
  real(wp), pointer, dimension(:,:) :: q_c    ! cloudborne aerosol mix fracs
  real(wp), pointer, dimension(:,:) :: q_i    ! interstitial aerosol mix fracs
  real(wp), pointer, dimension(:,:) :: q_g    ! gas mole fracs
  real(wp), pointer, dimension(:,:) :: n      ! modal number densities
  real(wp), pointer, dimension(:,:) :: dqdt_c ! cloudborne aerosol tends
  real(wp), pointer, dimension(:,:) :: dqdt_i ! interstitial aerosol tends
  real(wp), pointer, dimension(:,:) :: dqdt_g ! gas mole frac tends
  real(wp), pointer, dimension(:,:) :: dndt   ! modal number density tends

  ! Get Fortran data types from our C pointers.
  prognostics = prognostics_from_c_ptr(progs)
  atmosphere = atmosphere_from_c_ptr(atm)
  diagnostics = diagnostics_from_c_ptr(diags)
  tendencies = tendencies_from_c_ptr(tends)

  ! Iterate over modes and compute aerosol mix fraction tendencies.
  num_modes = size(modal_aero_config%aerosol_modes)
  q_c => prognostics%cloud_aerosols()
  q_i => prognostics%interstitial_aerosols()
  dqdt_c => tendencies%cloud_aerosols()
  dqdt_i => tendencies%interstitial_aerosols()

  ! Cloudborne aerosols decay exponentially into interstitial aerosols.
  do p=1,modal_aero_config%num_aerosol_populations
    do k=1,prognostics%num_levels
      dqdt_c(k, p) =  decay_rate * q_c(k, p)
      dqdt_i(k, p) = -dqdt_c(k, p)
    end do
  end do

  ! Gas mix ratio tendencies are zero.
  q_g => prognostics%gases()
  dqdt_g => tendencies%gases()
  dqdt_g(:,:) = 0.0_wp

  ! Modal number mix ratio tendencies are zero.
  n => prognostics%interstitial_num_mix_ratios()
  dndt => tendencies%interstitial_num_mix_ratios()
  dndt(:,:) = 0.0_wp

end subroutine

!> Disposes of the process-specific data allocated in process_stub_init.
!> You don't need this if you haven't allocated any process-specific data.
subroutine process_stub_finalize() bind(c)
  implicit none

  ! Deallocate any process-specific resources
end subroutine

! Parameter setters -- we only support setting "decay_rate" to a real value.
subroutine process_stub_set_integer_param(name, val) bind(c)
  use iso_c_binding, only: c_ptr, c_int
  implicit none

  type(c_ptr), value, intent(in) :: name
  integer(c_int), value, intent(in) :: val
end subroutine

subroutine process_stub_set_logical_param(name, val) bind(c)
  use iso_c_binding, only: c_ptr, c_bool
  implicit none

  type(c_ptr), value, intent(in) :: name
  logical(c_bool), value, intent(in) :: val
end subroutine

subroutine process_stub_set_real_param(name, val) bind(c)
  use iso_c_binding, only: c_ptr, c_real
  use haero, only: c_to_f_string
  implicit none

  type(c_ptr), value, intent(in) :: name
  real(c_real), value, intent(in) :: val

  if (trim(c_to_f_string(name)) == "decay_rate") then
    decay_rate = val
  end if
end subroutine

end module


