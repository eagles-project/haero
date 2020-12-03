!> This module implements a stub for a prognostic process written in Fortran.
!> It demonstrates how to use Haero's Fortran data types in a way that allows
!> such a process to communicate with the underlying C++ machinery.
!>
!> The prognostic process evolves a collection of sinusoidally-varying
!> signals in aerosols and gases. There is an accompanying diagnostic process
!> that diagnoses the frequencies of these sinusoids. All species in a mode
!> have the same frequency, as do all gas species.
module prog_process_stub

  use iso_c_binding, only: c_ptr
  use haero, only: wp, model, prognostics_t, diagnostics_t, tendencies_t, &
                   prognostics_from_c_ptr, diagnostics_from_c_ptr, &
                   tendencies_from_c_ptr

  implicit none
  private

  ! Process parameters
  real(wp) :: tau ! e-folding rate for cloudborne aerosols

  public :: prog_stub_init, &
            prog_stub_run, &
            prog_stub_finalize

contains

!> Performs initialization.
subroutine prog_stub_init() bind(c)
  implicit none

  integer :: num_modes, i

  ! Initialize process parameters.
  tau = 1.0_wp
end subroutine

!> Calls the update for the process, computing tendencies for each affected
!> prognostic variable.
subroutine prog_stub_run(t, dt, progs, diags, tends) bind(c)
  use iso_c_binding, only: c_ptr, c_f_pointer
  implicit none

  ! Arguments
  real(wp), value, intent(in) :: t     ! simulation time
  real(wp), value, intent(in) :: dt    ! simulation time step
  type(c_ptr), intent(in)     :: progs ! prognostic variables
  type(c_ptr), intent(inout)  :: diags ! diagnostic variables
  type(c_ptr), intent(inout)  :: tends ! tendencies

  ! Fortran prognostics and diagnostics types
  type(prognostics_t) :: prognostics
  type(diagnostics_t) :: diagnostics
  type(tendencies_t)  :: tendencies

  ! Other local variables.
  integer :: num_modes, m, i, k
  real(wp), pointer, dimension(:,:,:) :: q_i    ! interstitial aerosol mix fracs
  real(wp), pointer, dimension(:,:,:) :: q_c    ! cloudborne aerosol mix fracs
  real(wp), pointer, dimension(:,:,:) :: q_g    ! gas mole fracs
  real(wp), pointer, dimension(:,:,:) :: n      ! modal number densities
  real(wp), pointer, dimension(:,:,:) :: dqdt_i ! interstitial aerosol tends
  real(wp), pointer, dimension(:,:,:) :: dqdt_c ! cloudborne aerosol tends
  real(wp), pointer, dimension(:,:,:) :: dqdt_g ! gas mole frac tends
  real(wp), pointer, dimension(:,:,:) :: dndt   ! modal number density tends

  ! Get Fortran data types from our C pointers.
  prognostics = prognostics_from_c_ptr(progs)
  diagnostics = diagnostics_from_c_ptr(diags)
  tendencies = tendencies_from_c_ptr(tends)

  ! Iterate over modes and compute aerosol mix fraction tendencies.
  num_modes = size(model%modes)
  do m=1,num_modes
    q_i = prognostics%interstitial_aerosols(m)
    q_c = prognostics%cloudborne_aerosols(m)
    dqdt_i = tendencies%interstitial_aerosols(m)
    dqdt_c = tendencies%cloudborne_aerosols(m)

    ! Cloudborne aerosols decay exponentially into interstitial aerosols.
    ! Modal number densities are unchanged.
    do k=1,model%num_levels
      do i=1,model%num_columns
        dqdt_i(k, i, m) = -tau * q_i(k, i, m)
        dqdt_c(k, i, m) = +tau * q_i(k, i, m)
      end do
    end do
  end do

  ! Gas mole fraction tendencies are zero.
  q_g = prognostics%gas_mole_fractions()
  dqdt_g = tendencies%gas_mole_fractions()
  dqdt_g(:,:,:) = 0.0_wp

  ! Modal number density tendencies are zero.
  n = prognostics%modal_num_densities()
  dndt = tendencies%modal_num_densities()
  dndt(:,:,:) = 0.0_wp

end subroutine

!> Disposes of the process-specific data allocated in prog_process_init.
!> You don't need this if you haven't allocated any process-specific data.
subroutine prog_stub_finalize() bind(c)
  implicit none

  ! Deallocate any process-specific resources
end subroutine

end module


