!> This module implements a stub for a diagnostic process written in Fortran.
!> It demonstrates how to use Haero's Fortran data types in a way that allows
!> such a process to communicate with the underlying C++ machinery.

!> Given an ambient atmospheric temperature and a set of numbers describing
!> (modal) aerosol and gas populations, this process computes their partial
!> pressures using the ideal gas law.
module diag_fprocess_stub

  use iso_c_binding, only: c_ptr
  use haero, only: wp, model, prognostics_t, atmosphere_t, diagnostics_t, &
                   prognostics_from_c_ptr, atmosphere_from_c_ptr, &
                   diagnostics_from_c_ptr, var_not_found

  implicit none
  private

  ! Process parameters
  real(wp), parameter :: R = 8.314472 ! universal gas constant [J/K/mol]
  real(wp), parameter :: Avogadro = 6.022e23 ! Avogadro's number [#/mol]
  real(wp), parameter :: dry_air_density = 1.292 ! Density of dry air [kg/m**3]
  real(wp), parameter :: dry_air_mw = 28.9647 ! Molecular weight of dry air [g/mol]
  real(wp), parameter :: kg_per_g = 0.001 ! conversion factor

  public :: diag_stub_init, &
            diag_stub_update, &
            diag_stub_finalize

contains

!> Performs initialization
subroutine diag_stub_init() bind(c)
  implicit none

  ! Nothing to do here.
end subroutine

!> Calls the update for the process.
subroutine diag_stub_update(t, progs, atm, diags) bind(c)
  use iso_c_binding, only: c_ptr, c_f_pointer
  implicit none

  ! Arguments
  real(wp), value, intent(in)    :: t     ! simulation time
  type(c_ptr), value, intent(in) :: progs ! prognostic variables
  type(c_ptr), value, intent(in) :: atm   ! atmosphere state variables
  type(c_ptr), value, intent(in) :: diags ! diagnostic variables

  ! Fortran prognostics and diagnostics types
  type(prognostics_t) :: prognostics
  type(atmosphere_t) :: atmosphere
  type(diagnostics_t) :: diagnostics

  ! Other local variables.
  integer :: num_modes, num_gases, m, k, s, p, token
  real(wp) :: mu_p, nu_p
  real(wp), pointer, dimension(:,:) :: q_c  ! cloudborne aerosol mix fracs
  real(wp), pointer, dimension(:,:) :: q_i  ! interstitial aerosol mix fracs
  real(wp), pointer, dimension(:,:) :: q_g  ! gas mole fracs
  real(wp), pointer, dimension(:,:) :: n    ! modal number densities
  real(wp), pointer, dimension(:)   :: temp ! atmospheric temperature
  real(wp), pointer, dimension(:,:) :: p_g  ! gas partial pressure
  real(wp), pointer, dimension(:,:) :: p_a  ! modal aerosol partial pressure

  ! Get Fortran data types from our C pointers.
  prognostics = prognostics_from_c_ptr(progs)
  atmosphere = atmosphere_from_c_ptr(atm)
  diagnostics = diagnostics_from_c_ptr(diags)
  n => prognostics%modal_num_concs()
  q_g => prognostics%gases()
  token = diagnostics%find_var("temperature")
  if (token /= var_not_found) temp => diagnostics%var(token)

  ! Zero the partial pressures.
  token = diagnostics%find_modal_var("pressure")
  if (token /= var_not_found) p_a => diagnostics%modal_var(token)
  p_a(:, :) = 0
  token = diagnostics%find_gas_var("pressure")
  if (token /= var_not_found) p_g => diagnostics%gas_var(token)
  p_g(:, :) = 0

  ! Diagnose modal aerosol partial pressure using ideal gas law
  num_modes = size(model%modes)
  q_c => prognostics%cloudborne_aerosols()
  q_i => prognostics%interstitial_aerosols()

  do p=1,model%num_populations
    do k=1,model%num_levels
      call model%get_mode_and_species(p, m, s)
      mu_p = model%aero_species(m, s)%molecular_wt ! species molecular weight
      ! Compute the number of aerosol moles per unit volume
      nu_p = (q_c(k, p) + q_i(k, p)) * n(k, m) / Avogadro
      p_a(k, m) = p_a(k, m) + nu_p * R * temp(k)
    end do
  end do

  ! Diagnose gas partial pressure using ideal gas law
  num_gases = size(model%gas_species)
  do s=1,num_gases
    do k=1,model%num_levels
      ! Compute the number of gas moles per unit volume
      nu_p = q_g(k, s) * dry_air_density / (kg_per_g * dry_air_mw)
      p_g(k, s) = nu_p * R * temp(k)
    end do
  end do

end subroutine

!> Disposes of the process-specific data allocated in diag_process_init.
!> You don't need this if you haven't allocated any process-specific data.
subroutine diag_stub_finalize() bind(c)
  implicit none

  ! Nothing to do here either.
end subroutine

end module

