!> This module implements a stub for a diagnostic process written in Fortran.
!> It demonstrates how to use Haero's Fortran data types in a way that allows
!> such a process to communicate with the underlying C++ machinery.

!> Given an ambient atmospheric temperature and a set of numbers describing
!> (modal) aerosol and gas populations, this process computes their partial
!> pressures using the ideal gas law.
module diag_process_stub

  use iso_c_binding, only: c_ptr
  use haero, only: wp, model, prognostics_t, diagnostics_t, &
                   prognostics_from_c_ptr, diagnostics_from_c_ptr

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
subroutine diag_stub_update(t, progs, diags) bind(c)
  use iso_c_binding, only: c_ptr, c_f_pointer
  implicit none

  ! Arguments
  real(wp), value, intent(in)    :: t     ! simulation time
  type(c_ptr), value, intent(in) :: progs ! prognostic variables
  type(c_ptr), value, intent(in) :: diags ! diagnostic variables

  ! Fortran prognostics and diagnostics types
  type(prognostics_t) :: prognostics
  type(diagnostics_t) :: diagnostics

  ! Other local variables.
  integer :: num_modes, num_gas_species, m, i, k, s
  real(wp) :: mu_s, nu_s
  real(wp), pointer, dimension(:,:,:) :: q_c  ! cloudborne aerosol mix fracs
  real(wp), pointer, dimension(:,:,:) :: q_i  ! interstitial aerosol mix fracs
  real(wp), pointer, dimension(:,:,:) :: q_g  ! gas mole fracs
  real(wp), pointer, dimension(:,:,:) :: n    ! modal number densities
  real(wp), pointer, dimension(:,:)   :: temp ! atmospheric temperature
  real(wp), pointer, dimension(:,:,:) :: p_g  ! gas partial pressure
  real(wp), pointer, dimension(:,:,:) :: p_a  ! modal aerosol partial pressure

  ! Get Fortran data types from our C pointers.
  prognostics = prognostics_from_c_ptr(progs)
  diagnostics = diagnostics_from_c_ptr(diags)
  n => prognostics%modal_num_densities()
  q_g => prognostics%gas_mole_fractions()
  temp => diagnostics%var("temperature")

  ! Zero the partial pressures.
  p_a => diagnostics%modal_var("pressure")
  p_a(:, :, :) = 0
  p_g => diagnostics%gas_var("pressure")
  p_g(:, :, :) = 0

  ! Diagnose modal aerosol partial pressure using ideal gas law
  num_modes = size(model%modes)
  do m=1,num_modes
    q_c => prognostics%cloudborne_aerosols(m)
    q_i => prognostics%interstitial_aerosols(m)

    do i=1,model%num_columns
      do k=1,model%num_levels
        do s=1,model%num_mode_species(m)
          mu_s = model%aero_species(m, s)%molecular_wt ! species molecular weight
          ! Compute the number of aerosol moles per unit volume
          nu_s = (q_c(s, k, i) + q_i(s, k, i)) * n(k, i, m) / Avogadro
          p_a(k, i, m) = p_a(k, i, m) + nu_s * R * temp(k, i)
        end do
      end do
    end do
  end do

  ! Diagnose gas partial pressure using ideal gas law
  num_gas_species = size(model%gas_species)
  do i=1,model%num_columns
    do k=1,model%num_levels
      do s=1,num_gas_species
        ! Compute the number of gas moles per unit volume
        nu_s = q_g(s, k, i) * dry_air_density / (kg_per_g * dry_air_mw)
        p_g(s, k, i) = nu_s * R * temp(k, i)
      end do
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

