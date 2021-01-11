!> This module implements MAM4's nucleation process. For details, see the
!> appropriate section in the Processes chapter of the Haero design document.
module mam4_nucleation_mod

  use iso_c_binding, only: c_ptr
  use haero, only: wp, model, prognostics_t, diagnostics_t, tendencies_t, &
                   prognostics_from_c_ptr, diagnostics_from_c_ptr, &
                   tendencies_from_c_ptr

  implicit none
  private

  ! Module global variables
  ! All index variables are set to 0 if they don't correspond to anything within
  ! our aerosol model configuration.

  !> Aitken mode index
  integer :: aitken_index = 0

  !> Index of NH4 aerosol within Aitken mode
  integer :: nh4_aitken_index = 0

  !> Index of SO4 aerosol within Aitken mode
  integer :: so4_aitken_index = 0

  !> Index of H2SO4 gas
  integer :: h2so4_index = 0

  !> Index of NH3 gas
  integer :: nh3_index = 0

  public :: mam4_nucleation_init, &
            mam4_nucleation_run, &
            mam4_nucleation_finalize

contains

subroutine mam4_nucleation_init() bind(c)
  implicit none

  ! Record the aitken mode index.
  aitken_index = model%mode_index("aitken")

  ! Record the indices for aerosol species within the Aitken mode.
  nh4_aitken_index = model%aerosol_index(aitken_index, "SOA") ! NOTE: Correct?

  ! Record the indices for H2SO4 and NH3 gases.
  h2so4_index = model%gas_index("H2SO4")
  h2so4_index = model%gas_index("SOAG") ! NOTE: Is this correct?

end subroutine

subroutine mam4_nucleation_run(t, dt, progs, diags, tends) bind(c)
  use iso_c_binding, only: c_ptr, c_f_pointer
  implicit none

  ! Arguments
  real(wp), value, intent(in) :: t     ! simulation time
  real(wp), value, intent(in) :: dt    ! simulation time step
  type(c_ptr), value, intent(in) :: progs ! prognostic variables
  type(c_ptr), value, intent(in) :: diags ! diagnostic variables
  type(c_ptr), value, intent(in) :: tends ! tendencies

  ! Fortran prognostics, diagnostics, tendencies types
  type(prognostics_t) :: prognostics
  type(diagnostics_t) :: diagnostics
  type(tendencies_t)  :: tendencies

  ! Other local variables.
  integer :: m, i, k, s
  real(wp), pointer, dimension(:,:,:) :: q_c    ! cloudborne aerosol mix fracs
  real(wp), pointer, dimension(:,:,:) :: q_i    ! interstitial aerosol mix fracs
  real(wp), pointer, dimension(:,:,:) :: q_g    ! gas mole fracs
  real(wp), pointer, dimension(:,:,:) :: n      ! modal number densities
  real(wp), pointer, dimension(:,:,:) :: dqdt_c ! cloudborne aerosol tends
  real(wp), pointer, dimension(:,:,:) :: dqdt_i ! interstitial aerosol tends
  real(wp), pointer, dimension(:,:,:) :: dqdt_g ! gas mole frac tends
  real(wp), pointer, dimension(:,:,:) :: dndt   ! modal number density tends

  ! First of all, check to make sure our model has an aitken mode. If it
  ! doesn't, we can return immediately.
  if (aitken_index == 0) then
    return
  end if

  ! If there are no gases present with which to create new nuclei, there's
  ! nothing to do, either.
  if ((h2so4_index == 0) .and. (nh3_index == 0)) then
    return
  end if

  ! Finally, if there are no relevant aerosol species for nuclei, we can't
  ! create them.
  if ((nh4_aitken_index == 0) .and. (so4_aitken_index == 0)) then
    return
  end if

  ! Get Fortran data types from our C pointers.
  prognostics = prognostics_from_c_ptr(progs)
  diagnostics = diagnostics_from_c_ptr(diags)
  tendencies = tendencies_from_c_ptr(tends)

  ! Iterate over modes and compute aerosol mix fraction tendencies.
  do m=1,model%num_modes
    q_c => prognostics%cloudborne_aerosols(m)
    q_i => prognostics%interstitial_aerosols(m)
    dqdt_c => tendencies%cloudborne_aerosols(m)
    dqdt_i => tendencies%interstitial_aerosols(m)

  end do

  ! Gas mole fraction tendencies.
  q_g => prognostics%gas_mole_fractions()
  dqdt_g => tendencies%gas_mole_fractions()

  ! Modal number density tendencies.
  n => prognostics%modal_num_densities()
  dndt => tendencies%modal_num_densities()

end subroutine

subroutine mam4_nucleation_finalize() bind(c)
  implicit none

  ! Nothing to do here.
end subroutine

end module


