module mam_gasaerexch_test_bridge

  implicit none
  private

  ! Module functions
  public :: mam_soaexch_1subarea_bridge
  public :: mam_gasaerexch_1subarea_1gas_nonvolatile_bridge
  !public :: gas_aer_uptkrates_1box1gas_bridge

contains

subroutine mam_soaexch_1subarea_bridge( &
  lund, &
  dt, &
  temp, &
  pmid, &
  aircon, &
  n_mode, &
  ntot_amode, &
  max_mode, &
  qgas_cur, &
  qgas_avg, &
  qaer_cur, &
  qnum_cur, &
  uptkaer, &
  mode_aging_optaa, &
  lptr2_soa_a_amode) bind(c)

  use iso_c_binding, only: c_int
  use haero, only: modal_aero_config
  use haero_precision, only: wp
  use mam_gasaerexch, only: mam_soaexch_1subarea
  implicit none


! arguments
  integer(c_int), value, intent(in) :: lund                  ! logical unit for diagnostic output
  integer(c_int), value, intent(in) :: n_mode                ! current number of modes (including temporary)
  integer(c_int), value, intent(in) :: ntot_amode
  integer(c_int), value, intent(in) :: max_mode

  real(wp), value, intent(in) :: dt        ! current integration timestep (s)
  real(wp), value, intent(in) :: temp             ! temperature (K)
  real(wp), value, intent(in) :: pmid             ! pressure at model levels (Pa)
  real(wp), value, intent(in) :: aircon           ! air molar concentration (kmol/m3)

  real(wp), intent(inout), dimension( modal_aero_config%num_gases ) :: qgas_cur, qgas_avg
  real(wp), intent(inout), dimension( maxval(modal_aero_config%num_mode_species), modal_aero_config%num_aerosol_modes ) :: qaer_cur
  real(wp), intent(inout), dimension( modal_aero_config%num_aerosol_modes ) :: qnum_cur
  real(wp), intent(in   ), dimension( modal_aero_config%num_gases, modal_aero_config%num_aerosol_modes ) :: uptkaer
  integer(c_int),  intent(in   ), dimension( modal_aero_config%num_aerosol_modes) ::  mode_aging_optaa
  integer(c_int),  intent(in   ), dimension( modal_aero_config%num_aerosol_modes, *) ::  lptr2_soa_a_amode

  ! Call the actual subroutine.
  call  mam_soaexch_1subarea( &
  lund, &
  dt, &
  temp, &
  pmid, &
  aircon, &
  n_mode, &
  ntot_amode, &
  max_mode, &
  qgas_cur, &
  qgas_avg, &
  qaer_cur, &
  qnum_cur, &
  uptkaer, &
  mode_aging_optaa, &
  lptr2_soa_a_amode)
end subroutine

subroutine mam_gasaerexch_1subarea_1gas_nonvolatile_bridge( &
  dt, &
  qgas_netprod_otrproc, &
  n_mode, &
  uptkaer,  &
  qgas_cur,  &
  qgas_avg,  &
  qaer_cur) bind(c)

  use iso_c_binding, only: c_int
  use haero_precision, only: wp
  use mam_gasaerexch, only: mam_gasaerexch_1subarea_1gas_nonvolatile
  implicit none

  integer, parameter ::  max_mode = 5

  ! Arguments
  real(wp), intent(in) :: dt               ! current integration timestep (s)
  real(wp), intent(in) :: qgas_netprod_otrproc
                ! qgas_netprod_otrproc = gas net production rate from other processes
                !    such as gas-phase chemistry and emissions (mol/mol/s)
                ! this allows the condensation (gasaerexch) routine to apply production and condensation loss
                !    together, which is more accurate numerically
                ! NOTE - must be >= zero, as numerical method can fail when it is negative
                ! NOTE - currently only the values for h2so4 and nh3 should be non-zero

  integer(c_int),  intent(in) :: n_mode          ! current number of modes (including temporary)
  real(wp), intent(in)    :: uptkaer(1:max_mode) ! gas to aerosol mass transfer rate (1/s)
  real(wp), intent(inout) :: qgas_cur            ! current gas mix ratios (mol/mol)
  real(wp), intent(out)   :: qgas_avg            ! average gas mix ratios over the dt integration
  real(wp), intent(inout) :: qaer_cur(1:max_mode)! current aerosol mass mix ratios (mol/mol)


  ! Call the actual subroutine.
  ! But first set this public value on the module that the function will use.
  !call set_real_param("nuc_adjust_factor"   , factor_bin_tern_ratenucl)
  !call set_real_param("pbl_adjust_factor"   , factor_pbl_ratenucl)
  !call set_real_param("aitken_adjust_factor", aitken_adjust_factor)

  ! Call the actual subroutine.
  call mam_gasaerexch_1subarea_1gas_nonvolatile( &
  dt, &
  qgas_netprod_otrproc, &
  n_mode, &
  uptkaer,  &
  qgas_cur,  &
  qgas_avg,  &
  qaer_cur)
end subroutine

end module
