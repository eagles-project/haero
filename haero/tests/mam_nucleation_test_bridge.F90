!> This module is a bridge between the MAMNucleationFProcess C++ class and the
!> mam_nucleation Fortran module.
module mam_nucleation_test_bridge

  implicit none
  private

  ! Module functions
  public :: ternary_nuc_merik2007_bridge
  public :: binary_nuc_vehk2002_bridge

contains

subroutine ternary_nuc_merik2007_bridge(t, rh, c2, c3, j_log, ntot, nacid, namm, r) bind(c)
  use haero_precision, only: wp
  use mam_nucleation, only: ternary_nuc_merik2007
  implicit none

  ! Arguments
  real(wp), value, intent(in) :: t, rh, c2, c3
  real(wp), intent(out) :: j_log, ntot, nacid, namm, r

  ! Call the actual subroutine.
  call ternary_nuc_merik2007(t, rh, c2, c3, j_log, ntot, nacid, namm, r)
end subroutine

subroutine binary_nuc_vehk2002_bridge(temp, rh, so4vol, ratenucl, rateloge, cnum_h2so4, cnum_tot, radius_cluster) bind(c)
  use haero_precision, only: wp
  use mam_nucleation, only: binary_nuc_vehk2002
  implicit none

  ! Arguments
  real(wp), value, intent(in) :: temp, rh, so4vol
  real(wp), intent(out) :: ratenucl, rateloge, cnum_h2so4, cnum_tot, radius_cluster

  ! Call the actual subroutine.
  call binary_nuc_vehk2002(temp, rh, so4vol, ratenucl, rateloge, cnum_h2so4, cnum_tot, radius_cluster)
end subroutine

subroutine pbl_nuc_wang2008_bridge(factor_pbl_ratenucl, so4vol, flagaa, flagaa2, ratenucl, rateloge, &
  cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster) bind(c)
  use iso_c_binding, only: c_int
  use haero_precision, only: wp
  use mam_nucleation, only: pbl_nuc_wang2008
  use mam_nucleation, only: adjust_factor_pbl_ratenucl
  implicit none

  ! Arguments
  real(wp), value, intent(in) :: factor_pbl_ratenucl
  real(wp), value, intent(in) :: so4vol
  integer(c_int), value, intent(in) :: flagaa
  integer(c_int), intent(inout) :: flagaa2
  real(wp), intent(inout) :: ratenucl, rateloge, cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster

  ! Call the actual subroutine.
  ! But first set this public value on the module that the function will use.
  adjust_factor_pbl_ratenucl = factor_pbl_ratenucl
  call pbl_nuc_wang2008(so4vol, flagaa, flagaa2, ratenucl, rateloge, &
    cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster)
end subroutine

subroutine mer07_veh02_nuc_mosaic_1box_bridge(   &
  factor_bin_tern_ratenucl, &
  factor_pbl_ratenucl, &
  newnuc_method_flagaa, &
  dtnuc, &
  temp_in, &
  rh_in, &
  press_in,   &
  zm_in, &
  pblh_in,   &
  qh2so4_cur, &
  qh2so4_avg, &
  qnh3_cur, &
  h2so4_uptkrate,   &
  mw_so4a_host,   &
  nsize, &
  maxd_asize, &
  dplom_sect, &
  dphim_sect,   &
  isize_nuc, &
  qnuma_del, &
  qso4a_del, &
  qnh4a_del,   &
  qh2so4_del, qnh3_del, dens_nh4so4a, ldiagaa,   &
  dnclusterdt ) bind(c)

  use iso_c_binding, only: c_int, c_ptr
  use haero_precision, only: wp
  use mam_nucleation, only: mer07_veh02_nuc_mosaic_1box
  use mam_nucleation, only: adjust_factor_bin_tern_ratenucl
  use mam_nucleation, only: adjust_factor_pbl_ratenucl
  use iso_c_binding, only: c_ptr, c_f_pointer
  implicit none

  ! Arguments

  real(wp), value, intent(in) :: factor_bin_tern_ratenucl
  real(wp), value, intent(in) :: factor_pbl_ratenucl
  real(wp), value, intent(in) :: dtnuc             ! nucleation time step (s)
  real(wp), value, intent(in) :: temp_in           ! temperature, in k
  real(wp), value, intent(in) :: rh_in             ! relative humidity, as fraction
  real(wp), value, intent(in) :: press_in          ! air pressure (pa)
  real(wp), value, intent(in) :: zm_in             ! layer midpoint height (m)
  real(wp), value, intent(in) :: pblh_in           ! pbl height (m)
  real(wp), value, intent(in) :: qh2so4_cur
  real(wp), value, intent(in) :: qh2so4_avg
  ! gas h2so4 mixing ratios (mol/mol-air)
  real(wp), value, intent(in) :: qnh3_cur          ! gas nh3 mixing ratios (mol/mol-air)
  ! qxxx_cur = current value (after gas chem and condensation)
  ! qxxx_avg = estimated average value (for simultaneous source/sink calcs)
  real(wp), value, intent(in) :: h2so4_uptkrate    ! h2so4 uptake rate to aerosol (1/s)
  real(wp), value, intent(in) :: mw_so4a_host      ! mw of so4 aerosol in host code (g/mol)

  integer(c_int), value, intent(in) :: newnuc_method_flagaa     ! 1=merikanto et al (2007) ternary
                                                         ! 2=vehkamaki et al (2002) binary
  integer(c_int), value, intent(in) :: nsize                    ! number of aerosol size bins
  integer(c_int), value, intent(in) :: maxd_asize               ! dimension for dplom_sect, ...
  type(c_ptr), value, intent(in) :: dplom_sect  ! dry diameter at lower bnd of bin (m)
  type(c_ptr), value, intent(in) :: dphim_sect  ! dry diameter at upper bnd of bin (m)
  integer(c_int), value, intent(in) :: ldiagaa

  ! arguments (out)
  integer(c_int), intent(out) :: isize_nuc         ! size bin into which new particles go
  real(wp), intent(out) :: qnuma_del        ! change to aerosol number mixing ratio (#/mol-air)
  real(wp), intent(out) :: qso4a_del        ! change to aerosol so4 mixing ratio (mol/mol-air)
  real(wp), intent(out) :: qnh4a_del        ! change to aerosol nh4 mixing ratio (mol/mol-air)
  real(wp), intent(out) :: qh2so4_del       ! change to gas h2so4 mixing ratio (mol/mol-air)
  real(wp), intent(out) :: qnh3_del         ! change to gas nh3 mixing ratio (mol/mol-air)
  ! aerosol changes are > 0; gas changes are < 0
  real(wp), intent(out) :: dens_nh4so4a     ! dry-density of the new nh4-so4 aerosol mass (kg/m3)
  real(wp), intent(out), optional :: dnclusterdt ! cluster nucleation rate (#/m3/s)

  real(wp), pointer, dimension(:)  :: dplom_sect_
  real(wp), pointer, dimension(:)  :: dphim_sect_

  call c_f_pointer(dplom_sect, dplom_sect_, shape=[maxd_asize])
  call c_f_pointer(dphim_sect, dphim_sect_, shape=[maxd_asize])

  ! Call the actual subroutine.
  ! But first set this public value on the module that the function will use.
  adjust_factor_bin_tern_ratenucl = factor_bin_tern_ratenucl
  adjust_factor_pbl_ratenucl      = factor_pbl_ratenucl
  call mer07_veh02_nuc_mosaic_1box(newnuc_method_flagaa, dtnuc, temp_in, rh_in, press_in,   &
    zm_in, pblh_in,   &
    qh2so4_cur, qh2so4_avg, qnh3_cur, h2so4_uptkrate,   &
    mw_so4a_host,   &
    nsize, maxd_asize, dplom_sect_, dphim_sect_,   &
    isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
    qh2so4_del, qnh3_del, dens_nh4so4a, ldiagaa,   &
    dnclusterdt )
end subroutine

end module

