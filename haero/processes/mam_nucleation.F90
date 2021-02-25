!> This module implements MAM's nucleation process. For details, see the
!> appropriate section in the Processes chapter of the Haero design document.
module mam_nucleation

  use haero, only: wp, model_t, species_t, &
                   prognostics_t, atmosphere_t, diagnostics_t, tendencies_t
  use haero_constants, only: pi, R_gas, Avogadro

  implicit none
  private

  ! Module functions
  public :: init, run, finalize

  !-------------------
  ! Module parameters
  !-------------------

  real(wp), parameter :: onethird = 1.0_wp/3.0_wp

  ! accomodation coefficient for h2so4 condensation
  real(wp), parameter :: accom_coef_h2so4 = 0.65_wp

  ! dry densities (kg/m3) molecular weights of aerosol
  ! ammsulf, ammbisulf, and sulfacid (from mosaic  dens_electrolyte values)
  !       real(wp), parameter :: dens_ammsulf   = 1.769e3
  !       real(wp), parameter :: dens_ammbisulf = 1.78e3
  !       real(wp), parameter :: dens_sulfacid  = 1.841e3
  ! use following to match cam3 modal_aero densities
  real(wp), parameter :: dens_ammsulf   = 1.770e3_wp
  real(wp), parameter :: dens_ammbisulf = 1.770e3_wp
  real(wp), parameter :: dens_sulfacid  = 1.770e3_wp

  ! molecular weights (g/mol) of aerosol ammsulf, ammbisulf, and sulfacid
  !    for ammbisulf and sulfacid, use 114 & 96 here rather than 115 & 98
  !    because we don't keep track of aerosol hion mass
  real(wp), parameter :: mw_ammsulf   = 132.0_wp
  real(wp), parameter :: mw_ammbisulf = 114.0_wp
  real(wp), parameter :: mw_sulfacid  =  96.0_wp

  ! Nucleation method parameter.
  !  1=merikanto et al (2007) ternary   2=vehkamaki et al (2002) binary
  ! 11=merikanto ternary + first-order boundary layer
  ! 12=merikanto ternary + second-order boundary layer
  integer, parameter :: newnuc_method_flagaa = 11
  ! integer, parameter :: newnuc_method_flagaa = 12

  ! min h2so4 vapor for nuc calcs = 4.0e-16 mol/mol-air ~= 1.0e4 molecules/cm3,
  real(wp), public, parameter :: qh2so4_cutoff = 4.0e-16_wp

  ! adjustment factors
!!!  real(wp), parameter, public :: adjust_factor_dnaitdt = 1.0_wp       ! applied to final dnait/dt. Redefined in modal_aero_amicphys

  real(wp), parameter, public :: adjust_factor_bin_tern_ratenucl = 1.0_wp  !  applied to binary/ternary nucleation rate
  ! real(wp), parameter, public :: adjust_factor_pbl_ratenucl = 1.0_wp  ! applied to boundary layer nucleation rate
  real(wp),            public :: adjust_factor_pbl_ratenucl = 1.0_wp  ! applied to boundary layer nucleation rate
                                                                      ! value reassigned in amicphys
  !-------------------------
  ! Module global variables
  !-------------------------

  ! All index variables are set to 0 if they don't correspond to anything within
  ! our aerosol model configuration.

  !> Aitken mode index
  integer :: aitken_index = 0

  !> Index of SO4 aerosol within Aitken mode
  integer :: so4_aitken_index = 0

  !> Index of H2SO4 gas
  integer :: h2so4_index = 0

  !> Molecular weight of SO4 aerosol
  real(wp) :: mw_so4

  integer, public :: newnuc_h2so4_conc_flag = 1

contains

subroutine init(model)
  implicit none

  ! Arguments
  type(model_t), intent(in) :: model

  type(species_t) so4
  integer :: so4_index

  ! Record the aitken mode index.
  aitken_index = model%mode_index("aitken")

  ! Record the index for SO4 aerosol.
  so4_index = model%aerosol_index(aitken_index, "SO4")
  so4_aitken_index = model%population_index(aitken_index, so4_index)

  ! Record the index for H2SO4 gas (source of new nuclei).
  h2so4_index = model%gas_index("H2SO4")

  ! Jot down some aerosol properties.
  so4 = model%aero_species(aitken_index, so4_aitken_index)
  mw_so4 = so4%molecular_wt

end subroutine

subroutine run(model, t, dt, prognostics, atmosphere, diagnostics, tendencies)
  implicit none

  ! Arguments
  type(model_t), intent(in)         :: model
  real(wp), value, intent(in)       :: t
  real(wp), value, intent(in)       :: dt
  type(prognostics_t), intent(in)   :: prognostics
  type(atmosphere_t), intent(in)    :: atmosphere
  type(diagnostics_t), intent(in)   :: diagnostics
  type(tendencies_t), intent(inout) :: tendencies

  ! Other local variables.
  integer :: i, k
  real(wp), pointer, dimension(:,:) :: q_c    ! cloudborne aerosol mix ratios
  real(wp), pointer, dimension(:,:) :: q_i    ! interstitial aerosol mix ratios
  real(wp), pointer, dimension(:,:) :: q_g    ! gas mix ratios
  real(wp), pointer, dimension(:,:) :: n      ! modal number concentrations
  real(wp), pointer, dimension(:,:) :: dqdt_c ! cloudborne aerosol tends
  real(wp), pointer, dimension(:,:) :: dqdt_i ! interstitial aerosol tends
  real(wp), pointer, dimension(:,:) :: dqdt_g ! gas mole frac tends
  real(wp), pointer, dimension(:,:) :: dndt   ! modal number density tends


  real(wp) :: q_so4     ! Mix ratio for sulfate aerosol [kg aerosol / kg dry air]
  real(wp) :: q_h2so4   ! Mix ratio for sulphuric acid gas [kg gas/kg dry air]
  real(wp) :: n_aitken  ! Number concentration for Aitken mode [#/kg dry air]
  real(wp) :: J_nuc     ! nucleation rate [#/cc/s]
  real(wp) :: dqndt_so4 ! tendency for SO4 number mixing ratio
  real(wp) :: md_so4    ! Dry mass for SO4 nucleus
  real(wp) :: temp      ! atmospheric temperature [K]
  real(wp) :: press     ! atmospheric pressure [Pa]
  real(wp) :: c_air     ! molar concentration of air [mol/m^3]
  real(wp) :: z_k       ! Elevation at center of kth vertical level
  real(wp) :: h_pbl     ! Planetary boundary layer height [m]
  real(wp) :: rel_hum   ! Relative humidity (0-1) [-]
  real(wp) :: uptkrate_h2so4, del_h2so4_gasprod, del_h2so4_aeruptk

  ! First of all, check to make sure our model has an aitken mode. If it
  ! doesn't, we can return immediately.
  if (aitken_index == 0) then
    return
  end if

  ! If there's no H2SO4 present with which to create new nuclei, there's
  ! nothing to do, either.
  if (h2so4_index == 0) then
    return
  end if

  ! Finally, if there are no relevant aerosol species for nuclei, we can't
  ! create them.
  if (so4_aitken_index == 0) then
    return
  end if

  ! Gas mole fraction tendencies.
  q_g => prognostics%gases()
  dqdt_g => tendencies%gases()

  ! Mix fractions and tendencies for SO4 aerosol in the Aitken mode.
  ! All new nuclei are deposited into interstitial aerosols.
  q_i => prognostics%interstitial_aerosols()
  dqdt_i => tendencies%interstitial_aerosols()

  ! Modal number density tendencies.
  n => prognostics%modal_num_concs()
  dndt => tendencies%modal_num_concs()

  ! Traverse the vertical levels and compute tendencies from nucleation.
  do k = 1,model%num_levels
    ! Compute the molar concentration of air at the given pressure and
    ! temperature.
    c_air = press/(temp*R_gas)

    ! Compute the nucleation rate of H2SO4.
    q_h2so4 = q_g(k, h2so4_index)
    q_so4 = q_i(k, so4_aitken_index)
    n_aitken = n(k, aitken_index)
    J_nuc = h2so4_nucleation_rate(q_so4, q_h2so4, n_aitken, dt, &
                                  temp, press, c_air, z_k, h_pbl, rel_hum, &
                                  uptkrate_h2so4, del_h2so4_gasprod, &
                                  del_h2so4_aeruptk)

    ! Compute the dry mass of an SO4 nucleus.
    md_so4 = 1_wp ! FIXME

    ! Compute tendencies given J_nuc.
    dqndt_so4 = 1e6_wp * J_nuc / c_air
    dqdt_i(k, so4_aitken_index) = dqndt_so4 * (md_so4 / mw_so4)
    dqdt_g(k, h2so4_index) = -dqdt_i(k, so4_aitken_index)
    dndt(k, aitken_index) = dqndt_so4 * c_air
  end do
end subroutine

! Calculates new particle production from homogeneous nucleation over timestep
! dt, using nucleation rates from Vehkamaki et al. (2002) h2so4-h2o binary
! parameterization.
!
! The new particles are "grown" to the lower-bound size of the host code's
! smallest size bin. (this "growth" is somewhat ad hoc, and would not be
! necessary if the host code's size bins extended down to ~1 nm.)
!
! If the h2so4 mass mixing ratio of the grown new particles exceed the current
! gas mixing ratios, the new particle production is reduced so that the new
! particle mass ratios match the gas ratios.
!
! The correction of kerminen and kulmala (2002) is applied to account
! for loss of the new particles by coagulation as they are
! growing to the "host code mininum size".
!
! References:
!
!  * Vehkamäki, H., M. Kulmala, I. Napari, K.E.J. Lehtinen,
!       C. Timmreck, M. Noppel and A. Laaksonen, 2002,
!       An improved parameterization for sulfuric acid-water nucleation
!       rates for tropospheric and stratospheric conditions,
!       J. Geophys. Res., 107, 4622, doi:10.1029/2002jd002184
!
!  * Kerminen, V., and M. Kulmala, 2002,
!   Analytical formulae connecting the "real" and the "apparent"
!   nucleation rate and the nuclei number concentration
!   for atmospheric nucleation events
function h2so4_nucleation_rate(q_so4, q_h2so4, n_aitken, dt, temp, pmid, aircon, &
                               zmid, pblh, relhum, uptkrate_h2so4, &
                               del_h2so4_gasprod, del_h2so4_aeruptk) result(J_nuc)

  implicit none

  real(wp), intent(in) :: q_so4             ! SO4 aerosol mixing ratio [kg/kg]
  real(wp), intent(in) :: q_h2so4           ! H2SO4 gas mole fraction [mol/mol]
  real(wp), intent(in) :: n_aitken          ! Aitken mode number density [#/m^3]
  real(wp), intent(in) :: dt                ! simulation time step [s]
  real(wp), intent(in) :: temp              ! temperature [K]
  real(wp), intent(in) :: pmid              ! pressure at vertical levels [Pa]
  real(wp), intent(in) :: aircon            ! air molar concentration (kmol/m3)
  real(wp), intent(in) :: zmid              ! midpoint height above surface [m]
  real(wp), intent(in) :: pblh              ! planetary boundary layer height [m]
  real(wp), intent(in) :: relhum            ! relative humidity (0-1) [-]
  real(wp), intent(in) :: uptkrate_h2so4    ! H2SO4 gas uptake rate
  real(wp), intent(in) :: del_h2so4_gasprod
  real(wp), intent(in) :: del_h2so4_aeruptk

  real(wp) :: J_nuc

  integer :: itmp
  integer :: l
  integer :: m

  real(wp) :: dmdt_ait, dmdt_aitsv1, dmdt_aitsv2, dmdt_aitsv3
  real(wp) :: dndt_ait, dndt_aitsv1, dndt_aitsv2, dndt_aitsv3
  real(wp) :: dso4dt_ait
  real(wp) :: dpnuc
  real(wp) :: dplom_mode(1), dphim_mode(1)
  real(wp) :: mass1p
  real(wp) :: mass1p_aithi, mass1p_aitlo
  real(wp) :: qh2so4_cur, qh2so4_avg, qh2so4_del
  real(wp) :: qnuma_del
  real(wp) :: qso4a_del
  real(wp) :: relhumnn
  real(wp) :: tmp_q2, tmp_q3
  real(wp) :: tmp_q_del
  real(wp) :: tmp_frso4, tmp_uptkrate
  real(wp) :: temp_bb
  real(wp) :: rh_bb
  real(wp) :: c_h2so4_bb

  ! min h2so4 vapor for nuc calcs = 4.0e-16 mol/mol-air ~= 1.0e4 molecules/cm3,
  real(wp), parameter :: q_h2so4_cutoff = 4.0e-16_wp

#if 0
  J_nuc = 0_wp

  ! skip if h2so4 vapor < qh2so4_cutoff
  if (q_h2so4 <= q_h2so4_cutoff) then
    return
  end if

  ! estimate the average mixing ratio for H2SO4 using the method in standard
  ! cam5.2 modal_aero_newnuc
  tmpa = max( 0.0_wp, del_h2so4_gasprod )
  tmp_q3 = qh2so4_cur
  ! tmp_q2 = qh2so4 before aeruptk
  ! (note tmp_q3, tmp_q2 both >= 0.0)
  tmp_q2 = tmp_q3 + max( 0.0_wp, -del_h2so4_aeruptk )

  ! tmpb = log( tmp_q2/tmp_q3 ) BUT with some checks added
  if (tmp_q2 <= tmp_q3) then
    tmpb = 0.0_wp
  else
    tmpc = tmp_q2 * exp( -20.0_wp )
    if (tmp_q3 <= tmpc) then
      tmp_q3 = tmpc
      tmpb = 20.0_wp
    else
      tmpb = log( tmp_q2/tmp_q3 )
    end if
  end if
  ! d[ln(qh2so4)]/dt [1/s] from uptake (condensation) to aerosol
  tmp_uptkrate = tmpb/dt

  ! qh2so4_avg = estimated average qh2so4
  ! when production & loss are done simultaneously
  if (tmpb <= 0.1_wp) then
    qh2so4_avg = tmp_q3*(1.0_wp + 0.5_wp*tmpb) - 0.5_wp*tmpa
  else
    tmpc = tmpa/tmpb
    qh2so4_avg = (tmp_q3 - tmpc)*((exp(tmpb)-1.0_wp)/tmpb) + tmpc
  end if

  ! If the average is too low, there is no nucleation.
  if (qh2so4_avg <= qh2so4_cutoff) then
    return
  end if

  ! dry-diameter limits for "grown" new particles
  dplom_mode(1) = exp(0.67_wp*log(dgnumlo_aer(nait))  &
                      + 0.33_wp*log(dgnum_aer(nait)))
  dphim_mode(1) = dgnumhi_aer(nait)

  ! mass1p_... = mass (kg) of so4 in a single particle of diameter ...
  ! (assuming same dry density for so4)
  ! mass1p_aitlo - dp = dplom_mode(1)
  ! mass1p_aithi - dp = dphim_mode(1)
  tmpa = dens_so4a_host*pi/6.0_wp
  mass1p_aitlo = tmpa*(dplom_mode(1)**3)
  mass1p_aithi = tmpa*(dphim_mode(1)**3)

  ! limit RH to between 0.1% and 99%
  relhumnn = max( 0.01_wp, min( 0.99_wp, relhum ) )

  ! Calc h2so4 in molecules/cm3
  cair = press_in/(temp_in*R_gas)
  c_h2so4_in = qh2so4_avg * cair * Avogadro * 1.0e-6_wp

  ! Compute the intermediate nucleation rate (Vehkamaki 2002).
  J_star = 1.0e-38_wp
  log_J_star = log(J_star)
  if (so4vol_in >= 1.0e4_wp) then
    temp_bb = max( 230.15_wp, min( 305.15_wp, temp_in ) )
    rh_bb = max( 1.0e-4_wp, min( 1.0_wp, rh_in ) )
    c_h2so4_bb = max( 1.0e4_wp, min( 1.0e11_wp, c_h2so4_in ) )
    call binary_nuc_vehk2002(temp_bb, rh_bb, c_h2so4_bb, &
      J_star, N_h2so4, N_tot, r_star)
  end if

  ! Update log(J_star)
  log_J_star = log(J_star) + log( max( 1.0e-38_wp, adjust_factor_bin_ratenucl ) )

  ! Apply the planetary boundary layer correction where needed
  if ( zm_in <= max(pblh_in, 100.0_wp) ) then
    c_h2so4_bb = c_h2so4_in
    call pbl_nuc_wang2008( c_h2so4_bb, J_star, log_J_star, N_tot, &
                           N_h2so4, r_star )
    log_J_star = log(J_star)
  end if

  ! if nucleation rate is less than 1e-6 #/cm3/s ~= 0.1 #/cm3/day,
  ! exit with new particle formation = 0
  if (log_J_star <= -13.82_wp) then
    return
  end if
  J_star = exp(log_J_star)

  ! Now grow the nuclei and adjust the nucleation rate accordingly.
  J_nuc = growth_adjusted_nucleation_rate(J_star)

  ! max production of aerosol dry mass (kg-aero/m3-air)
  tmpa = max( 0.0_wp, (ratenuclt_kk*dt*mass_part) )
  ! max production of aerosol so4 (mol-so4a/mol-air)
  tmpe = tmpa/(kgaero_per_moleso4a*cair)
  ! max production of aerosol so4 (mol/mol-air)
  ! based on ratenuclt_kk and mass_part
  qmolso4a_del_max = tmpe

  ! check if max production exceeds available h2so4 vapor
  freduce = 1.0_wp
  if (qmolso4a_del_max > qh2so4_cur) then
    freduce = qh2so4_cur/qmolso4a_del_max
  end if

  ! if adjusted nucleation rate is less than 1e-12 #/m3/s ~= 0.1 #/cm3/day,
  ! exit with new particle formation = 0
  if (J_nuc*freduce <= 1.0e-12_wp) then
    J_nuc = 0_wp
    return
  end if

  ! convert qnuma_del from (#/mol-air) to (#/kmol-air)
  qnuma_del = qnuma_del*1.0e3_wp

  ! number nuc rate (#/kmol-air/s) from number nuc amt
  dndt_aitken = qnuma_del/dt

  ! fraction of mass nuc going to so4
  tmpa = qso4a_del*mw_so4a_host
  tmpb = tmpa
  tmp_frso4 = 1.0_wp

  ! mass nuc rate (kg/kmol-air/s) from mass nuc amts
  dmdt_ait = max( 0.0_wp, (tmpb/dt) )

  dndt_aitsv1 = dndt_aitken
  dmdt_aitsv1 = dmdt_ait
  dndt_aitsv2 = 0.0_wp
  dmdt_aitsv2 = 0.0_wp
  dndt_aitsv3 = 0.0_wp
  dmdt_aitsv3 = 0.0_wp

  if (dndt_aitken < 1.0e2_wp) then
    ! ignore newnuc if number rate < 100 #/kmol-air/s ~= 0.3 #/mg-air/d
    dndt_aitken = 0.0_wp
    dmdt_ait = 0.0_wp
  else
    dndt_aitsv2 = dndt_ait
    dmdt_aitsv2 = dmdt_ait

    ! mirage2 code checked for complete h2so4 depletion here,
    ! but this is now done in veh02_nuc_mosaic_1box
    mass1p = dmdt_ait/dndt_ait
    dndt_aitsv3 = dndt_ait
    dmdt_aitsv3 = dmdt_ait

    ! apply particle size constraints
    if (mass1p < mass1p_aitlo) then
      ! reduce dndt to increase new particle size
      dndt_aitken = dmdt_ait/mass1p_aitlo
    else if (mass1p > mass1p_aithi) then
      ! reduce dmdt to decrease new particle size
      dmdt_ait = dndt_ait*mass1p_aithi
    end if
  end if

  ! *** apply adjustment factor to avoid unrealistically high
  ! aitken number concentrations in mid and upper troposphere
  dndt_aitken = dndt_aitken * newnuc_adjust_factor_dnaitdt
  dmdt_ait = dmdt_ait * newnuc_adjust_factor_dnaitdt

  ! Compute the tendency for new SO4 aerosol
  dqdt_so4 = dmdt_ait*tmp_frso4/mw_so4a_host

  ! Compute the tendency for expended H2SO4 gas
  dqdt_h2so4 = dmdt_ait*tmp_frso4/mw_so4a_host
#endif
end function

subroutine finalize(model)
  implicit none

  ! Arguments
  type(model_t), intent(in) :: model

  ! Nothing to do here.
end subroutine

subroutine veh02_nuc_mosaic_1box(dt, temp_in, rh_in, press_in,   &
  zm_in, pblh_in, qh2so4_cur, qh2so4_avg, h2so4_uptkrate,   &
  mw_so4a_host, nsize, maxd_asize, dplom_sect, dphim_sect,   &
  isize_nuc, qnuma_del, qso4a_del, qh2so4_del)

  use haero_constants, only: &
    R_gas, &
    mw_so4a => molec_weight_so4

  implicit none

  ! subr arguments (in)
  real(wp), intent(in) :: dt                ! nucleation time step (s)
  real(wp), intent(in) :: temp_in           ! temperature, in k
  real(wp), intent(in) :: rh_in             ! relative humidity, as fraction
  real(wp), intent(in) :: press_in          ! air pressure (pa)
  real(wp), intent(in) :: zm_in             ! layer midpoint height (m)
  real(wp), intent(in) :: pblh_in           ! pbl height (m)

  ! gas h2so4 mixing ratios (mol/mol-air)
  real(wp), intent(in) :: qh2so4_cur, qh2so4_avg
  ! qxxx_cur = current value (after gas chem and condensation)
  ! qxxx_avg = estimated average value (for simultaneous source/sink calcs)
  real(wp), intent(in) :: h2so4_uptkrate    ! h2so4 uptake rate to aerosol (1/s)
  real(wp), intent(in) :: mw_so4a_host      ! mw of so4 aerosol in host code (g/mol)

  integer, intent(in) :: nsize                    ! number of aerosol size bins
  integer, intent(in) :: maxd_asize               ! dimension for dplom_sect, ...
  real(wp), intent(in) :: dplom_sect(maxd_asize)  ! dry diameter at lower bnd of bin (m)
  real(wp), intent(in) :: dphim_sect(maxd_asize)  ! dry diameter at upper bnd of bin (m)

  ! subr arguments (out)
  integer, intent(out) :: isize_nuc   ! size bin into which new particles go
  real(wp), intent(out) :: qnuma_del  ! change to aerosol number mixing ratio (#/mol-air)
  real(wp), intent(out) :: qso4a_del  ! change to aerosol so4 mixing ratio (mol/mol-air)
  real(wp), intent(out) :: qh2so4_del ! change to gas h2so4 mixing ratio (mol/mol-air)
  ! aerosol changes are > 0; gas changes are < 0

  real(wp) :: J_star     ! binary nucleation rate [#/cc/s]
  real(wp) :: log_J_star ! ln (J_star)
  real(wp) :: N_h2so4    ! # of H2SO4 molecules in the critical cluster [#]
  real(wp) :: N_tot      ! total # of molecules in the critical cluster [#]
  real(wp) :: r_star     ! radius of critical cluster [nm]

  ! local variables
  integer :: i
  integer :: igrow

  real(wp) cair                     ! dry-air molar density (mol/m3)
  real(wp) cs_prime_kk              ! kk2002 "cs_prime" parameter (1/m2)
  real(wp) cs_kk                    ! kk2002 "cs" parameter (1/s)
  real(wp) dens_part                ! "grown" single-particle dry density (kg/m3)
  real(wp) dfin_kk, dnuc_kk         ! kk2002 final/initial new particle wet diameter (nm)
  real(wp) tmpa, tmpb, tmpc, tmpe, tmpq
  real(wp) tmpa1, tmpb1
  real(wp) tmp_m1, tmp_m2, tmp_m3, tmp_n1, tmp_n2, tmp_n3
  real(wp) tmp_spd                  ! h2so4 vapor molecular speed (m/s)
  real(wp) factor_kk
  real(wp) fogas, foso4a, fonuma
  real(wp) freduce                  ! reduction factor applied to nucleation
                                    ! rate due to limited availability of
                                    ! h2so4 gas
  real(wp) gr_kk                    ! kk2002 "gr" parameter (nm/h)
  real(wp) kgaero_per_moleso4a      ! (kg dry aerosol)/(mol aerosol so4)
  real(wp) mass_part                ! "grown" single-particle dry mass (kg)
  real(wp) nu_kk                    ! kk2002 "nu" parameter (nm)
  real(wp) qmolso4a_del_max         ! max production of aerosol so4 over dt (mol/mol-air)
  real(wp) ratenuclt_bb             ! nucleation rate (#/m3/s)
  real(wp) ratenuclt_kk             ! nucleation rate after kk2002 adjustment (#/m3/s)
  real(wp) rh_bb                    ! bounded value of rh_in
  real(wp) so4vol_in                ! concentration of h2so4 for nucl. calc., molecules cm-3
  real(wp) so4vol_bb                ! bounded value of so4vol_in
  real(wp) temp_bb                  ! bounded value of temp_in
  real(wp) voldry_part              ! "grown" single-particle dry volume (m3)
  real(wp) wet_volfrac_so4a         ! grown particle (dry-volume-from-so4)/(wet-volume)
end subroutine

! calculates binary nucleation rate and critical cluster size
! using the parameterization in
!     vehkamäki, h., m. kulmala, i. napari, k.e.j. lehtinen,
!        c. timmreck, m. noppel and a. laaksonen, 2002,
!        an improved parameterization for sulfuric acid-water nucleation
!        rates for tropospheric and stratospheric conditions,
!        j. geophys. res., 107, 4622, doi:10.1029/2002jd002184
subroutine binary_nuc_vehk2002(temp, rh, c_h2so4, &
                               J_star, N_h2so4, N_tot, r_star)

  implicit none

  ! subr arguments (in)
  real(wp), intent(in) :: temp      ! temperature [K]
  real(wp), intent(in) :: rh        ! relative humidity (0-1) [-]
  real(wp), intent(in) :: c_h2so4   ! concentration of h2so4 [#/cc]

  ! subr arguments (out)
  real(wp), intent(out) :: J_star   ! binary nucleation rate [#/cc/s]
  real(wp), intent(out) :: N_h2so4  ! # h2so4 molecules in critical cluster [#]
  real(wp), intent(out) :: N_tot    ! total # molecules in critical cluster [#]
  real(wp), intent(out) :: r_star   ! radius of critical cluster [nm]

  ! local variables
  real(wp) :: crit_x
  real(wp) :: acoe, bcoe, ccoe, dcoe, ecoe, fcoe, gcoe, hcoe, icoe, jcoe
  real(wp) :: log_J_star

  ! calc sulfuric acid mole fraction in critical cluster
  crit_x = 0.740997_wp - 0.00266379_wp * temp   &
    - 0.00349998_wp * log (c_h2so4)   &
    + 0.0000504022_wp * temp * log (c_h2so4)   &
    + 0.00201048_wp * log (rh)   &
    - 0.000183289_wp * temp * log (rh)   &
    + 0.00157407_wp * (log (rh)) ** 2.0_wp   &
    - 0.0000179059_wp * temp * (log (rh)) ** 2.0_wp   &
    + 0.000184403_wp * (log (rh)) ** 3.0_wp   &
    - 1.50345e-6_wp * temp * (log (rh)) ** 3.0_wp

  ! calc nucleation rate
  acoe    = 0.14309_wp+2.21956_wp*temp   &
    - 0.0273911_wp * temp**2.0_wp   &
    + 0.0000722811_wp * temp**3.0_wp + 5.91822_wp/crit_x

  bcoe    = 0.117489_wp + 0.462532_wp *temp   &
    - 0.0118059_wp * temp**2.0_wp   &
    + 0.0000404196_wp * temp**3.0_wp + 15.7963_wp/crit_x

  ccoe    = -0.215554_wp-0.0810269_wp * temp   &
    + 0.00143581_wp * temp**2.0_wp   &
    - 4.7758e-6_wp * temp**3.0_wp   &
    - 2.91297_wp/crit_x

  dcoe    = -3.58856_wp+0.049508_wp * temp   &
    - 0.00021382_wp * temp**2.0_wp   &
    + 3.10801e-7_wp * temp**3.0_wp   &
    - 0.0293333_wp/crit_x

  ecoe    = 1.14598_wp - 0.600796_wp * temp   &
    + 0.00864245_wp * temp**2.0_wp   &
    - 0.0000228947_wp * temp**3.0_wp   &
    - 8.44985_wp/crit_x

  fcoe    = 2.15855_wp + 0.0808121_wp * temp   &
    -0.000407382_wp * temp**2.0_wp   &
    -4.01957e-7_wp * temp**3.0_wp   &
    + 0.721326_wp/crit_x

  gcoe    = 1.6241_wp - 0.0160106_wp * temp   &
    + 0.0000377124_wp * temp**2.0_wp   &
    + 3.21794e-8_wp * temp**3.0_wp   &
    - 0.0113255_wp/crit_x

  hcoe    = 9.71682_wp - 0.115048_wp * temp   &
    + 0.000157098_wp * temp**2.0_wp   &
    + 4.00914e-7_wp * temp**3.0_wp   &
    + 0.71186_wp/crit_x

  icoe    = -1.05611_wp + 0.00903378_wp * temp   &
    - 0.0000198417_wp * temp**2.0_wp   &
    + 2.46048e-8_wp  * temp**3.0_wp   &
    - 0.0579087_wp/crit_x

  jcoe    = -0.148712_wp + 0.00283508_wp * temp   &
    - 9.24619e-6_wp  * temp**2.0_wp   &
    + 5.00427e-9_wp * temp**3.0_wp   &
    - 0.0127081_wp/crit_x

  log_J_star = (   &
    acoe   &
    + bcoe * log (rh)   &
    + ccoe * ( log (rh))**2.0_wp   &
    + dcoe * ( log (rh))**3.0_wp   &
    + ecoe * log (c_h2so4)   &
    + fcoe * (log (rh)) * (log (c_h2so4))   &
    + gcoe * ((log (rh) ) **2.0_wp)   &
    * (log (c_h2so4))   &
    + hcoe * (log (c_h2so4)) **2.0_wp   &
    + icoe * log (rh)   &
    * ((log (c_h2so4)) **2.0_wp)   &
    + jcoe * (log (c_h2so4)) **3.0_wp   &
    )
  log_J_star = min( log_J_star, log(1.0e38_wp) )
  J_star = exp ( log_J_star )

  ! calc number of molecules in critical cluster
  acoe    = -0.00295413_wp - 0.0976834_wp*temp   &
    + 0.00102485_wp * temp**2.0_wp   &
    - 2.18646e-6_wp * temp**3.0_wp - 0.101717_wp/crit_x

  bcoe    = -0.00205064_wp - 0.00758504_wp*temp   &
    + 0.000192654_wp * temp**2.0_wp   &
    - 6.7043e-7_wp * temp**3.0_wp - 0.255774_wp/crit_x

  ccoe    = +0.00322308_wp + 0.000852637_wp * temp   &
    - 0.0000154757_wp * temp**2.0_wp   &
    + 5.66661e-8_wp * temp**3.0_wp   &
    + 0.0338444_wp/crit_x

  dcoe    = +0.0474323_wp - 0.000625104_wp * temp   &
    + 2.65066e-6_wp * temp**2.0_wp   &
    - 3.67471e-9_wp * temp**3.0_wp   &
    - 0.000267251_wp/crit_x

  ecoe    = -0.0125211_wp + 0.00580655_wp * temp   &
    - 0.000101674_wp * temp**2.0_wp   &
    + 2.88195e-7_wp * temp**3.0_wp   &
    + 0.0942243_wp/crit_x

  fcoe    = -0.038546_wp - 0.000672316_wp * temp   &
    + 2.60288e-6_wp * temp**2.0_wp   &
    + 1.19416e-8_wp * temp**3.0_wp   &
    - 0.00851515_wp/crit_x

  gcoe    = -0.0183749_wp + 0.000172072_wp * temp   &
    - 3.71766e-7_wp * temp**2.0_wp   &
    - 5.14875e-10_wp * temp**3.0_wp   &
    + 0.00026866_wp/crit_x

  hcoe    = -0.0619974_wp + 0.000906958_wp * temp   &
    - 9.11728e-7_wp * temp**2.0_wp   &
    - 5.36796e-9_wp * temp**3.0_wp   &
    - 0.00774234_wp/crit_x

  icoe    = +0.0121827_wp - 0.00010665_wp * temp   &
    + 2.5346e-7_wp * temp**2.0_wp   &
    - 3.63519e-10_wp * temp**3.0_wp   &
    + 0.000610065_wp/crit_x

  jcoe    = +0.000320184_wp - 0.0000174762_wp * temp   &
    + 6.06504e-8_wp * temp**2.0_wp   &
    - 1.4177e-11_wp * temp**3.0_wp   &
    + 0.000135751_wp/crit_x

  N_tot = exp (   &
    acoe   &
    + bcoe * log (rh)   &
    + ccoe * ( log (rh))**2.0_wp   &
    + dcoe * ( log (rh))**3.0_wp   &
    + ecoe * log (c_h2so4)   &
    + fcoe * (log (rh)) * (log (c_h2so4))   &
    + gcoe * ((log (rh) ) **2.0_wp)   &
    * (log (c_h2so4))   &
    + hcoe * (log (c_h2so4)) **2.0_wp   &
    + icoe * log (rh)   &
    * ((log (c_h2so4)) **2.0_wp)   &
    + jcoe * (log (c_h2so4)) **3.0_wp   &
    )

  N_h2so4 = N_tot * crit_x

  !   calc radius (nm) of critical cluster
  r_star = exp( -1.6524245_wp + 0.42316402_wp*crit_x   &
    + 0.3346648_wp*log(N_tot) )

end subroutine

! calculates boundary nucleation nucleation rate
! using the first or second-order parameterization in
!     wang, m., and j.e. penner, 2008,
!        aerosol indirect forcing in a global model with particle nucleation,
!        atmos. chem. phys. discuss., 8, 13943-13998
subroutine pbl_nuc_wang2008(c_h2so4, J_star, N_tot, N_h2so4, r_star)

  implicit none

  real(wp), intent(in)    :: c_h2so4 ! concentration of h2so4 [#/cc]
  real(wp), intent(inout) :: J_star  ! binary nucleation rate [#/cc/s]
  real(wp), intent(inout) :: N_tot   ! total # of molecules in critical cluster [#]
  real(wp), intent(inout) :: N_h2so4 ! # of h2so4 molecules in critical cluster [#]
  real(wp), intent(inout) :: r_star  ! radius of critical cluster [nm]

  ! local variables
  real(wp) :: tmp_diam, tmp_mass, tmp_volu
  real(wp) :: J_pbl, log_J_pbl

  ! empirical estimation of nucleation rate
  J_pbl = 1.0e-6_wp * c_h2so4
  J_pbl = J_pbl * adjust_factor_pbl_ratenucl
  log_J_pbl = log(max( 1.0e-38_wp, J_pbl))

  ! exit if pbl nuc rate is lower than (incoming) binary rate
  if (log_J_pbl <= log(J_star)) then
    return
  end if

  J_star = J_pbl

  ! following wang 2002, assume fresh nuclei are 1 nm diameter
  !    subsequent code will "grow" them to aitken mode size
  r_star = 0.5_wp

  ! assume fresh nuclei are pure h2so4
  !    since aitken size >> initial size, the initial composition
  !    has very little impact on the results
  tmp_diam = r_star * 2.0e-7_wp   ! diameter in cm
  tmp_volu = (tmp_diam**3) * (pi/6.0_wp)  ! volume in cm^3
  tmp_mass = tmp_volu * 1.8_wp            ! mass in g
  N_h2so4 = (tmp_mass / 98.0_wp) * 6.023e23_wp   ! no. of h2so4 molec assuming pure h2so4
  N_tot = N_h2so4
end subroutine

! Computes the final nucleation rate J_nuc from the intermediate nucleation rate
! J_star after having grown the SO4 nuclei so that they can fit into the
! Aitken mode.
function growth_adjusted_nucleation_rate(J_star) result (J_nuc)
  real(wp), intent(in) :: J_star
  real(wp)             :: J_nuc

#if 0
  real(wp) :: J_star_bb
  real(wp) :: tmpa, tmpb, tmpc
  real(wp) :: wetvol_dryvol    ! grown particle (wet-volume)/(dry-volume)
  real(wp) :: voldry_clus      ! critical-cluster dry volume (m3)
  real(wp) :: dpdry_clus       ! critical cluster diameter (m)
  real(wp) :: dpdry_part       ! "grown" single-particle dry diameter (m)
  real(wp) :: gamma_kk         ! kk2002 "gamma" parameter (nm2*m2/h)

  J_star_bb = J_star*1.0e6_wp  ! J_star_bb is #/m3/s; J_star is #/cm3/s

  ! wet/dry volume ratio - use simple kohler approx for amm_sulf/amm_bisulf
  tmpa = max( 0.10_wp, min( 0.95_wp, rh_in ) )
  wetvol_dryvol = 1.0_wp - 0.56_wp/log(tmpa)

  ! determine size bin into which the new particles go
  ! (probably it will always be bin #1, but ...)
  voldry_clus = ( max(N_h2so4,1.0_wp)*mw_so4a) /   &
    (1.0e3_wp*dens_sulf_acid*Avogadro)
  ! correction when host code sulfate is really ammonium bisulfate/sulfate
  voldry_clus = voldry_clus * (mw_so4a_host/mw_so4a)
  dpdry_clus = (voldry_clus*6.0_wp/pi)**one_third

  isize_nuc = 1
  dpdry_part = dplom_sect(1)
  if (dpdry_clus <= dplom_sect(1)) then
    igrow = 1   ! need to clusters to larger size
  else if (dpdry_clus >= dphim_sect(nsize)) then
    igrow = 0
    isize_nuc = nsize
    dpdry_part = dphim_sect(nsize)
  else
    igrow = 0
    do i = 1, nsize
      if (dpdry_clus < dphim_sect(i)) then
        isize_nuc = i
        dpdry_part = dpdry_clus
        dpdry_part = min( dpdry_part, dphim_sect(i) )
        dpdry_part = max( dpdry_part, dplom_sect(i) )
        exit
      end if
    end do
  end if
  voldry_part = (pi/6.0_wp)*(dpdry_part**3)

  ! All "grown particles" are SO4.
  if (igrow <= 0) then
    ! no "growing" so pure sulfuric acid
    tmp_n1 = 0.0_wp
    tmp_n2 = 0.0_wp
    tmp_n3 = 1.0_wp
  else
    ! combination of ammonium bisulfate and sulfuric acid
    ! tmp_n2 & tmp_n3 = mole fractions of the amm_bisulf & sulf_acid
    tmp_n1 = 0.0_wp
    tmp_n2 = 0.0_wp ! (qnh3_cur/qh2so4_cur)
    tmp_n2 = max( 0.0_wp, min( 1.0_wp, tmp_n2 ) )
    tmp_n3 = 1.0_wp - tmp_n2
  end if

  tmp_m1 = tmp_n1*mw_amm_sulf
  tmp_m2 = tmp_n2*mw_amm_bisulf
  tmp_m3 = tmp_n3*mw_sulf_acid
  dens_part = (tmp_m1 + tmp_m2 + tmp_m3)/   &
    ((tmp_m1/dens_amm_sulf) + (tmp_m2/dens_amm_bisulf)   &
    + (tmp_m3/dens_sulf_acid))
  mass_part  = voldry_part*dens_part
  ! (kg dry aerosol)/(mol aerosol so4)
  kgaero_per_moleso4a = 1.0e-3_wp*(tmp_m1 + tmp_m2 + tmp_m3)
  ! correction when host code sulfate is really ammonium bisulfate/sulfate
  kgaero_per_moleso4a = kgaero_per_moleso4a * (mw_so4a_host/mw_so4a)

  ! fraction of wet volume due to so4a
  tmpb = 1.0_wp
  wet_volfrac_so4a = 1.0_wp / wetvol_dryvol

  ! calc kerminen & kulmala (2002) correction
  if (igrow <=  0) then
    factor_kk = 1.0_wp
  else
    ! "gr" parameter (nm/h) = condensation growth rate of new particles
    ! use kk2002 eqn 21 for h2so4 uptake, and correct for nh3 & h2o uptake
    tmp_spd = 14.7_wp*sqrt(temp_in)   ! h2so4 molecular speed (m/s)
    gr_kk = 3.0e-9_wp*tmp_spd*mw_sulf_acid*so4vol_in/(dens_part*wet_volfrac_so4a)

    ! "gamma" parameter (nm2/m2/h)
    ! use kk2002 eqn 22

    ! dfin_kk = wet diam (nm) of grown particle having dry dia = dpdry_part (m)
    dfin_kk = 1.0e9_wp * dpdry_part * (wetvol_dryvol**one_third)
    ! dnuc_kk = wet diam (nm) of cluster
    dnuc_kk = 2.0_wp*radius_cluster
    dnuc_kk = max( dnuc_kk, 1.0_wp )
    ! neglect (dmean/150)**0.048 factor,
    ! which should be very close to 1.0 because of small exponent
    gamma_kk = 0.23_wp * (dnuc_kk)**0.2_wp   &
      * (dfin_kk/3.0_wp)**0.075_wp   &
      * (dens_part*1.0e-3_wp)**(-0.33_wp)   &
      * (temp_in/293.0_wp)**(-0.75_wp)

    ! "cs_prime parameter" (1/m2)
    ! instead kk2002 eqn 3, use
    !     cs_prime ~= tmpa / (4*pi*tmpb * h2so4_accom_coef)
    ! where
    !     tmpa = -d(ln(h2so4))/dt by conden to particles   (1/h units)
    !     tmpb = h2so4 vapor diffusivity (m2/h units)
    ! this approx is generally within a few percent of the cs_prime
    !     calculated directly from eqn 2,
    !     which is acceptable, given overall uncertainties
    ! tmpa = -d(ln(h2so4))/dt by conden to particles   (1/h units)
    tmpa = h2so4_uptkrate * 3600.0_wp
    tmpa1 = tmpa
    tmpa = max( tmpa, 0.0_wp )
    ! tmpb = h2so4 gas diffusivity (m2/s, then m2/h)
    tmpb = 6.7037e-6_wp * (temp_in**0.75_wp) / cair
    tmpb1 = tmpb         ! m2/s
    tmpb = tmpb*3600.0_wp   ! m2/h
    cs_prime_kk = tmpa/(4.0_wp*pi*tmpb*accom_coef_h2so4)
    cs_kk = cs_prime_kk*4.0_wp*pi*tmpb1

    ! "nu" parameter (nm) -- kk2002 eqn 11
    nu_kk = gamma_kk*cs_prime_kk/gr_kk
    ! nucleation rate adjustment factor (--) -- kk2002 eqn 13
    factor_kk = exp( (nu_kk/dfin_kk) - (nu_kk/dnuc_kk) )
  end if
  J_star_kk = J_star_bb*factor_kk
#endif

end function

end module

