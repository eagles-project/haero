!> This module implements MAM's nucleation process. For details, see the
!> appropriate section in the Processes chapter of the Haero design document.
module mam_nucleation

  use haero_precision, only: wp
  use haero, only: model_t, aerosol_species_t, gas_species_t, &
                   prognostics_t, atmosphere_t, diagnostics_t, tendencies_t, &
                   var_not_found
  use haero_constants, only: pi, R_gas, Avogadro

  implicit none
  private

  ! Module functions
  public :: init, run, finalize, ternary_nuc_merik2007, binary_nuc_vehk2002, &
    pbl_nuc_wang2008, mer07_veh02_nuc_mosaic_1box
  public :: adjust_factor_pbl_ratenucl
  public :: adjust_factor_bin_tern_ratenucl

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

  !-------------------------
  ! Module global variables
  !-------------------------

  ! These variables are adapted from MAM's modal_aero_microp_control module
  ! and initialized in the init subroutine. Here they are variables, not
  ! properties, because their values depend on the particular model
  ! configuration.

  !> The index of the Aitken mode
  integer :: nait

  !> The geometric mean particle diameters for all aerosol modes
  real(wp), dimension(:), allocatable :: dgnum_aer

  !> The minimum particle diameters for all aerosol modes
  real(wp), dimension(:), allocatable :: dgnumlo_aer

  !> The maximum particle diameters for all aerosol modes
  real(wp), dimension(:), allocatable :: dgnumhi_aer

  !> Gas mixing ratios buffer
  real(wp), dimension(:), allocatable :: qgas_cur

  !> Time-averaged gas mixing ratios buffer
  real(wp), dimension(:), allocatable :: qgas_avg

  !> Modal number concentration buffer
  real(wp), dimension(:), allocatable :: qnum_cur

  !> Modal aerosol mixing ratios buffer.
  real(wp), dimension(:,:), allocatable :: qaer_cur

  !> Modal water content buffer.
  real(wp), dimension(:), allocatable :: qwtr_cur

  !> This is an option flag for H2SO4 uptake.
  integer :: gaexch_h2so4_uptake_optaa

  !> Index of H2SO4 gas
  integer :: igas_h2so4

  !> Index of NH3 gas
  integer :: igas_nh3

  !> Index of NH4 aerosol within the Aitken mode
  integer :: iaer_nh4

  !> Index of SO4 aerosol within the Aitken mode
  integer :: iaer_so4

  !> The mass density of SO4 aerosol as assumed by the host atm model
  real(wp) :: dens_so4a_host

  !> The molecular weight of NH4 aerosol as assumed by the host atm model
  real(wp) :: mw_nh4a_host

  !> The molecular weight of SO4 aerosol as assumed by the host atm model
  real(wp) :: mw_so4a_host

  !> Controls treatment of H2SO4 condensation
  !>    1 = sequential   calc. of gas-chem prod then condensation loss
  !>    2 = simultaneous calc. of gas-chem prod and  condensation loss
  integer :: newnuc_h2so4_conc_optaa

  !> Adjustment factor for Aitken number concentration tendency
  real(wp) :: newnuc_adjust_factor_dnaitdt

  !> Adjustment factor for nucleation rate with binary/ternary nucleation.
  real(wp) :: adjust_factor_bin_tern_ratenucl

  !> Adjustment factor for nucleation rate corrected for the planetary boundary
  !> layer.
  real(wp) :: adjust_factor_pbl_ratenucl

contains

subroutine init(model)

  implicit none

  ! Arguments
  type(model_t), intent(in) :: model

  type(aerosol_species_t) so4, nh4
  integer :: m

  ! Extract mode properties.
  allocate(dgnum_aer(model%num_modes))
  allocate(dgnumlo_aer(model%num_modes))
  allocate(dgnumhi_aer(model%num_modes))
  do m = 1,model%num_modes
    dgnum_aer(m) = model%modes(m)%mean_std_dev
    dgnumlo_aer(m) = model%modes(m)%min_diameter
    dgnumhi_aer(m) = model%modes(m)%max_diameter
  end do

  ! Allocate gas and aerosol state buffers.
  allocate(qgas_cur(model%num_gases))
  allocate(qgas_avg(model%num_gases))
  allocate(qnum_cur(model%num_modes))
  allocate(qaer_cur(maxval(model%num_mode_species), model%num_modes))
  allocate(qwtr_cur(model%num_modes))

  ! Record the aitken mode index.
  nait = model%mode_index("aitken")

  ! Record the index for SO4 aerosol within the Aitken mode and fetch some
  ! properties.
  iaer_so4 = model%aerosol_index(nait, "SO4")
  if (iaer_so4 > 0) then
    so4 = model%aero_species(nait, iaer_so4)
    mw_so4a_host = so4%molecular_wt
  end if

  ! Record the index for NH4 aerosol within the Aitken mode and fetch some
  ! properties.
  iaer_nh4 = model%aerosol_index(nait, "NH4")
  if (iaer_nh4 > 0) then
    nh4 = model%aero_species(nait, iaer_nh4)
    mw_nh4a_host = nh4%molecular_wt
  end if

  ! Record the index for H2SO4 gas (source of new nuclei).
  igas_h2so4 = model%gas_index("H2SO4")

  ! Set some defaults
  gaexch_h2so4_uptake_optaa = 2
  newnuc_h2so4_conc_optaa = 2
  dens_so4a_host = 1.0_wp ! FIXME

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
  integer :: k, p, m, s, token
  real(wp), pointer, dimension(:,:) :: q_i    ! interstitial aerosol mix ratios
  real(wp), pointer, dimension(:,:) :: q_g    ! gas mix ratios
  real(wp), pointer, dimension(:,:) :: n      ! modal number concentrations

  real(wp), pointer, dimension(:) :: temp     ! atmospheric temperature
  real(wp), pointer, dimension(:) :: press    ! atmospheric pressure
  real(wp), pointer, dimension(:) :: rel_hum  ! atmospheric relative humidity
  real(wp), pointer, dimension(:) :: height   ! atmospheric height

  real(wp), pointer, dimension(:,:) :: dqdt_i ! interstitial aerosol tends
  real(wp), pointer, dimension(:,:) :: dqdt_g ! gas mole frac tends
  real(wp), pointer, dimension(:,:) :: dndt   ! modal number density tends

  ! Diagnostics computed by other processes.
  real(wp), pointer, dimension(:,:) :: qgas_averaged
  real(wp), pointer, dimension(:) :: uptkrate_h2so4
  real(wp), pointer, dimension(:) :: del_h2so4_gasprod
  real(wp), pointer, dimension(:) :: del_h2so4_aeruptk

  real(wp) :: aircon    ! molar concentration of air [mol/m^3]
  real(wp) :: pblh      ! Planetary boundary layer height [m]

  real(wp) :: h2so4_uptake_rate, h2so4_gasprod_change, h2so4_aeruptk_change
  real(wp) :: dndt_ait, dmdt_ait, dso4dt_ait, dnh4dt_ait
  real(wp) :: dnclusterdt ! diagnostic cluster nucleation rate (#/m3/s)

  ! First of all, check to make sure our model has an aitken mode. If it
  ! doesn't, we can return immediately.
  if (nait == 0) then
    return
  end if

  ! If there's no gas present with which to create new nuclei, there's
  ! nothing to do, either.
  if ((igas_h2so4 == 0) .and. (igas_nh3 == 0)) then
    return
  end if

  ! Finally, if there are no relevant aerosol species for nuclei, we can't
  ! create them.
  if ((iaer_so4 == 0) .and. (iaer_nh4 == 0)) then
    return
  end if

  ! Gas mole fraction tendencies
  q_g => prognostics%gases()
  dqdt_g => tendencies%gases()

  ! Mix fractions and tendencies for SO4 aerosol in the Aitken mode
  ! All new nuclei are deposited into interstitial aerosols.
  q_i => prognostics%interstitial_aerosols()
  dqdt_i => tendencies%interstitial_aerosols()

  ! Modal number density and tendencies
  n => prognostics%modal_num_concs()
  dndt => tendencies%modal_num_concs()

  ! Atmospheric state variables
  press => atmosphere%pressure()
  temp => atmosphere%temperature()
  rel_hum => atmosphere%relative_humidity()
  height => atmosphere%height()
  pblh = atmosphere%planetary_boundary_height()

  ! Diagnostics
  token = diagnostics%find_gas_var("qgas_averaged")
  if (token == var_not_found) then
    qgas_averaged => null()
  else
    qgas_averaged => diagnostics%gas_var(token)
  end if

  token = diagnostics%find_var("uptkrate_h2so4")
  if (token == var_not_found) then
    uptkrate_h2so4 => null()
  else
    uptkrate_h2so4 => diagnostics%var(token)
  end if

  token = diagnostics%find_var("del_h2so4_gasprod")
  if (token == var_not_found) then
    del_h2so4_gasprod => null()
  else
    del_h2so4_gasprod => diagnostics%var(token)
  end if

  token = diagnostics%find_var("del_h2so4_aeruptk")
  if (token == var_not_found) then
    del_h2so4_aeruptk => null()
  else
    del_h2so4_aeruptk => diagnostics%var(token)
  end if

  ! Traverse the vertical levels and compute tendencies from nucleation.
  do k = 1,model%num_levels
    ! Compute the molar concentration of air at the given pressure and
    ! temperature.
    aircon = press(k)/(temp(k)*R_gas)

    ! Extract prognostic state data.
    qgas_cur(:) = q_g(k, :)
    qnum_cur(:) = n(k, :)
    do p = 1,model%num_populations
      call model%get_mode_and_species(p, m, s)
      qaer_cur(s,m) = q_i(k, p)
    end do

    ! Extract diagnostic state data.
    if (associated(qgas_averaged)) then
      qgas_avg(:) = qgas_averaged(k, :)
    else
      qgas_avg(:) = 0_wp
    end if
    if (associated(uptkrate_h2so4)) then
      h2so4_uptake_rate = uptkrate_h2so4(k)
    else
      h2so4_uptake_rate = 0_wp
    end if
    if (associated(del_h2so4_gasprod)) then
      h2so4_gasprod_change = del_h2so4_gasprod(k)
    else
      h2so4_gasprod_change = 0_wp
    end if
    if (associated(del_h2so4_aeruptk)) then
      h2so4_aeruptk_change = del_h2so4_aeruptk(k)
    else
      h2so4_aeruptk_change = 0_wp
    end if
    qwtr_cur(:) = 0_wp !n(k, :) ! FIXME: Need to compute water content.

    call compute_tendencies(dt, &
      temp(k), press(k), aircon, height(k), pblh, rel_hum(k), &
      h2so4_uptake_rate, h2so4_gasprod_change, h2so4_aeruptk_change, &
      qgas_cur, qgas_avg, qnum_cur, qaer_cur, qwtr_cur, &
      dndt_ait, dmdt_ait, dso4dt_ait, dnh4dt_ait, &
      dnclusterdt)

    dqdt_i(k, :) = 0_wp
    dqdt_i(k, iaer_so4) = dso4dt_ait
    dqdt_g(k, :) = 0_wp
    dqdt_g(k, igas_h2so4) = -dso4dt_ait
    dndt(k, :) = 0_wp
    dndt(k, nait) = dndt_ait
  end do
end subroutine

subroutine finalize(model)
  implicit none

  ! Arguments
  type(model_t), intent(in) :: model

  ! Deallocate gas and aerosol state buffers
  deallocate(qgas_cur)
  deallocate(qgas_avg)
  deallocate(qnum_cur)
  deallocate(qaer_cur)
  deallocate(qwtr_cur)

  ! Deallocate mode metadata.
  deallocate(dgnum_aer)
  deallocate(dgnumlo_aer)
  deallocate(dgnumhi_aer)

end subroutine

! Computes tendencies due to aerosol nucleation (new particle formation).
! Treats both nucleation and subsequent growth of new particles to aitken
! mode size. Uses the following parameterizations:
! * Vehkamaki et al. (2002), parameterization for binary homogeneous nucleation
!   (h2so4-h2o)
! * Kerminen and Kulmala (2002), parameterization for new particle loss during
!   growth to aitken size
subroutine compute_tendencies(deltat, &
  temp, pmid, aircon, zmid, pblh, relhum, &
  uptkrate_h2so4, del_h2so4_gasprod, del_h2so4_aeruptk, &
  qgas_cur, qgas_avg, qnum_cur, qaer_cur, qwtr_cur, &
  dndt_ait, dmdt_ait, dso4dt_ait, dnh4dt_ait, &
  dnclusterdt)

  implicit none

  ! arguments
  real(wp), intent(in) :: deltat           ! model timestep (s)
  real(wp), intent(in) :: temp             ! temperature (K)
  real(wp), intent(in) :: pmid             ! pressure at model levels (Pa)
  real(wp), intent(in) :: aircon           ! air molar concentration (kmol/m3)
  real(wp), intent(in) :: zmid             ! midpoint height above surface (m)
  real(wp), intent(in) :: pblh             ! pbl height (m)
  real(wp), intent(in) :: relhum           ! relative humidity (0-1)
  real(wp), intent(in) :: uptkrate_h2so4
  real(wp), intent(in) :: del_h2so4_gasprod
  real(wp), intent(in) :: del_h2so4_aeruptk

  real(wp), intent(in), dimension(:) :: qgas_cur
  real(wp), intent(in), dimension(:) :: qgas_avg
  real(wp), intent(in), dimension(:) :: qnum_cur
  real(wp), intent(in), dimension(:,:) :: qaer_cur
  real(wp), intent(in), dimension(:) :: qwtr_cur
  real(wp), intent(out) :: dndt_ait, dmdt_ait, dnh4dt_ait, dso4dt_ait
  real(wp), intent(inout) :: dnclusterdt   ! cluster nucleation rate (#/m3/s)

  integer :: itmp
  integer :: ldiagveh02

  real(wp) :: dens_nh4so4a
  real(wp) :: dmdt_aitsv1, dmdt_aitsv2, dmdt_aitsv3
  real(wp) :: dndt_aitsv1, dndt_aitsv2, dndt_aitsv3
  real(wp) :: dplom_mode(1), dphim_mode(1)
  real(wp) :: mass1p
  real(wp) :: mass1p_aithi, mass1p_aitlo
  real(wp) :: qh2so4_cur, qh2so4_avg, qh2so4_del
  real(wp) :: qnh3_cur, qnh3_del, qnh4a_del
  real(wp) :: qnuma_del
  real(wp) :: qso4a_del
  real(wp) :: relhumnn
  real(wp) :: tmpa, tmpb, tmpc
  real(wp) :: tmp_q2, tmp_q3
  real(wp) :: tmp_frso4, tmp_uptkrate

  dndt_ait = 0.0_wp
  dmdt_ait = 0.0_wp
  dnh4dt_ait = 0.0_wp
  dso4dt_ait = 0.0_wp
  dnclusterdt = 0.0_wp

  ! qh2so4_cur = current qh2so4, after aeruptk
  ! qh2so4_avg = average qh2so4 over time-step
  qh2so4_cur = qgas_cur(igas_h2so4)

  if ((gaexch_h2so4_uptake_optaa == 1) .and. &
      (newnuc_h2so4_conc_optaa   == 1)) then
    ! estimate qh2so4_avg using the method in standard cam5.2 modal_aero_newnuc

    ! skip if h2so4 vapor < qh2so4_cutoff
    if (qh2so4_cur <= qh2so4_cutoff) goto 80000

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
    ! d[ln(qh2so4)]/dt (1/s) from uptake (condensation) to aerosol
    tmp_uptkrate = tmpb/deltat

    ! qh2so4_avg = estimated average qh2so4
    ! when production & loss are done simultaneously
    if (tmpb <= 0.1_wp) then
      qh2so4_avg = tmp_q3*(1.0_wp + 0.5_wp*tmpb) - 0.5_wp*tmpa
    else
      tmpc = tmpa/tmpb
      qh2so4_avg = (tmp_q3 - tmpc)*((exp(tmpb)-1.0_wp)/tmpb) + tmpc
    end if
  else
    ! use qh2so4_avg and first-order loss rate calculated in mam_gasaerexch_1subarea
    qh2so4_avg = qgas_avg(igas_h2so4)
    tmp_uptkrate = uptkrate_h2so4
  end if

  if (qh2so4_avg <= qh2so4_cutoff) goto 80000

  if (igas_nh3 > 0) then
    qnh3_cur = max( 0.0_wp, qgas_cur(igas_nh3) )
  else
    qnh3_cur = 0.0_wp
  end if

  ! dry-diameter limits for "grown" new particles
  dplom_mode(1) = exp( 0.67_wp*log(dgnumlo_aer(nait))   &
                     + 0.33_wp*log(dgnum_aer(nait)) )
  dphim_mode(1) = dgnumhi_aer(nait)

  ! mass1p_... = mass (kg) of so4 & nh4 in a single particle of diameter ...
  ! (assuming same dry density for so4 & nh4)
  ! mass1p_aitlo - dp = dplom_mode(1)
  ! mass1p_aithi - dp = dphim_mode(1)
  tmpa = dens_so4a_host*pi/6.0_wp
  mass1p_aitlo = tmpa*(dplom_mode(1)**3)
  mass1p_aithi = tmpa*(dphim_mode(1)**3)

  ! limit RH to between 0.1% and 99%
  relhumnn = max( 0.01_wp, min( 0.99_wp, relhum ) )

  ! call routine to get nucleation rates
  ldiagveh02 = -1
  call mer07_veh02_nuc_mosaic_1box(   &
    newnuc_method_flagaa,   &
    deltat, temp, relhumnn, pmid,   &
    zmid, pblh,   &
    qh2so4_cur, qh2so4_avg, qnh3_cur, tmp_uptkrate,   &
    mw_so4a_host,   &
    1, 1, dplom_mode, dphim_mode,   &
    itmp, qnuma_del, qso4a_del, qnh4a_del,   &
    qh2so4_del, qnh3_del, dens_nh4so4a,   &
    ldiagveh02, dnclusterdt )

  ! convert qnuma_del from (#/mol-air) to (#/kmol-air)
  qnuma_del = qnuma_del*1.0e3_wp

  ! number nuc rate (#/kmol-air/s) from number nuc amt
  dndt_ait = qnuma_del/deltat

  ! fraction of mass nuc going to so4
  tmpa = qso4a_del*mw_so4a_host
  if (igas_nh3 > 0) then
    tmpb = tmpa + qnh4a_del*mw_nh4a_host
    tmp_frso4 = max( tmpa, 1.0e-35_wp )/max( tmpb, 1.0e-35_wp )
  else
    tmpb = tmpa
    tmp_frso4 = 1.0_wp
  end if

  ! mass nuc rate (kg/kmol-air/s) from mass nuc amts
  dmdt_ait = max( 0.0_wp, (tmpb/deltat) )

  dndt_aitsv1 = dndt_ait
  dmdt_aitsv1 = dmdt_ait
  dndt_aitsv2 = 0.0_wp
  dmdt_aitsv2 = 0.0_wp
  dndt_aitsv3 = 0.0_wp
  dmdt_aitsv3 = 0.0_wp

  if (dndt_ait < 1.0e2_wp) then
    ! ignore newnuc if number rate < 100 #/kmol-air/s ~= 0.3 #/mg-air/d
    dndt_ait = 0.0_wp
    dmdt_ait = 0.0_wp
  else
    dndt_aitsv2 = dndt_ait
    dmdt_aitsv2 = dmdt_ait

    ! mirage2 code checked for complete h2so4 depletion here,
    ! but this is now done in mer07_veh02_nuc_mosaic_1box
    mass1p = dmdt_ait/dndt_ait
    dndt_aitsv3 = dndt_ait
    dmdt_aitsv3 = dmdt_ait

    ! apply particle size constraints
    if (mass1p < mass1p_aitlo) then
      ! reduce dndt to increase new particle size
      dndt_ait = dmdt_ait/mass1p_aitlo
    else if (mass1p > mass1p_aithi) then
      ! reduce dmdt to decrease new particle size
      dmdt_ait = dndt_ait*mass1p_aithi
    end if
  end if

  ! *** apply adjustment factor to avoid unrealistically high
  ! aitken number concentrations in mid and upper troposphere
  dndt_ait = dndt_ait * newnuc_adjust_factor_dnaitdt
  dmdt_ait = dmdt_ait * newnuc_adjust_factor_dnaitdt

  80000 continue

  ! dso4dt_ait, dnh4dt_ait are (kmol/kmol-air/s)
  dso4dt_ait = dmdt_ait*tmp_frso4/mw_so4a_host
  if (0.0_wp < mw_nh4a_host) then
    dnh4dt_ait = dmdt_ait*(1.0_wp - tmp_frso4)/mw_nh4a_host
  end if
end subroutine

! Calculates new particle production from homogeneous nucleation
! over timestep dtnuc, using nucleation rates from either
! Merikanto et al. (2007) h2so4-nh3-h2o ternary parameterization
! Vehkamaki et al. (2002) h2so4-h2o binary parameterization
!
! the new particles are "grown" to the lower-bound size of the host code's
!    smallest size bin.  (this "growth" is somewhat ad hoc, and would not be
!    necessary if the host code's size bins extend down to ~1 nm.)
!
!    if the h2so4 and nh3 mass mixing ratios (mixrats) of the grown new
!    particles exceed the current gas mixrats, the new particle production
!    is reduced so that the new particle mass mixrats match the gas mixrats.
!
!    the correction of kerminen and kulmala (2002) is applied to account
!    for loss of the new particles by coagulation as they are
!    growing to the "host code mininum size"
!
! References:
! * merikanto, j., i. napari, h. vehkamaki, t. anttila,
!   and m. kulmala, 2007, new parameterization of
!   sulfuric acid-ammonia-water ternary nucleation
!   rates at tropospheric conditions,
!   j. geophys. res., 112, d15207, doi:10.1029/2006jd0027977
!
! * vehkamäki, h., m. kulmala, i. napari, k.e.j. lehtinen,
!   c. timmreck, m. noppel and a. laaksonen, 2002,
!   an improved parameterization for sulfuric acid-water nucleation
!   rates for tropospheric and stratospheric conditions,
!   j. geophys. res., 107, 4622, doi:10.1029/2002jd002184
!
! * kerminen, v., and m. kulmala, 2002,
!   analytical formulae connecting the "real" and the "apparent"
!   nucleation rate and the nuclei number concentration
!   for atmospheric nucleation events
subroutine mer07_veh02_nuc_mosaic_1box(   &
  newnuc_method_flagaa, dtnuc, temp_in, rh_in, press_in,   &
  zm_in, pblh_in,   &
  qh2so4_cur, qh2so4_avg, qnh3_cur, h2so4_uptkrate,   &
  mw_so4a_host,   &
  nsize, maxd_asize, dplom_sect, dphim_sect,   &
  isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
  qh2so4_del, qnh3_del, dens_nh4so4a, ldiagaa,   &
  dnclusterdt )

  use haero_constants, only: &
    rgas => r_gas, &                ! Gas constant (J/K/kmol)
    avogad => avogadro, &           ! Avogadro's number (1/kmol)
    mw_so4a => molec_weight_so4, &  ! Molecular weight of sulfate
    mw_nh4a => molec_weight_nh4     ! Molecular weight of ammonium

  implicit none

  ! arguments (in)
  real(wp), intent(in) :: dtnuc             ! nucleation time step (s)
  real(wp), intent(in) :: temp_in           ! temperature, in k
  real(wp), intent(in) :: rh_in             ! relative humidity, as fraction
  real(wp), intent(in) :: press_in          ! air pressure (pa)
  real(wp), intent(in) :: zm_in             ! layer midpoint height (m)
  real(wp), intent(in) :: pblh_in           ! pbl height (m)
  real(wp), intent(in) :: qh2so4_cur, qh2so4_avg
  ! gas h2so4 mixing ratios (mol/mol-air)
  real(wp), intent(in) :: qnh3_cur          ! gas nh3 mixing ratios (mol/mol-air)
  ! qxxx_cur = current value (after gas chem and condensation)
  ! qxxx_avg = estimated average value (for simultaneous source/sink calcs)
  real(wp), intent(in) :: h2so4_uptkrate    ! h2so4 uptake rate to aerosol (1/s)
  real(wp), intent(in) :: mw_so4a_host      ! mw of so4 aerosol in host code (g/mol)

  integer, intent(in) :: newnuc_method_flagaa     ! 1=merikanto et al (2007) ternary
  ! 2=vehkamaki et al (2002) binary
  integer, intent(in) :: nsize                    ! number of aerosol size bins
  integer, intent(in) :: maxd_asize               ! dimension for dplom_sect, ...
  real(wp), intent(in) :: dplom_sect(maxd_asize)  ! dry diameter at lower bnd of bin (m)
  real(wp), intent(in) :: dphim_sect(maxd_asize)  ! dry diameter at upper bnd of bin (m)
  integer, intent(in) :: ldiagaa

  ! arguments (out)
  integer, intent(out) :: isize_nuc         ! size bin into which new particles go
  real(wp), intent(out) :: qnuma_del        ! change to aerosol number mixing ratio (#/mol-air)
  real(wp), intent(out) :: qso4a_del        ! change to aerosol so4 mixing ratio (mol/mol-air)
  real(wp), intent(out) :: qnh4a_del        ! change to aerosol nh4 mixing ratio (mol/mol-air)
  real(wp), intent(out) :: qh2so4_del       ! change to gas h2so4 mixing ratio (mol/mol-air)
  real(wp), intent(out) :: qnh3_del         ! change to gas nh3 mixing ratio (mol/mol-air)
  ! aerosol changes are > 0; gas changes are < 0
  real(wp), intent(out) :: dens_nh4so4a     ! dry-density of the new nh4-so4 aerosol mass (kg/m3)
  real(wp), intent(out), optional :: dnclusterdt ! cluster nucleation rate (#/m3/s)

  ! subr arguments (out) passed via common block
  !    these are used to duplicate the outputs of yang zhang's original test driver
  !    they are not really needed in wrf-chem
  real(wp) :: ratenuclt        ! j = ternary nucleation rate from napari param. (cm-3 s-1)
  real(wp) :: rateloge         ! ln (j)
  real(wp) :: cnum_h2so4       ! number of h2so4 molecules in the critical nucleus
  real(wp) :: cnum_nh3         ! number of nh3   molecules in the critical nucleus
  real(wp) :: cnum_tot         ! total number of molecules in the critical nucleus
  real(wp) :: radius_cluster   ! the radius of cluster (nm)

  ! local variables
  integer :: i
  integer :: igrow
  integer :: newnuc_method_flagaa2

  real(wp) cair                     ! dry-air molar density (mol/m3)
  real(wp) cs_prime_kk              ! kk2002 "cs_prime" parameter (1/m2)
  real(wp) cs_kk                    ! kk2002 "cs" parameter (1/s)
  real(wp) dens_part                ! "grown" single-particle dry density (kg/m3)
  real(wp) dfin_kk, dnuc_kk         ! kk2002 final/initial new particle wet diameter (nm)
  real(wp) dpdry_clus               ! critical cluster diameter (m)
  real(wp) dpdry_part               ! "grown" single-particle dry diameter (m)
  real(wp) tmpa, tmpb, tmpc, tmpe
  real(wp) tmpa1, tmpb1
  real(wp) tmp_m1, tmp_m2, tmp_m3, tmp_n1, tmp_n2, tmp_n3
  real(wp) tmp_spd                  ! h2so4 vapor molecular speed (m/s)
  real(wp) factor_kk
  real(wp) freduce                  ! reduction factor applied to nucleation rate
  ! due to limited availability of h2so4 & nh3 gases
  real(wp) freducea, freduceb
  real(wp) gamma_kk                 ! kk2002 "gamma" parameter (nm2*m2/h)
  real(wp) gr_kk                    ! kk2002 "gr" parameter (nm/h)
  real(wp) kgaero_per_moleso4a      ! (kg dry aerosol)/(mol aerosol so4)
  real(wp) mass_part                ! "grown" single-particle dry mass (kg)
  real(wp) molenh4a_per_moleso4a    ! (mol aerosol nh4)/(mol aerosol so4)
  real(wp) nh3ppt, nh3ppt_bb        ! actual and bounded nh3 (ppt)
  real(wp) nu_kk                    ! kk2002 "nu" parameter (nm)
  real(wp) qmolnh4a_del_max         ! max production of aerosol nh4 over dtnuc (mol/mol-air)
  real(wp) qmolso4a_del_max         ! max production of aerosol so4 over dtnuc (mol/mol-air)
  real(wp) ratenuclt_bb             ! nucleation rate (#/m3/s)
  real(wp) ratenuclt_kk             ! nucleation rate after kk2002 adjustment (#/m3/s)
  real(wp) rh_bb                    ! bounded value of rh_in
  real(wp) so4vol_in                ! concentration of h2so4 for nucl. calc., molecules cm-3
  real(wp) so4vol_bb                ! bounded value of so4vol_in
  real(wp) temp_bb                  ! bounded value of temp_in
  real(wp) voldry_clus              ! critical-cluster dry volume (m3)
  real(wp) voldry_part              ! "grown" single-particle dry volume (m3)
  real(wp) wetvol_dryvol            ! grown particle (wet-volume)/(dry-volume)
  real(wp) wet_volfrac_so4a         ! grown particle (dry-volume-from-so4)/(wet-volume)

  ! if h2so4 vapor < qh2so4_cutoff
  ! exit with new particle formation = 0
  isize_nuc = 1
  qnuma_del = 0.0_wp
  qso4a_del = 0.0_wp
  qnh4a_del = 0.0_wp
  qh2so4_del = 0.0_wp
  qnh3_del = 0.0_wp
  if ( present ( dnclusterdt ) ) dnclusterdt = 0.0_wp

  if ((newnuc_method_flagaa /=  1) .and. &
      (newnuc_method_flagaa /=  2) .and. &
      (newnuc_method_flagaa /= 11) .and. &
      (newnuc_method_flagaa /= 12)) return

  ! make call to parameterization routine

  ! calc h2so4 in molecules/cm3 and nh3 in ppt
  cair = press_in/(temp_in*rgas)
  so4vol_in  = qh2so4_avg * cair * avogad * 1.0e-6_wp
  nh3ppt    = qnh3_cur * 1.0e12_wp
  ratenuclt = 1.0e-38_wp
  rateloge = log( ratenuclt )

  if ( (newnuc_method_flagaa /=  2) .and. (nh3ppt >= 0.1_wp) ) then
    ! make call to merikanto ternary parameterization routine
    ! (when nh3ppt < 0.1, use binary param instead)

    if (so4vol_in >= 5.0e4_wp) then
      temp_bb = max( 235.0_wp, min( 295.0_wp, temp_in ) )
      rh_bb = max( 0.05_wp, min( 0.95_wp, rh_in ) )
      so4vol_bb = max( 5.0e4_wp, min( 1.0e9_wp, so4vol_in ) )
      nh3ppt_bb = max( 0.1_wp, min( 1.0e3_wp, nh3ppt ) )
      call ternary_nuc_merik2007(   &
        temp_bb, rh_bb, so4vol_bb, nh3ppt_bb, rateloge, &
        cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster)
    end if
    newnuc_method_flagaa2 = 1

  else
    ! make call to vehkamaki binary parameterization routine

    if (so4vol_in >= 1.0e4_wp) then
      temp_bb = max( 230.15_wp, min( 305.15_wp, temp_in ) )
      rh_bb = max( 1.0e-4_wp, min( 1.0_wp, rh_in ) )
      so4vol_bb = max( 1.0e4_wp, min( 1.0e11_wp, so4vol_in ) )
      call binary_nuc_vehk2002(temp_bb, rh_bb, so4vol_bb, &
        ratenuclt, rateloge, cnum_h2so4, cnum_tot, radius_cluster )
    end if
    cnum_nh3 = 0.0_wp
    newnuc_method_flagaa2 = 2
  end if

  rateloge  = rateloge + log(max(1.0e-38_wp, adjust_factor_bin_tern_ratenucl))

  ! do boundary layer nuc
  if ((newnuc_method_flagaa == 11) .or.   &
      (newnuc_method_flagaa == 12)) then
    if ( zm_in <= max(pblh_in,100.0_wp) ) then
      so4vol_bb = so4vol_in
      call pbl_nuc_wang2008( so4vol_bb, newnuc_method_flagaa, &
        newnuc_method_flagaa2, ratenuclt, rateloge, cnum_tot, cnum_h2so4, &
        cnum_nh3, radius_cluster )
    end if
  end if

  ! if nucleation rate is less than 1e-6 #/cm3/s ~= 0.1 #/cm3/day,
  ! exit with new particle formation = 0
  if (rateloge <= -13.82_wp) return

  ratenuclt = exp( rateloge )
  ratenuclt_bb = ratenuclt*1.0e6_wp  ! ratenuclt_bb is #/m3/s; ratenuclt is #/cm3/s
  if ( present ( dnclusterdt ) ) dnclusterdt = ratenuclt_bb

  ! wet/dry volume ratio - use simple kohler approx for ammsulf/ammbisulf
  tmpa = max( 0.10_wp, min( 0.95_wp, rh_in ) )
  wetvol_dryvol = 1.0_wp - 0.56_wp/log(tmpa)

  ! determine size bin into which the new particles go
  ! (probably it will always be bin #1, but ...)
  voldry_clus = ( max(cnum_h2so4,1.0_wp)*mw_so4a + cnum_nh3*mw_nh4a ) /   &
    (1.0e3_wp*dens_sulfacid*avogad)
  ! correction when host code sulfate is really ammonium bisulfate/sulfate
  voldry_clus = voldry_clus * (mw_so4a_host/mw_so4a)
  dpdry_clus = (voldry_clus*6.0_wp/pi)**onethird

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

  ! determine composition and density of the "grown particles"
  ! the grown particles are assumed to be liquid
  !    (since critical clusters contain water)
  !    so any (nh4/so4) molar ratio between 0 and 2 is allowed
  ! assume that the grown particles will have
  !    (nh4/so4 molar ratio) = min( 2, (nh3/h2so4 gas molar ratio) )
  if (igrow <= 0) then
    ! no "growing" so pure sulfuric acid
    tmp_n1 = 0.0_wp
    tmp_n2 = 0.0_wp
    tmp_n3 = 1.0_wp
  else if (qnh3_cur >= qh2so4_cur) then
    ! combination of ammonium sulfate and ammonium bisulfate
    ! tmp_n1 & tmp_n2 = mole fractions of the ammsulf & ammbisulf
    tmp_n1 = (qnh3_cur/qh2so4_cur) - 1.0_wp
    tmp_n1 = max( 0.0_wp, min( 1.0_wp, tmp_n1 ) )
    tmp_n2 = 1.0_wp - tmp_n1
    tmp_n3 = 0.0_wp
  else
    ! combination of ammonium bisulfate and sulfuric acid
    ! tmp_n2 & tmp_n3 = mole fractions of the ammbisulf & sulfacid
    tmp_n1 = 0.0_wp
    tmp_n2 = (qnh3_cur/qh2so4_cur)
    tmp_n2 = max( 0.0_wp, min( 1.0_wp, tmp_n2 ) )
    tmp_n3 = 1.0_wp - tmp_n2
  end if

  tmp_m1 = tmp_n1*mw_ammsulf
  tmp_m2 = tmp_n2*mw_ammbisulf
  tmp_m3 = tmp_n3*mw_sulfacid
  dens_part = (tmp_m1 + tmp_m2 + tmp_m3)/   &
    ((tmp_m1/dens_ammsulf) + (tmp_m2/dens_ammbisulf)   &
    + (tmp_m3/dens_sulfacid))
  dens_nh4so4a = dens_part
  mass_part  = voldry_part*dens_part
  ! (mol aerosol nh4)/(mol aerosol so4)
  molenh4a_per_moleso4a = 2.0_wp*tmp_n1 + tmp_n2
  ! (kg dry aerosol)/(mol aerosol so4)
  kgaero_per_moleso4a = 1.0e-3_wp*(tmp_m1 + tmp_m2 + tmp_m3)
  ! correction when host code sulfate is really ammonium bisulfate/sulfate
  kgaero_per_moleso4a = kgaero_per_moleso4a * (mw_so4a_host/mw_so4a)

  ! fraction of wet volume due to so4a
  tmpb = 1.0_wp + molenh4a_per_moleso4a*17.0_wp/98.0_wp
  wet_volfrac_so4a = 1.0_wp / ( wetvol_dryvol * tmpb )

  ! calc kerminen & kulmala (2002) correction
  if (igrow <= 0) then
    factor_kk = 1.0_wp
  else
    ! "gr" parameter (nm/h) = condensation growth rate of new particles
    ! use kk2002 eqn 21 for h2so4 uptake, and correct for nh3 & h2o uptake
    tmp_spd = 14.7_wp*sqrt(temp_in)   ! h2so4 molecular speed (m/s)
    gr_kk = 3.0e-9_wp*tmp_spd*mw_sulfacid*so4vol_in/(dens_part*wet_volfrac_so4a)

    ! "gamma" parameter (nm2/m2/h)
    ! use kk2002 eqn 22
    ! dfin_kk = wet diam (nm) of grown particle having dry dia = dpdry_part (m)
    dfin_kk = 1.0e9_wp * dpdry_part * (wetvol_dryvol**onethird)
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
  ratenuclt_kk = ratenuclt_bb*factor_kk

  ! max production of aerosol dry mass (kg-aero/m3-air)
  tmpa = max( 0.0_wp, (ratenuclt_kk*dtnuc*mass_part) )
  ! max production of aerosol so4 (mol-so4a/mol-air)
  tmpe = tmpa/(kgaero_per_moleso4a*cair)
  ! max production of aerosol so4 (mol/mol-air)
  ! based on ratenuclt_kk and mass_part
  qmolso4a_del_max = tmpe

  ! check if max production exceeds available h2so4 vapor
  freducea = 1.0_wp
  if (qmolso4a_del_max > qh2so4_cur) then
    freducea = qh2so4_cur/qmolso4a_del_max
  end if

  ! check if max production exceeds available nh3 vapor
  freduceb = 1.0_wp
  if (molenh4a_per_moleso4a >= 1.0e-10_wp) then
    ! max production of aerosol nh4 (ppm) based on ratenuclt_kk and mass_part
    qmolnh4a_del_max = qmolso4a_del_max*molenh4a_per_moleso4a
    if (qmolnh4a_del_max > qnh3_cur) then
      freduceb = qnh3_cur/qmolnh4a_del_max
    end if
  end if
  freduce = min( freducea, freduceb )

  ! if adjusted nucleation rate is less than 1e-12 #/m3/s ~= 0.1 #/cm3/day,
  ! exit with new particle formation = 0
  if (freduce*ratenuclt_kk <= 1.0e-12_wp) return

  ! note:  suppose that at this point, freduce < 1.0 (no gas-available
  !    constraints) and molenh4a_per_moleso4a < 2.0
  ! if the gas-available constraints is do to h2so4 availability,
  !    then it would be possible to condense "additional" nh3 and have
  !    (nh3/h2so4 gas molar ratio) < (nh4/so4 aerosol molar ratio) <= 2
  ! one could do some additional calculations of
  !    dens_part & molenh4a_per_moleso4a to realize this
  ! however, the particle "growing" is a crude approximate way to get
  !    the new particles to the host code's minimum particle size,
  ! are such refinements worth the effort?

  ! changes to h2so4 & nh3 gas (in mol/mol-air), limited by amounts available
  tmpa = 0.9999_wp
  qh2so4_del = min( tmpa*qh2so4_cur, freduce*qmolso4a_del_max )
  qnh3_del   = min( tmpa*qnh3_cur, qh2so4_del*molenh4a_per_moleso4a )
  qh2so4_del = -qh2so4_del
  qnh3_del   = -qnh3_del

  ! changes to so4 & nh4 aerosol (in mol/mol-air)
  qso4a_del = -qh2so4_del
  qnh4a_del =   -qnh3_del
  ! change to aerosol number (in #/mol-air)
  qnuma_del = 1.0e-3_wp*(qso4a_del*mw_so4a + qnh4a_del*mw_nh4a)/mass_part

  ! do the following (tmpa, tmpb, tmpc) calculations as a check
  ! max production of aerosol number (#/mol-air)
  tmpa = max( 0.0_wp, (ratenuclt_kk*dtnuc/cair) )
  ! adjusted production of aerosol number (#/mol-air)
  tmpb = tmpa*freduce
  ! relative difference from qnuma_del
  tmpc = (tmpb - qnuma_del)/max(tmpb, qnuma_del, 1.0e-35_wp)
end subroutine

! calculates boundary nucleation nucleation rate
! using the first or second-order parameterization in
!     wang, m., and j.e. penner, 2008,
!        aerosol indirect forcing in a global model with particle nucleation,
!        atmos. chem. phys. discuss., 8, 13943-13998
subroutine pbl_nuc_wang2008( so4vol,   &
  newnuc_method_flagaa, newnuc_method_flagaa2,   &
  ratenucl, rateloge,   &
  cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster ) 

  implicit none

  ! subr arguments (in)
  real(wp), intent(in) :: so4vol            ! concentration of h2so4 (molecules cm-3)
  integer, intent(in)  :: newnuc_method_flagaa
  ! [11,12] value selects [first,second]-order parameterization

  ! subr arguments (inout)
  integer, intent(inout)  :: newnuc_method_flagaa2
  real(wp), intent(inout) :: ratenucl         ! binary nucleation rate, j (# cm-3 s-1)
  real(wp), intent(inout) :: rateloge         ! log( ratenucl )
  real(wp), intent(inout) :: cnum_tot         ! total number of molecules
  ! in the critical nucleus
  real(wp), intent(inout) :: cnum_h2so4       ! number of h2so4 molecules
  real(wp), intent(inout) :: cnum_nh3         ! number of nh3 molecules
  real(wp), intent(inout) :: radius_cluster   ! the radius of cluster (nm)

  ! local variables
  real(wp) :: tmp_diam, tmp_mass, tmp_volu
  real(wp) :: tmp_rateloge, tmp_ratenucl

  ! nucleation rate
  if (newnuc_method_flagaa == 11) then
    tmp_ratenucl = 1.0e-6_wp * so4vol
  else if (newnuc_method_flagaa == 12) then
    tmp_ratenucl = 1.0e-12_wp * (so4vol**2)
  else
    return
  end if

  tmp_ratenucl = tmp_ratenucl * adjust_factor_pbl_ratenucl
  tmp_rateloge = log( max( 1.0e-38_wp, tmp_ratenucl ) )

  ! exit if pbl nuc rate is lower than (incoming) ternary/binary rate
  if (tmp_rateloge <= rateloge) return

  rateloge = tmp_rateloge
  ratenucl = tmp_ratenucl
  newnuc_method_flagaa2 = newnuc_method_flagaa

  ! following wang 2002, assume fresh nuclei are 1 nm diameter
  !    subsequent code will "grow" them to aitken mode size
  radius_cluster = 0.5_wp

  ! assume fresh nuclei are pure h2so4
  !    since aitken size >> initial size, the initial composition
  !    has very little impact on the results
  tmp_diam = radius_cluster * 2.0e-7_wp   ! diameter in cm
  tmp_volu = (tmp_diam**3) * (pi/6.0_wp)  ! volume in cm^3
  tmp_mass = tmp_volu * 1.8_wp            ! mass in g
  cnum_h2so4 = (tmp_mass / 98.0_wp) * 6.023e23_wp   ! no. of h2so4 molec assuming pure h2so4
  cnum_tot = cnum_h2so4
  cnum_nh3 = 0.0_wp
end subroutine

! calculates binary nucleation rate and critical cluster size
! using the parameterization in
!     vehkamäki, h., m. kulmala, i. napari, k.e.j. lehtinen,
!        c. timmreck, m. noppel and a. laaksonen, 2002,
!        an improved parameterization for sulfuric acid-water nucleation
!        rates for tropospheric and stratospheric conditions,
!        j. geophys. res., 107, 4622, doi:10.1029/2002jd002184
subroutine binary_nuc_vehk2002( temp, rh, so4vol,   &
  ratenucl, rateloge,   &
  cnum_h2so4, cnum_tot, radius_cluster )

  implicit none

  ! arguments (in)
  real(wp), intent(in) :: temp              ! temperature (k)
  real(wp), intent(in) :: rh                ! relative humidity (0-1)
  real(wp), intent(in) :: so4vol            ! concentration of h2so4 (molecules cm-3)

  ! arguments (out)
  real(wp), intent(out) :: ratenucl         ! binary nucleation rate, j (# cm-3 s-1)
  real(wp), intent(out) :: rateloge         ! log( ratenucl )

  real(wp), intent(out) :: cnum_h2so4       ! number of h2so4 molecules
  ! in the critical nucleus
  real(wp), intent(out) :: cnum_tot         ! total number of molecules
  ! in the critical nucleus
  real(wp), intent(out) :: radius_cluster   ! the radius of cluster (nm)


  ! local variables
  real(wp) :: crit_x
  real(wp) :: acoe, bcoe, ccoe, dcoe, ecoe, fcoe, gcoe, hcoe, icoe, jcoe
  real(wp) :: tmpa

  ! calc sulfuric acid mole fraction in critical cluster
  crit_x = 0.740997_wp - 0.00266379_wp * temp   &
    - 0.00349998_wp * log (so4vol)   &
    + 0.0000504022_wp * temp * log (so4vol)   &
    + 0.00201048_wp * log (rh)   &
    - 0.000183289_wp * temp * log (rh)   &
    + 0.00157407_wp * (log (rh)) ** 2.0_wp   &
    - 0.0000179059_wp * temp * (log (rh)) ** 2.0_wp   &
    + 0.000184403_wp * (log (rh)) ** 3.0_wp   &
    - 1.50345e-6_wp * temp * (log (rh)) ** 3.0_wp

  ! calc nucleation rate
  acoe = 0.14309_wp+2.21956_wp*temp   &
    - 0.0273911_wp * temp**2.0_wp   &
    + 0.0000722811_wp * temp**3.0_wp + 5.91822_wp/crit_x

  bcoe = 0.117489_wp + 0.462532_wp *temp   &
    - 0.0118059_wp * temp**2.0_wp   &
    + 0.0000404196_wp * temp**3.0_wp + 15.7963_wp/crit_x

  ccoe = -0.215554_wp-0.0810269_wp * temp   &
    + 0.00143581_wp * temp**2.0_wp   &
    - 4.7758e-6_wp * temp**3.0_wp   &
    - 2.91297_wp/crit_x

  dcoe = -3.58856_wp+0.049508_wp * temp   &
    - 0.00021382_wp * temp**2.0_wp   &
    + 3.10801e-7_wp * temp**3.0_wp   &
    - 0.0293333_wp/crit_x

  ecoe = 1.14598_wp - 0.600796_wp * temp   &
    + 0.00864245_wp * temp**2.0_wp   &
    - 0.0000228947_wp * temp**3.0_wp   &
    - 8.44985_wp/crit_x

  fcoe = 2.15855_wp + 0.0808121_wp * temp   &
    - 0.000407382_wp * temp**2.0_wp   &
    - 4.01957e-7_wp * temp**3.0_wp   &
    + 0.721326_wp/crit_x

  gcoe = 1.6241_wp - 0.0160106_wp * temp   &
    + 0.0000377124_wp * temp**2.0_wp   &
    + 3.21794e-8_wp * temp**3.0_wp   &
    - 0.0113255_wp/crit_x

  hcoe = 9.71682_wp - 0.115048_wp * temp   &
    + 0.000157098_wp * temp**2.0_wp   &
    + 4.00914e-7_wp * temp**3.0_wp   &
    + 0.71186_wp/crit_x

  icoe = -1.05611_wp + 0.00903378_wp * temp   &
    - 0.0000198417_wp * temp**2.0_wp   &
    + 2.46048e-8_wp  * temp**3.0_wp   &
    - 0.0579087_wp/crit_x

  jcoe = -0.148712_wp + 0.00283508_wp * temp   &
    - 9.24619e-6_wp  * temp**2.0_wp   &
    + 5.00427e-9_wp * temp**3.0_wp   &
    - 0.0127081_wp/crit_x

  tmpa = (   &
    acoe   &
    + bcoe * log (rh)   &
    + ccoe * ( log (rh))**2.0_wp   &
    + dcoe * ( log (rh))**3.0_wp   &
    + ecoe * log (so4vol)   &
    + fcoe * (log (rh)) * (log (so4vol))   &
    + gcoe * ((log (rh) ) **2.0_wp)   &
    * (log (so4vol))   &
    + hcoe * (log (so4vol)) **2.0_wp   &
    + icoe * log (rh)   &
    * ((log (so4vol)) **2.0_wp)   &
    + jcoe * (log (so4vol)) **3.0_wp   &
  )
  rateloge = tmpa
  if (log(1.0e38_wp) < tmpa) then
    print*,"Error: tmpa exceeds limit of about 87."
    stop
  end if
  ratenucl = exp ( tmpa )

  ! calc number of molecules in critical cluster
  acoe = -0.00295413_wp - 0.0976834_wp*temp   &
    + 0.00102485_wp * temp**2.0_wp   &
    - 2.18646e-6_wp * temp**3.0_wp - 0.101717_wp/crit_x

  bcoe = -0.00205064_wp - 0.00758504_wp*temp   &
    + 0.000192654_wp * temp**2.0_wp   &
    - 6.7043e-7_wp * temp**3.0_wp - 0.255774_wp/crit_x

  ccoe = +0.00322308_wp + 0.000852637_wp * temp   &
    - 0.0000154757_wp * temp**2.0_wp   &
    + 5.66661e-8_wp * temp**3.0_wp   &
    + 0.0338444_wp/crit_x

  dcoe = +0.0474323_wp - 0.000625104_wp * temp   &
    + 2.65066e-6_wp * temp**2.0_wp   &
    - 3.67471e-9_wp * temp**3.0_wp   &
    - 0.000267251_wp/crit_x

  ecoe = -0.0125211_wp + 0.00580655_wp * temp   &
    - 0.000101674_wp * temp**2.0_wp   &
    + 2.88195e-7_wp * temp**3.0_wp   &
    + 0.0942243_wp/crit_x

  fcoe = -0.038546_wp - 0.000672316_wp * temp   &
    + 2.60288e-6_wp * temp**2.0_wp   &
    + 1.19416e-8_wp * temp**3.0_wp   &
    - 0.00851515_wp/crit_x

  gcoe = -0.0183749_wp + 0.000172072_wp * temp   &
    - 3.71766e-7_wp * temp**2.0_wp   &
    - 5.14875e-10_wp * temp**3.0_wp   &
    + 0.00026866_wp/crit_x

  hcoe = -0.0619974_wp + 0.000906958_wp * temp   &
    - 9.11728e-7_wp * temp**2.0_wp   &
    - 5.36796e-9_wp * temp**3.0_wp   &
    - 0.00774234_wp/crit_x

  icoe = +0.0121827_wp - 0.00010665_wp * temp   &
    + 2.5346e-7_wp * temp**2.0_wp   &
    - 3.63519e-10_wp * temp**3.0_wp   &
    + 0.000610065_wp/crit_x

  jcoe = +0.000320184_wp - 0.0000174762_wp * temp   &
    + 6.06504e-8_wp * temp**2.0_wp   &
    - 1.4177e-11_wp * temp**3.0_wp   &
    + 0.000135751_wp/crit_x

  cnum_tot = exp (   &
    acoe   &
    + bcoe * log (rh)   &
    + ccoe * ( log (rh))**2.0_wp   &
    + dcoe * ( log (rh))**3.0_wp   &
    + ecoe * log (so4vol)   &
    + fcoe * (log (rh)) * (log (so4vol))   &
    + gcoe * ((log (rh) ) **2.0_wp)   &
    * (log (so4vol))   &
    + hcoe * (log (so4vol)) **2.0_wp   &
    + icoe * log (rh)   &
    * ((log (so4vol)) **2.0_wp)   &
    + jcoe * (log (so4vol)) **3.0_wp   &
  )

  cnum_h2so4 = cnum_tot * crit_x

  ! calc radius (nm) of critical cluster
  radius_cluster = exp( -1.6524245_wp + 0.42316402_wp*crit_x   &
    + 0.3346648_wp*log(cnum_tot) )
end subroutine

! fortran 90 subroutine that calculates the parameterized composition
! and nucleation rate of critical clusters in h2o-h2so4-nh3 vapor
!
! warning: the fit should not be used outside its limits of validity
! (limits indicated below)
!
! in:
! t:     temperature (k), limits 235-295 k
! rh:    relative humidity as fraction (eg. 0.5=50%) limits 0.05-0.95
! c2:    sulfuric acid concentration (molecules/cm3) limits 5x10^4 - 10^9 molecules/cm3
! c3:    ammonia mixing ratio (ppt) limits 0.1 - 1000 ppt
!
! out:
! j_log: logarithm of nucleation rate (1/(s cm3))
! ntot:  total number of molecules in the critical cluster
! nacid: number of sulfuric acid molecules in the critical cluster
! namm:  number of ammonia molecules in the critical cluster
! r:     radius of the critical cluster (nm)
subroutine ternary_nuc_merik2007(t, rh, c2, c3, j_log, ntot, nacid, namm, r )
  implicit none

  real(wp), intent(in) :: t, rh, c2, c3
  real(wp), intent(out) :: j_log, ntot, nacid, namm, r
  real(wp) :: t_onset

  t_onset=143.6002929064716_wp + 1.0178856665693992_wp*rh + &
    10.196398812974294_wp*log(c2) - &
    0.1849879416839113_wp*log(c2)**2 - 17.161783213150173_wp*log(c3) + &
    (109.92469248546053_wp*log(c3))/log(c2) + &
    0.7734119613144357_wp*log(c2)*log(c3) - 0.15576469879527022_wp*log(c3)**2

  if (t_onset > t) then
    j_log = -12.861848898625231_wp + 4.905527742256349_wp*c3 - 358.2337705052991_wp*rh -&
      0.05463019231872484_wp*c3*t + 4.8630382337426985_wp*rh*t + &
      0.00020258394697064567_wp*c3*t**2 - 0.02175548069741675_wp*rh*t**2 - &
      2.502406532869512e-7_wp*c3*t**3 + 0.00003212869941055865_wp*rh*t**3 - &
      4.39129415725234e6_wp/log(c2)**2 + (56383.93843154586_wp*t)/log(c2)**2 -&
      (239.835990963361_wp*t**2)/log(c2)**2 + &
      (0.33765136625580167_wp*t**3)/log(c2)**2 - &
      (629.7882041830943_wp*rh)/(c3**3*log(c2)) + &
      (7.772806552631709_wp*rh*t)/(c3**3*log(c2)) - &
      (0.031974053936299256_wp*rh*t**2)/(c3**3*log(c2)) + &
      (0.00004383764128775082_wp*rh*t**3)/(c3**3*log(c2)) + &
      1200.472096232311_wp*log(c2) - 17.37107890065621_wp*t*log(c2) + &
      0.08170681335921742_wp*t**2*log(c2) - 0.00012534476159729881_wp*t**3*log(c2) - &
      14.833042158178936_wp*log(c2)**2 + 0.2932631303555295_wp*t*log(c2)**2 - &
      0.0016497524241142845_wp*t**2*log(c2)**2 + &
      2.844074805239367e-6_wp*t**3*log(c2)**2 - 231375.56676032578_wp*log(c3) - &
      100.21645273730675_wp*rh*log(c3) + 2919.2852552424706_wp*t*log(c3) + &
      0.977886555834732_wp*rh*t*log(c3) - 12.286497122264588_wp*t**2*log(c3) - &
      0.0030511783284506377_wp*rh*t**2*log(c3) + &
      0.017249301826661612_wp*t**3*log(c3) + 2.967320346100855e-6_wp*rh*t**3*log(c3) + &
      (2.360931724951942e6_wp*log(c3))/log(c2) - &
      (29752.130254319443_wp*t*log(c3))/log(c2) + &
      (125.04965118142027_wp*t**2*log(c3))/log(c2) - &
      (0.1752996881934318_wp*t**3*log(c3))/log(c2) + &
      5599.912337254629_wp*log(c2)*log(c3) - 70.70896612937771_wp*t*log(c2)*log(c3) + &
      0.2978801613269466_wp*t**2*log(c2)*log(c3) - &
      0.00041866525019504_wp*t**3*log(c2)*log(c3) + 75061.15281456841_wp*log(c3)**2 - &
      931.8802278173565_wp*t*log(c3)**2 + 3.863266220840964_wp*t**2*log(c3)**2 - &
      0.005349472062284983_wp*t**3*log(c3)**2 - &
      (732006.8180571689_wp*log(c3)**2)/log(c2) + &
      (9100.06398573816_wp*t*log(c3)**2)/log(c2) - &
      (37.771091915932004_wp*t**2*log(c3)**2)/log(c2) + &
      (0.05235455395566905_wp*t**3*log(c3)**2)/log(c2) - &
      1911.0303773001353_wp*log(c2)*log(c3)**2 + &
      23.6903969622286_wp*t*log(c2)*log(c3)**2 - &
      0.09807872005428583_wp*t**2*log(c2)*log(c3)**2 + &
      0.00013564560238552576_wp*t**3*log(c2)*log(c3)**2 - &
      3180.5610833308_wp*log(c3)**3 + 39.08268568672095_wp*t*log(c3)**3 - &
      0.16048521066690752_wp*t**2*log(c3)**3 + &
      0.00022031380023793877_wp*t**3*log(c3)**3 + &
      (40751.075322248245_wp*log(c3)**3)/log(c2) - &
      (501.66977622013934_wp*t*log(c3)**3)/log(c2) + &
      (2.063469732254135_wp*t**2*log(c3)**3)/log(c2) - &
      (0.002836873785758324_wp*t**3*log(c3)**3)/log(c2) + &
      2.792313345723013_wp*log(c2)**2*log(c3)**3 - &
      0.03422552111802899_wp*t*log(c2)**2*log(c3)**3 + &
      0.00014019195277521142_wp*t**2*log(c2)**2*log(c3)**3 - &
      1.9201227328396297e-7_wp*t**3*log(c2)**2*log(c3)**3 - &
      980.923146020468_wp*log(rh) + 10.054155220444462_wp*t*log(rh) - &
      0.03306644502023841_wp*t**2*log(rh) + 0.000034274041225891804_wp*t**3*log(rh) + &
      (16597.75554295064_wp*log(rh))/log(c2) - &
      (175.2365504237746_wp*t*log(rh))/log(c2) + &
      (0.6033215603167458_wp*t**2*log(rh))/log(c2) - &
      (0.0006731787599587544_wp*t**3*log(rh))/log(c2) - &
      89.38961120336789_wp*log(c3)*log(rh) + 1.153344219304926_wp*t*log(c3)*log(rh) - &
      0.004954549700267233_wp*t**2*log(c3)*log(rh) + &
      7.096309866238719e-6_wp*t**3*log(c3)*log(rh) + &
      3.1712136610383244_wp*log(c3)**3*log(rh) - &
      0.037822330602328806_wp*t*log(c3)**3*log(rh) + &
      0.0001500555743561457_wp*t**2*log(c3)**3*log(rh) - &
      1.9828365865570703e-7_wp*t**3*log(c3)**3*log(rh)

    ntot = 57.40091052369212_wp - 0.2996341884645408_wp*t + &
      0.0007395477768531926_wp*t**2 - &
      5.090604835032423_wp*log(c2) + 0.011016634044531128_wp*t*log(c2) + &
      0.06750032251225707_wp*log(c2)**2 - 0.8102831333223962_wp*log(c3) + &
      0.015905081275952426_wp*t*log(c3) - 0.2044174683159531_wp*log(c2)*log(c3) + &
      0.08918159167625832_wp*log(c3)**2 - 0.0004969033586666147_wp*t*log(c3)**2 + &
      0.005704394549007816_wp*log(c3)**3 + 3.4098703903474368_wp*j_log - &
      0.014916956508210809_wp*t*j_log + 0.08459090011666293_wp*log(c3)*j_log - &
      0.00014800625143907616_wp*t*log(c3)*j_log + 0.00503804694656905_wp*j_log**2

    r = 3.2888553966535506e-10_wp - 3.374171768439839e-12_wp*t + &
      1.8347359507774313e-14_wp*t**2 + 2.5419844298881856e-12_wp*log(c2) - &
      9.498107643050827e-14_wp*t*log(c2) + 7.446266520834559e-13_wp*log(c2)**2 + &
      2.4303397746137294e-11_wp*log(c3) + 1.589324325956633e-14_wp*t*log(c3) - &
      2.034596219775266e-12_wp*log(c2)*log(c3) - 5.59303954457172e-13_wp*log(c3)**2 - &
      4.889507104645867e-16_wp*t*log(c3)**2 + 1.3847024107506764e-13_wp*log(c3)**3 + &
      4.141077193427042e-15_wp*j_log - 2.6813110884009767e-14_wp*t*j_log + &
      1.2879071621313094e-12_wp*log(c3)*j_log - &
      3.80352446061867e-15_wp*t*log(c3)*j_log - 1.8790172502456827e-14_wp*j_log**2

    nacid = -4.7154180661803595_wp + 0.13436423483953885_wp*t - &
      0.00047184686478816176_wp*t**2 - &
      2.564010713640308_wp*log(c2) + 0.011353312899114723_wp*t*log(c2) + &
      0.0010801941974317014_wp*log(c2)**2 + 0.5171368624197119_wp*log(c3) - &
      0.0027882479896204665_wp*t*log(c3) + 0.8066971907026886_wp*log(c3)**2 - &
      0.0031849094214409335_wp*t*log(c3)**2 - 0.09951184152927882_wp*log(c3)**3 + &
      0.00040072788891745513_wp*t*log(c3)**3 + 1.3276469271073974_wp*j_log - &
      0.006167654171986281_wp*t*j_log - 0.11061390967822708_wp*log(c3)*j_log + &
      0.0004367575329273496_wp*t*log(c3)*j_log + 0.000916366357266258_wp*j_log**2

    namm = 71.20073903979772_wp - 0.8409600103431923_wp*t + &
      0.0024803006590334922_wp*t**2 + &
      2.7798606841602607_wp*log(c2) - 0.01475023348171676_wp*t*log(c2) + &
      0.012264508212031405_wp*log(c2)**2 - 2.009926050440182_wp*log(c3) + &
      0.008689123511431527_wp*t*log(c3) - 0.009141180198955415_wp*log(c2)*log(c3) + &
      0.1374122553905617_wp*log(c3)**2 - 0.0006253227821679215_wp*t*log(c3)**2 + &
      0.00009377332742098946_wp*log(c3)**3 + 0.5202974341687757_wp*j_log - &
      0.002419872323052805_wp*t*j_log + 0.07916392322884074_wp*log(c3)*j_log - &
      0.0003021586030317366_wp*t*log(c3)*j_log + 0.0046977006608603395_wp*j_log**2

  else
    ! nucleation rate less that 5e-6, setting j_log arbitrary small
    j_log = -300._wp
  end if
end subroutine

end module

