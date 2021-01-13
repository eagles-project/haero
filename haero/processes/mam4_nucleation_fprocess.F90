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

  !> min h2so4 vapor for nuc calcs = 4.0e-16 mol/mol-air ~= 1.0e4 molecules/cm3,
  real(wp), parameter :: qh2so4_cutoff = 4.0e-16_wp

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

! calculates new particle production from homogeneous nucleation
!    over timestep dtnuc, using nucleation rates from either
!    merikanto et al. (2007) h2so4-nh3-h2o ternary parameterization
!    vehkamaki et al. (2002) h2so4-h2o binary parameterization
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
! references:
!
!    vehkamäki, h., m. kulmala, i. napari, k.e.j. lehtinen,
!       c. timmreck, m. noppel and a. laaksonen, 2002,
!       an improved parameterization for sulfuric acid-water nucleation
!       rates for tropospheric and stratospheric conditions,
!       j. geophys. res., 107, 4622, doi:10.1029/2002jd002184
!
!    kerminen, v., and m. kulmala, 2002,
!  analytical formulae connecting the "real" and the "apparent"
!  nucleation rate and the nuclei number concentration
!  for atmospheric nucleation events
subroutine mer07_veh02_nuc_mosaic_1box(   &
  newnuc_method_flagaa, dtnuc, temp_in, rh_in, press_in,   &
  zm_in, pblh_in,   &
  qh2so4_cur, qh2so4_avg, qnh3_cur, h2so4_uptkrate,   &
  mw_so4a_host,   &
  nsize, maxd_asize, dplom_sect, dphim_sect,   &
  isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
  qh2so4_del, qnh3_del, dens_nh4so4a, ldiagaa,   &
  dnclusterdt )

  use mo_constants, only: rgas, &               ! Gas constant (J/K/kmol)
  avogad => avogadro    ! Avogadro's number (1/kmol)
  use physconst,    only: mw_so4a => mwso4, &   ! Molecular weight of sulfate
  mw_nh4a => mwnh4      ! Molecular weight of ammonium

  implicit none

  ! subr arguments (in)
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

  ! subr arguments (out)
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
  integer, save :: icase = 0, icase_reldiffmax = 0
  !       integer, parameter :: ldiagaa = -1
  integer :: lun
  integer :: newnuc_method_flagaa2

  real(wp), parameter :: onethird = 1.0_wp/3.0_wp

  real(wp), parameter :: accom_coef_h2so4 = 0.65_wp   ! accomodation coef for h2so4 conden

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

  real(wp), save :: reldiffmax = 0.0_wp

  real(wp) cair                     ! dry-air molar density (mol/m3)
  real(wp) cs_prime_kk              ! kk2002 "cs_prime" parameter (1/m2)
  real(wp) cs_kk                    ! kk2002 "cs" parameter (1/s)
  real(wp) dens_part                ! "grown" single-particle dry density (kg/m3)
  real(wp) dfin_kk, dnuc_kk         ! kk2002 final/initial new particle wet diameter (nm)
  real(wp) dpdry_clus               ! critical cluster diameter (m)
  real(wp) dpdry_part               ! "grown" single-particle dry diameter (m)
  real(wp) tmpa, tmpb, tmpc, tmpe, tmpq
  real(wp) tmpa1, tmpb1
  real(wp) tmp_m1, tmp_m2, tmp_m3, tmp_n1, tmp_n2, tmp_n3
  real(wp) tmp_spd                  ! h2so4 vapor molecular speed (m/s)
  real(wp) factor_kk
  real(wp) fogas, foso4a, fonh4a, fonuma
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
  if ( present ( dnclusterdt ) ) then
    dnclusterdt = 0.0_wp
  end if

  if ((newnuc_method_flagaa /=  1) .and. &
      (newnuc_method_flagaa /=  2) .and. &
      (newnuc_method_flagaa /= 11) .and. &
      (newnuc_method_flagaa /= 12)) then
    return
  end if

  ! call parameterization routine

  ! calc h2so4 in molecules/cm3 and nh3 in ppt
  cair = press_in/(temp_in*rgas)
  so4vol_in  = qh2so4_avg * cair * avogad * 1.0e-6_wp
  nh3ppt    = qnh3_cur * 1.0e12_wp
  ratenuclt = 1.0e-38_wp
  rateloge = log( ratenuclt )

  ! call vehkamaki binary parameterization routine

  if (so4vol_in >= 1.0e4_wp) then
    temp_bb = max( 230.15_wp, min( 305.15_wp, temp_in ) )
    rh_bb = max( 1.0e-4_wp, min( 1.0_wp, rh_in ) )
    so4vol_bb = max( 1.0e4_wp, min( 1.0e11_wp, so4vol_in ) )
    call binary_nuc_vehk2002(   &
      temp_bb, rh_bb, so4vol_bb,   &
      ratenuclt, rateloge,   &
      cnum_h2so4, cnum_tot, radius_cluster )
  end if

  cnum_nh3 = 0.0_wp
  newnuc_method_flagaa2 = 2

  rateloge = rateloge + log( max( 1.0e-38_wp, adjust_factor_bin_tern_ratenucl ) )

  ! do boundary layer nuc
  if ((newnuc_method_flagaa == 11) .or.   &
    (newnuc_method_flagaa == 12)) then
    if ( zm_in <= max(pblh_in,100.0_wp) ) then
      so4vol_bb = so4vol_in
      call pbl_nuc_wang2008( so4vol_bb,   &
        newnuc_method_flagaa, newnuc_method_flagaa2,   &
        ratenuclt, rateloge,   &
        cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster )
    end if
  end if

  ! if nucleation rate is less than 1e-6 #/cm3/s ~= 0.1 #/cm3/day,
  ! exit with new particle formation = 0
  if (rateloge  .le. -13.82_wp) then
    return
  end if

  ratenuclt = exp( rateloge )
  ratenuclt_bb = ratenuclt*1.0e6_wp  ! ratenuclt_bb is #/m3/s; ratenuclt is #/cm3/s
  if ( present ( dnclusterdt ) ) then
    dnclusterdt = ratenuclt_bb
  end if

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
  if (igrow .le. 0) then
    ! no "growing" so pure sulfuric acid
    tmp_n1 = 0.0_wp
    tmp_n2 = 0.0_wp
    tmp_n3 = 1.0_wp
  else if (qnh3_cur .ge. qh2so4_cur) then
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
  if (igrow <=  0) then
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
  if (qmolso4a_del_max .gt. qh2so4_cur) then
    freducea = qh2so4_cur/qmolso4a_del_max
  end if

  ! check if max production exceeds available nh3 vapor
  freduceb = 1.0_wp
  if (molenh4a_per_moleso4a .ge. 1.0e-10_wp) then
    ! max production of aerosol nh4 (ppm) based on ratenuclt_kk and mass_part
    qmolnh4a_del_max = qmolso4a_del_max*molenh4a_per_moleso4a
    if (qmolnh4a_del_max .gt. qnh3_cur) then
      freduceb = qnh3_cur/qmolnh4a_del_max
    end if
  end if
  freduce = min( freducea, freduceb )

  ! if adjusted nucleation rate is less than 1e-12 #/m3/s ~= 0.1 #/cm3/day,
  ! exit with new particle formation = 0
  if (freduce*ratenuclt_kk .le. 1.0e-12_wp) then
    return
  end if

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
end subroutine mer07_veh02_nuc_mosaic_1box

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
  if (tmp_rateloge <= rateloge) then
    return
  end if

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

end subroutine pbl_nuc_wang2008

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

  ! subr arguments (in)
  real(wp), intent(in) :: temp              ! temperature (k)
  real(wp), intent(in) :: rh                ! relative humidity (0-1)
  real(wp), intent(in) :: so4vol            ! concentration of h2so4 (molecules cm-3)

  ! subr arguments (out)
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
  real(wp) :: tmpa, tmpb

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

  tmpa     =     (   &
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
  tmpa = min( tmpa, log(1.0e38_wp) )
  ratenucl = exp ( tmpa )

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

  !   calc radius (nm) of critical cluster
  radius_cluster = exp( -1.6524245_wp + 0.42316402_wp*crit_x   &
    + 0.3346648_wp*log(cnum_tot) )

end subroutine binary_nuc_vehk2002

end module

  ! !PUBLIC DATA MEMBERS:
  integer, public :: newnuc_h2so4_conc_flag = 1

  ! adjustment factors
!!!  real(wp), parameter, public :: adjust_factor_dnaitdt = 1.0_wp       ! applied to final dnait/dt. Redefined in modal_aero_amicphys


  real(wp), parameter, public :: adjust_factor_bin_tern_ratenucl = 1.0_wp  !  applied to binary/ternary nucleation rate
  ! real(wp), parameter, public :: adjust_factor_pbl_ratenucl = 1.0_wp  ! applied to boundary layer nucleation rate
  real(wp),            public :: adjust_factor_pbl_ratenucl = 1.0_wp  ! applied to boundary layer nucleation rate
                                                                      ! value reassigned in amicphys


  ! !NON-PUBLIC DATA MEMBERS:
  integer, parameter  :: pcnstxx = gas_pcnst
  integer  :: l_h2so4_sv, l_nh3_sv, lnumait_sv, lnh4ait_sv, lso4ait_sv

  ! !DESCRIPTION: This module implements ...
  !
  ! !REVISION HISTORY:
  !
  !   R.Easter 2007.09.14:  Adapted from MIRAGE2 code
  !
  !EOP
  !----------------------------------------------------------------------
  !BOC

  ! list private module data here

  !EOC
  !----------------------------------------------------------------------


contains

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  !BOP
  ! !ROUTINE:  modal_aero_newnuc_sub --- ...
  !
  ! !INTERFACE:


!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine mam_newnuc_1subarea(                               &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
         latndx,            lonndx,           lund,                 &
         deltat,                                                    &
         temp,              pmid,             aircon,               &
         zmid,              pblh,             relhum,               &
         uptkrate_h2so4,   del_h2so4_gasprod, del_h2so4_aeruptk,    &
         n_mode,                                                    &
         qgas_cur,          qgas_avg,                               &
         qnum_cur,                                                  &
         qaer_cur,                                                  &
         qwtr_cur,                                                  &
         dnclusterdt                                                )

! uses
      use chem_mods,     only: adv_mass

      implicit none

! arguments
      integer,  intent(in) :: nstep                 ! model time-step number
      integer,  intent(in) :: lchnk                 ! chunk identifier
      integer,  intent(in) :: i, k                  ! column and level indices
      integer,  intent(in) :: jsub                  ! sub-area index
      integer,  intent(in) :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in) :: lund                  ! logical unit for diagnostic output
      integer,  intent(in) :: n_mode                ! current number of modes (including temporary)

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

      real(wp), intent(inout) :: dnclusterdt   ! cluster nucleation rate (#/m3/s)

      real(wp), intent(inout), dimension( 1:max_gas ) :: &
         qgas_cur
      real(wp), intent(in   ), dimension( 1:max_gas ) :: &
         qgas_avg
      real(wp), intent(inout), dimension( 1:max_mode ) :: &
         qnum_cur
      real(wp), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_cur
      real(wp), intent(inout), dimension( 1:max_mode ) :: &
         qwtr_cur

! DESCRIPTION:
!   computes changes due to aerosol nucleation (new particle formation)
!       treats both nucleation and subsequent growth of new particles
!          to aitken mode size
!   uses the following parameterizations
!       vehkamaki et al. (2002) parameterization for binary
!           homogeneous nucleation (h2so4-h2o) plus
!       kerminen and kulmala (2002) parameterization for
!           new particle loss during growth to aitken size
!
! REVISION HISTORY:
!   R.Easter 2007.09.14:  Adapted from MIRAGE2 code and CMAQ V4.6 code
!

! local variables
      integer, parameter :: ldiag1=-1, ldiag2=-1, ldiag3=-1, ldiag4=-1
      integer, parameter :: newnuc_method_flagaa = 11
!     integer, parameter :: newnuc_method_flagaa = 12
        !  1=merikanto et al (2007) ternary   2=vehkamaki et al (2002) binary
        ! 11=merikanto ternary + first-order boundary layer
        ! 12=merikanto ternary + second-order boundary layer

      integer :: itmp
      integer :: l
      integer :: ldiagveh02
      integer :: m

      real(wp) :: dens_nh4so4a
      real(wp) :: dmdt_ait, dmdt_aitsv1, dmdt_aitsv2, dmdt_aitsv3
      real(wp) :: dndt_ait, dndt_aitsv1, dndt_aitsv2, dndt_aitsv3
      real(wp) :: dnh4dt_ait, dso4dt_ait
      real(wp) :: dpnuc
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
      real(wp) :: tmp_q_del
      real(wp) :: tmp_frso4, tmp_uptkrate

      character(len=1) :: tmpch1, tmpch2, tmpch3


! begin
      dnclusterdt = 0.0_wp

! qh2so4_cur = current qh2so4, after aeruptk
! qh2so4_avg = average qh2so4 over time-step
      qh2so4_cur = qgas_cur(igas_h2so4)

      if ( (gaexch_h2so4_uptake_optaa == 1) .and. &
           (newnuc_h2so4_conc_optaa   == 1) ) then
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

!   qh2so4_avg = estimated average qh2so4
!   when production & loss are done simultaneously
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

!   dry-diameter limits for "grown" new particles
      dplom_mode(1) = exp( 0.67_wp*log(dgnumlo_aer(nait))   &
                         + 0.33_wp*log(dgnum_aer(nait)) )
      dphim_mode(1) = dgnumhi_aer(nait)

!   mass1p_... = mass (kg) of so4 & nh4 in a single particle of diameter ...
!                (assuming same dry density for so4 & nh4)
!      mass1p_aitlo - dp = dplom_mode(1)
!      mass1p_aithi - dp = dphim_mode(1)
      tmpa = dens_so4a_host*pi/6.0_wp
      mass1p_aitlo = tmpa*(dplom_mode(1)**3)
      mass1p_aithi = tmpa*(dphim_mode(1)**3)

!   limit RH to between 0.1% and 99%
      relhumnn = max( 0.01_wp, min( 0.99_wp, relhum ) )


!   call ... routine to get nucleation rates
      ldiagveh02 = -1
#ifdef l_debug_print
#if ( defined( MAM_STANDALONE ) )
      if (ldiag2 > 0) then
      if ((lonndx == 37) .and. (latndx == 23)) then
      if ((k >= 24) .or. (mod(k,4) == 0)) then
         ldiagveh02 = +1
         write(lund,'(/a,i8,3i4,f8.2,1p,4e10.2)')   &
            'veh02 call - nstep,lat,lon,k; tk,rh,p,cair',   &
            nstep, latndx, lonndx, k,   &
            temp, relhumnn, pmid, aircon*1.0e3_wp
         ! output aircon at (mol/m3)
      end if
      end if
      end if ! (ldiag2 > 0)
#endif
#endif

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
!----------------------------------------------------------------------
!       subr mer07_veh02_nuc_mosaic_1box(   &
!          newnuc_method_flagaa,   &
!          dtnuc, temp_in, rh_in, press_in,   &
!          qh2so4_cur, qh2so4_avg, qnh3_cur, h2so4_uptkrate,   &
!          nsize, maxd_asize, dplom_sect, dphim_sect,   &
!          isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
!          qh2so4_del, qnh3_del, dens_nh4so4a )
!
!! subr arguments (in)
!        real(wp), intent(in) :: dtnuc             ! nucleation time step (s)
!        real(wp), intent(in) :: temp_in           ! temperature, in k
!        real(wp), intent(in) :: rh_in             ! relative humidity, as fraction
!        real(wp), intent(in) :: press_in          ! air pressure (pa)
!
!        real(wp), intent(in) :: qh2so4_cur, qh2so4_avg
!                                                  ! gas h2so4 mixing ratios (mol/mol-air)
!        real(wp), intent(in) :: qnh3_cur          ! gas nh3 mixing ratios (mol/mol-air)
!             ! qxxx_cur = current value (after gas chem and condensation)
!             ! qxxx_avg = estimated average value (for simultaneous source/sink calcs)
!        real(wp), intent(in) :: h2so4_uptkrate    ! h2so4 uptake rate to aerosol (1/s)

!
!        integer, intent(in) :: nsize                    ! number of aerosol size bins
!        integer, intent(in) :: maxd_asize               ! dimension for dplom_sect, ...
!        real(wp), intent(in) :: dplom_sect(maxd_asize)  ! dry diameter at lower bnd of bin (m)
!        real(wp), intent(in) :: dphim_sect(maxd_asize)  ! dry diameter at upper bnd of bin (m)
!
!! subr arguments (out)
!        integer, intent(out) :: isize_nuc         ! size bin into which new particles go
!        real(wp), intent(out) :: qnuma_del        ! change to aerosol number mixing ratio (#/mol-air)
!        real(wp), intent(out) :: qso4a_del        ! change to aerosol so4 mixing ratio (mol/mol-air)
!        real(wp), intent(out) :: qnh4a_del        ! change to aerosol nh4 mixing ratio (mol/mol-air)
!        real(wp), intent(out) :: qh2so4_del       ! change to gas h2so4 mixing ratio (mol/mol-air)
!        real(wp), intent(out) :: qnh3_del         ! change to gas nh3 mixing ratio (mol/mol-air)
!                                                  ! aerosol changes are > 0; gas changes are < 0
!        real(wp), intent(out) :: dens_nh4so4a     ! dry-density of the new nh4-so4 aerosol mass (kg/m3)
!----------------------------------------------------------------------


!   convert qnuma_del from (#/mol-air) to (#/kmol-air)
      qnuma_del = qnuma_del*1.0e3_wp
#ifdef l_debug_print
#if ( defined( MAM_STANDALONE ) )
        if ( ldiag97 ) then
        write(lun97,'(/a,2i5,1p,3e11.3,2x,3e11.3)') &
           'newnuc - i, k, qavg s/n, uprt, rh 1/2', &
           i, k, qh2so4_avg, qh2so4_cur, qnh3_cur, &
           tmp_uptkrate, relhum, relhumnn
        write(lun97,'( a,10x,1p,3e11.3,2x,3e11.3)') &
           '               del qn, qso4a, qnh4a  ', &
           qnuma_del, qso4a_del, qnh4a_del
        end if
#endif
#endif

!   number nuc rate (#/kmol-air/s) from number nuc amt
      dndt_ait = qnuma_del/deltat

!   fraction of mass nuc going to so4
      tmpa = qso4a_del*mw_so4a_host
      if (igas_nh3 > 0) then
         tmpb = tmpa + qnh4a_del*mw_nh4a_host
         tmp_frso4 = max( tmpa, 1.0e-35_wp )/max( tmpb, 1.0e-35_wp )
      else
         tmpb = tmpa
         tmp_frso4 = 1.0_wp
      end if

!   mass nuc rate (kg/kmol-air/s) from mass nuc amts
      dmdt_ait = max( 0.0_wp, (tmpb/deltat) )

      dndt_aitsv1 = dndt_ait
      dmdt_aitsv1 = dmdt_ait
      dndt_aitsv2 = 0.0_wp
      dmdt_aitsv2 = 0.0_wp
      dndt_aitsv3 = 0.0_wp
      dmdt_aitsv3 = 0.0_wp
      tmpch1 = ' '
      tmpch2 = ' '

      if (dndt_ait < 1.0e2_wp) then
!   ignore newnuc if number rate < 100 #/kmol-air/s ~= 0.3 #/mg-air/d
         dndt_ait = 0.0_wp
         dmdt_ait = 0.0_wp
         tmpch1 = 'A'

      else
         dndt_aitsv2 = dndt_ait
         dmdt_aitsv2 = dmdt_ait
         tmpch1 = 'B'

!   mirage2 code checked for complete h2so4 depletion here,
!   but this is now done in mer07_veh02_nuc_mosaic_1box
         mass1p = dmdt_ait/dndt_ait
         dndt_aitsv3 = dndt_ait
         dmdt_aitsv3 = dmdt_ait

!   apply particle size constraints
         if (mass1p < mass1p_aitlo) then
!   reduce dndt to increase new particle size
            dndt_ait = dmdt_ait/mass1p_aitlo
            tmpch1 = 'C'
         else if (mass1p > mass1p_aithi) then
!   reduce dmdt to decrease new particle size
            dmdt_ait = dndt_ait*mass1p_aithi
            tmpch1 = 'E'
         end if
      end if

! *** apply adjustment factor to avoid unrealistically high
!     aitken number concentrations in mid and upper troposphere
      dndt_ait = dndt_ait * newnuc_adjust_factor_dnaitdt
      dmdt_ait = dmdt_ait * newnuc_adjust_factor_dnaitdt

      tmp_q_del = dndt_ait*deltat
      qnum_cur(     nait) = qnum_cur(     nait) + tmp_q_del

!   dso4dt_ait, dnh4dt_ait are (kmol/kmol-air/s)
      dso4dt_ait = dmdt_ait*tmp_frso4/mw_so4a_host
      dnh4dt_ait = dmdt_ait*(1.0_wp - tmp_frso4)/mw_nh4a_host

      if (dso4dt_ait > 0.0_wp) then
         tmp_q_del = dso4dt_ait*deltat
         qaer_cur(     iaer_so4,nait) = qaer_cur(     iaer_so4,nait) + tmp_q_del

         tmp_q_del = min( tmp_q_del, qgas_cur(igas_h2so4) )
         qgas_cur(     igas_h2so4) = qgas_cur(     igas_h2so4) - tmp_q_del
      end if

      if ((igas_nh3 > 0) .and. (dnh4dt_ait > 0.0_wp)) then
         tmp_q_del = dnh4dt_ait*deltat
         qaer_cur(     iaer_nh4,nait) = qaer_cur(     iaer_nh4,nait) + tmp_q_del

         tmp_q_del = min( tmp_q_del, qgas_cur(igas_nh3) )
         qgas_cur(     igas_nh3) = qgas_cur(     igas_nh3) - tmp_q_del
      end if

!!   temporary diagnostic
!        if (ldiag3 > 0) then
!        if ((dndt_ait /= 0.0_wp) .or. (dmdt_ait /= 0.0_wp)) then
!           write(lund,'(3a,1x,i7,3i5,1p,5e12.4)')   &
!              'newnucxx', tmpch1, tmpch2, nstep, lchnk, i, k,   &
!              dndt_ait, dmdt_ait, cldx
!!          call endrun( 'modal_aero_newnuc_sub' )
!        end if
!        end if


#ifdef l_debug_print
#if ( defined( CAMBOX_NEVER_ACTIVATE_THIS ) )
!   diagnostic output start ----------------------------------------
       if (ldiag4 > 0) then
       if ((lonndx == 37) .and. (latndx == 23)) then
       if ((k >= 24) .or. (mod(k,4) == 0)) then
        write(lund,97010) nstep, latndx, lonndx, k, temp, aircon*1.0e3_wp
        write(lund,97020) 'pmid                         ',   &
                pmid
        write(lund,97030) 'qv,qvsw,      rh_av, rh_clr  ',   &
                qv(i,k), qvswtr,       relhumav, relhum
        write(lund,97020) 'h2so4_cur,       _av, nh3_cur',   &
             qh2so4_cur,         qh2so4_avg, qnh3_cur
        write(lund,97020) 'del_h2so4_gasprod, _aeruptk  ',   &
             del_h2so4_gasprod(i,k), del_h2so4_aeruptk(i,k),   &
             tmp_uptkrate*3600.0
        write(lund,97020) ' '
        write(lund,97050) 'tmpch1, tmpch2               ', tmpch1, tmpch2
        write(lund,97020) 'dndt_, dmdt_aitsv1           ',   &
                          dndt_aitsv1, dmdt_aitsv1
        write(lund,97020) 'dndt_, dmdt_aitsv2           ',   &
                          dndt_aitsv2, dmdt_aitsv2
        write(lund,97020) 'dndt_, dmdt_aitsv3           ',   &
                          dndt_aitsv3, dmdt_aitsv3
        write(lund,97020) 'dndt_, dmdt_ait              ',   &
                          dndt_ait, dmdt_ait
        write(lund,97020) 'dso4dt_, dnh4dt_ait          ',   &
                          dso4dt_ait, dnh4dt_ait
        write(lund,97020) 'qso4a_del, qh2so4_del        ',   &
                          qso4a_del, qh2so4_del
        write(lund,97020) 'qnh4a_del, qnh3_del          ',   &
                          qnh4a_del, qnh3_del
        write(lund,97020) 'dqdt(h2so4), (nh3)           ',   &
              dqdt(i,k,l_h2so4), dqdt(i,k,l_nh3)
        write(lund,97020) 'dqdt(so4a), (nh4a), (numa)   ',   &
              dqdt(i,k,lso4ait), dqdt(i,k,lnh4ait), dqdt(i,k,lnumait)

       dpnuc = 0.0_wp
       if (dndt_aitsv1 > 1.0e-5_wp) dpnuc = (6.0_wp*dmdt_aitsv1/   &
                   (pi*dens_so4a_host*dndt_aitsv1))**0.3333333_wp
       if (dpnuc > 0.0_wp) then
       write(lund,97020) 'dpnuc,      dp_aitlo, _aithi ',   &
                    dpnuc, dplom_mode(1), dphim_mode(1)
       write(lund,97020) 'mass1p, mass1p_aitlo, _aithi ',   &
                    mass1p, mass1p_aitlo, mass1p_aithi
       end if

97010  format( / 'NEWNUC nstep,lat,lon,k,tk,cair', i8, 3i4, f8.2, 1pe12.4 )
97020  format( a, 1p, 6e12.4 )
97030  format( a, 1p, 2e12.4, 0p, 5f10.6 )
97040  format( 29x, 1p, 6e12.4 )
97050  format( a, 2(3x,a) )
       end if ! ((k >= 24) .or. (mod(k,4) == 0))
       end if ! ((lonndx == 37) .and. (latndx == 23))
       end if ! (ldiag4 > 0)
!   diagnostic output end   ------------------------------------------
#endif
#endif


80000 continue


      return
      end subroutine mam_newnuc_1subarea

