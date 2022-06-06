module mam_gasaerexch
!-----------------------------------------------------------------------------------------
! Purpose:
!  Parameterization of gas-aerosol exchange. This microphysical process is
!  also referred to as the condensation/evaporation of chemical species
!  or the uptake of gas-phase species by aerosol particles.
!
!  This module contains
!    - condensation-specific parameters and their initialization,
!    - a "main" subroutine that calls child subroutines to
!      * calculate the condensation (uptake) coefficient
!        for all gas species and aerosol modes, and
!      * solve the condensation equations for nonvolatile species
!        and SOAs
!
!  The child subroutines that calculate the condensation (uptake) coefficients
!  and solve the condensation equations for nonvolatile species are also
!  included here.
!
!  The subroutines that solve the SOA condensation equations have been
!  (temporarily) moved to modal_aero_soaexch.F90.
!
!  Aging of the primary carbon mode as a result of condensation has been
!  moved to modal_aero_aging.F90
!
! History:
!   Original version by R. C. Easter
!   Separated from modal_aero_amicphys.F90 in EAMv1 and refactored by Hui Wan, 2020-2021
!-----------------------------------------------------------------------------------------
   use haero_precision,   only: wp
   use haero,             only: modal_aerosol_config_t, aerosol_species_t, &
                                prognostics_t, atmosphere_t, diagnostics_t, tendencies_t
   use haero_constants,   only: pi, pstd => pressure_stp, &
                                r_universal => r_gas, &
                                mw_air => molec_weight_dry_air, &  ! kg/mole
                                vol_molar_air => molec_diffusion_dry_air

   use mam_soaexch, only: mam_soaexch_1subarea

   implicit none
   private

   public :: init           ! init subroutine
   public :: run            ! main subroutine for condensation
   public :: finalize       ! end  subroutine for condensation

   public :: gas_aer_uptkrates_1box1gas                ! not ment to be called directly
   public :: mam_gasaerexch_1subarea_1gas_nonvolatile
   !-------------------------
   ! Module global variables
   !-------------------------

   ! These variables are adapted from MAM's modal_aero_microp_control module
   ! and initialized in the init subroutine. Here they are variables, not
   ! parameters, because their values depend on the particular model
   ! configuration. Accordingly, their values must be retained between calls.
   save

   !> Settable parameters:
   !> * aitken_adjust_factor (real) - Adjustment factor for Aitken number
   !>     concentration tendency
   !> * h2so4_uptake_opt (integer) - Controls treatment of H2SO4 uptake:
   !> * h2so4_condense_opt (integer) - Controls treatment of H2SO4 condensation:
   !>     1 = sequential   calc. of gas-chem prod, then condensation loss
   !>     2 = simultaneous calc. of gas-chem prod and condensation loss
   !> * pbl_adjust_factor (real) - adjustment factor for nucleation rate
   !>     corrected for the planetary boundary layer.
   !> * nuc_adjust_factor (real) - adjustment factor for nucleation rate with
   !>     binary/ternary nucleation.
   public :: set_integer_param, set_logical_param, set_real_param

   !> Type of condensation equation and numerical methods to apply to different gas species
   integer, parameter :: NA   = 0  ! no condensation, no calculation
   integer, parameter :: ANAL = 1  ! quasi-analytical solution for nonvolatile species
   integer, parameter :: IMPL = 2  ! implicit method with adaptive time steps for SOAG

   !> array containing flags for all gas species
   integer, dimension(:), allocatable :: eqn_and_numerics_category

   !> array containg flags indicating whether a specific gas can condense to a specific aerosol mode
   logical ,dimension(:,:), allocatable  ::  l_gas_condense_to_mode

   !> The index of the Aitken mode
   integer :: nait

   !> The index of the primary carbon mode
   integer :: npca

   !> Index for secondary-organic aerosol species.
   integer :: nsoa

   !> Index for primary-organic aerosol species.
   integer :: npoa

   !> Index of h2so4 aerosol within the Aitken mode
   integer :: iaer_h2so4

   !> The molecular weight of h2so4 aerosol as assumed by the host atm model
   real(wp) :: mw_h2so4

   !> The molecular weight of h2so4 aerosol as assumed by the host atm model
   real(wp) :: vol_molar_h2so4

   !> save off config%num_aerosol_modes for use in subroutines
   integer ::  max_mode

   !> Save off config%num_gases
   integer :: max_gas

   !> Save off maxval(config%num_mode_species)
   integer :: max_aer

   !> accomodation coefficient for h2so4 condensation
   real(wp), parameter :: accom_coef_h2so4 = 0.65_wp

   !> ratio of gas uptake coeff w.r.t. that of h2so4
   real(wp),parameter :: soag_h2so4_uptake_coeff_ratio = 0.81_wp  ! for SOAG
   real(wp),parameter ::  nh3_h2so4_uptake_coeff_ratio = 2.08_wp  ! for NH3

   !> number geometric_mean diameter
   real(wp), dimension(:), allocatable :: alnsg_aer

   !> array containing the ratios of all gas species
   real(wp), dimension(:), allocatable :: uptk_rate_factor

   !> for a limiter applied to NH4
   real(wp),parameter :: aer_nh4_so4_molar_ratio_max = 2._wp

   !> current gas mix ratios (mol/mol)
   real(wp), dimension(:), allocatable :: qgas_cur

   !> average gas mix ratios over the dt integration
   real(wp), dimension(:), allocatable :: qgas_avg

   !> qgas_netprod_otrproc = gas net production rate from other processes
   !>    such as gas-phase chemistry and emissions (mol/mol/s)
   !> this allows the condensation (gasaerexch) routine to apply production and condensation loss
   !>    together, which is more accurate numerically
   !> NOTE - must be >= zero, as numerical method can fail when it is negative
   !> NOTE - currently only the values for h2so4 and nh3 should be non-zero
   real(wp), dimension(:), allocatable :: qgas_netprod_otrproc


   ! gas to aerosol mass transfer rate (1/s)
   real(wp), dimension(:,:), allocatable :: uptkaer

   ! quadrature parameter related to uptake rate
   real(wp), parameter :: beta_inp = 0._wp     ! quadrature parameter (--)

   ! dimensions: (aerosol species, modes)
   logical, dimension(:,:), allocatable :: l_mode_can_contain_species

   ! dimension:  (aerosol species)
   logical, dimension(:), allocatable   :: l_mode_can_age

   integer :: ngas                             ! total # of gases handled by the parameterization
   integer :: ntot_amode                       ! total # of aerosol modes handled by the parameterization
   integer :: igas_h2so4
   integer :: igas_nh3
   integer :: iaer_so4
   integer :: iaer_pom
   ! integer :: igas_soag_bgn  Not sure what SOAG is and how these were set in MAMRefactor
   ! integer :: igas_soag_end  but I think for Haero these should be ignored

   integer, dimension(:), allocatable :: idx_gas_to_aer                ! dimension: (gas species)

   integer, dimension(:), allocatable ::  mode_aging_optaa
  !----------------------------------------------------------------------

contains

  subroutine init(config)
   type(modal_aerosol_config_t), intent(in) :: config

   ! local variables
   integer :: igas, iaer

   type(aerosol_species_t) h2so4

   max_mode= config%num_aerosol_modes
   max_gas = config%num_gases
   ngas    = max_gas
   max_aer = maxval(config%num_mode_species)
   nait = config%aerosol_mode_index("aitken");
   npca = config%aerosol_mode_index("primary_carbon");

   allocate(l_mode_can_contain_species(maxval(config%num_mode_species), max_mode))
   allocate(l_gas_condense_to_mode(config%num_gases, max_mode))
   allocate(l_mode_can_age(maxval(config%num_mode_species)))
   allocate(idx_gas_to_aer(config%num_gases))
   allocate(eqn_and_numerics_category(config%num_gases))
   allocate(uptk_rate_factor(config%num_gases))

   allocate(qgas_cur(config%num_gases))
   allocate(qgas_avg(config%num_gases))
   allocate(qgas_netprod_otrproc(config%num_gases))

   allocate(uptkaer(config%num_gases,max_mode)) ! gas to aerosol mass transfer rate (1/s)

   allocate(alnsg_aer(max_mode))
   allocate(mode_aging_optaa(max_mode))

   !------------------------------------------------------------------
   ! MAM currently assumes that the uptake rate of other gases
   ! are proportional to the uptake rate of sulfuric acid gas (H2SO4).
   ! Here the array uptk_rate_factor contains the uptake rate ratio
   ! w.r.t. that of H2SO4.
   !------------------------------------------------------------------
   ! Set default to 0, meaning no uptake
   uptk_rate_factor(:) = 0._wp

   igas_h2so4 = config%gas_index("H2SO4")
   igas_nh3   = config%gas_index("NH3")
   ! H2SO4 is the ref species, so the ratio is 1
   uptk_rate_factor(igas_h2so4) = 1._wp

   ! For SOAG. (igas_soag_bgn and igas_soag_end are the start- and end-indices)
   ! Remove use of igas_soag_bgn but keep comments in case need to be added back
   ! uptk_rate_factor(igas_soag_bgn:igas_soag_end) = soag_h2so4_uptake_coeff_ratio
   uptk_rate_factor(:) = soag_h2so4_uptake_coeff_ratio

   ! For NH3 (if implemented)
   if (igas_nh3 > 0) then
      uptk_rate_factor(igas_nh3) = nh3_h2so4_uptake_coeff_ratio
   end if

   ! Record the index for h2so4 aerosol within the Aitken mode and fetch some
   ! properties.
   iaer_h2so4 = config%aerosol_species_index(nait, "H2SO4")
   iaer_so4   = config%aerosol_species_index(nait, "SO4")
   iaer_pom   = config%aerosol_species_index(nait, "POM")
   if (iaer_h2so4 > 0) then
     h2so4 = config%aerosol_species(nait, iaer_h2so4)
     mw_h2so4 = h2so4%molecular_wt
   end if

   nsoa = config%aerosol_species_index(nait, "SOA")
   npoa = config%aerosol_species_index(nait, "POA")

   !-------------------------------------------------------------------
   ! MAM currently uses a splitting method to deal with gas-aerosol
   ! mass exchange. A quasi-analytical solution assuming timestep-wise
   ! constant uptake rate is applied to nonvolatile species, while
   ! an implicit time stepping method with adaptive step size
   ! is applied to gas-phase SOA species which are assumed semi-volatile.
   ! There are two different subroutines in this modules to deal
   ! with the two categories of cases. The array
   ! eqn_and_numerics_category flags which gas species belongs to
   ! which category.
   !-------------------------------------------------------------------
   eqn_and_numerics_category(:)                           = NA
   eqn_and_numerics_category(igas_h2so4)                  = ANAL
   ! Not sure what to do about igas_soag_bgn.  Seems not needed in HAERO
   !eqn_and_numerics_category(igas_soag_bgn:igas_soag_end) = IMPL

   if (igas_nh3 > 0 ) eqn_and_numerics_category(igas_nh3) = ANAL

   !-------------------------------------------------------------------
   ! Determine whether specific gases will condense to specific modes
   !-------------------------------------------------------------------
   l_gas_condense_to_mode(:,:) = .false.

   do igas = 1, ngas    !loop through all registered gas species
      if (eqn_and_numerics_category(igas) /= NA) then ! this gas species can condense

         iaer = idx_gas_to_aer(igas)   ! what aerosol species does the gas become when condensing?
         l_gas_condense_to_mode(igas,1:ntot_amode) = l_mode_can_contain_species(iaer,1:ntot_amode) &
                                                .or. l_mode_can_age(1:ntot_amode)
      end if
   end do ! igas

  end subroutine init

  subroutine finalize()
    implicit none

    ! Deallocate mode metadata.
    deallocate(l_mode_can_contain_species)
    deallocate(l_gas_condense_to_mode)
    deallocate(l_mode_can_age)
    deallocate(idx_gas_to_aer)
    deallocate(eqn_and_numerics_category)
    deallocate(uptk_rate_factor)

    deallocate(qgas_cur)
    deallocate(qgas_avg)
    deallocate(qgas_netprod_otrproc)

    deallocate(uptkaer)

    deallocate(alnsg_aer)
    deallocate(mode_aging_optaa)

  end subroutine

  subroutine set_integer_param(name, val)
    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: val

  end subroutine

  subroutine set_logical_param(name, val)
    implicit none

    character(len=*), intent(in) :: name
    logical, intent(in) :: val

    ! No logical parameters to set!
  end subroutine

  subroutine set_real_param(name, val)
    implicit none

    character(len=*), intent(in) :: name
    real(wp), intent(in) :: val

  end subroutine

  !----------------------------------------------------------------------
  subroutine run(t, dt, prognostics, atmosphere, diagnostics, tendencies)

     !   igas_h2so4, igas_nh3, idx_gas_to_aer, iaer_so4,       &! in
     !   l_calc_gas_uptake_coeff,                              &! in
     !   lund,                                                 &! in
     !   dt,                                                   &! in
     !   temp,              pmid,             aircon,          &! in
     !   ngas, n_mode, ntot_amode,                             &! in
     !   qgas_cur,          qgas_avg,                          &! inout, inout
     !   qgas_netprod_otrproc,                                 &! in
     !   qaer_cur,                                             &! inout
     !   qnum_cur,                                             &! inout
     !   dgn_awet,                                             &! in
     !   uptkaer,           uptkrate_h2so4                     )! inout, inout

      implicit none

      ! Arguments
      real(wp), value, intent(in)       :: t
      real(wp), value, intent(in)       :: dt    ! integration timestep (s)
      type(prognostics_t), intent(in)   :: prognostics
      type(atmosphere_t), intent(in)    :: atmosphere
      type(diagnostics_t), intent(in)   :: diagnostics
      type(tendencies_t), intent(inout) :: tendencies


      logical   :: l_calc_gas_uptake_coeff
      integer   :: lund        ! logical unit for diagnostic output
      !integer   :: ngas        ! number of potentially condensable gases
      integer   :: n_mode      ! current number of modes (including temporary) !Q: what are temporary modes?
      integer   :: nghq        ! number of ghq points
      integer   :: ntot_amode  ! total of modes (not including temporary?)

      real(wp)  :: temp        ! air temperature (K)
      real(wp)  :: pmid        ! air pressure at model levels (Pa)
      real(wp)  :: aircon      ! air molar concentration (kmol/m3)


     ! qgas/aer values are updated during the dt integration
      real(wp) :: qaer_cur(max_aer,max_mode)     ! current aerosol mass mix ratios (mol/mol)
      real(wp) :: qnum_cur(max_mode)             ! current aerosol number mix ratios (#/kmol)

      real(wp) :: dgn_awet(max_mode)             ! wet geo. mean dia. (m) of number distrib.

      real(wp) :: uptkrate_h2so4
                      ! h2so4(g) to aerosol mass transfer rate, summed over all modes (1/s)
                      ! this is needed by the nucleation routine (mam_newnuc_1subarea)

      integer  ::  lptr2_soa_a_amode(max_mode, nsoa)
      ! local variables

      integer :: iaer, igas
      integer :: n

      real(wp) :: mass_excess  !amount of nh4 mass exceeding aer_nh4_so4_molar_ratio_max*so4 on
                               !molar basis after condensation of nh3. This
                               !amount is then reduced to 0 by effectively evaporating nh4 back to nh3

      logical :: l_condense_to_mode(max_mode)

      real(wp) :: uptkaer_ref(1:max_mode)

      real(wp) :: r_pi = pi

      n_mode = max_mode
      nghq   = 2  !set number of ghq points for direct ghq
      ntot_amode = max_mode
      !=========================================================================
      ! Hui Wan's note from March 2021:
      ! I have not fully understood the mixed use of n_mode and ntot_amode
      ! in this module. Need to watch out for possible bugs.
      ! (It looks like ntot_amode is the # of modes tracked by the host model
      ! while n_mode also includes the # of temporary modes. Does this mean
      ! n_mode could be larger than ntot_amode? If we indeed run into a case
      ! where n_mode > ntot_amode, should the uptake rates be calculated
      ! for the temporary modes, too, or is it assumed that the gases will
      ! not condense to the temporary modes?)
      !=========================================================================
      if (n_mode/=ntot_amode) then
        print *, 'Value of n_mode not equal to ntot_amode. See comment ', &
                 'in mam_gasaerexch.F90::run(). n_mode:', n_mode, &
                 ' ntot_amode:', ntot_amode
        stop 1
      endif

      !=========================================================================
      ! Initialize the time-step mean gas concentration (explain why?)
      !=========================================================================
      qgas_avg(1:ngas) = 0.0_wp

      if (l_calc_gas_uptake_coeff) then
      !=========================================================================
      ! Calculate the reference uptake coefficient for all aerosol modes using
      ! properties of the H2SO4 gas
      !=========================================================================
        uptkaer_ref(:) = 0._wp           ! initialize with zero (-> no uptake)
        l_condense_to_mode(:) = .true.   ! do calcullation for ALL modes
        igas = igas_h2so4                ! use properties of the H2SO4 gas

        call gas_aer_uptkrates_1box1gas( l_condense_to_mode,                              &! in
                                         temp, pmid, pstd, mw_h2so4, 1000*mw_air,     &! in
                                         vol_molar_h2so4, vol_molar_air,              &! in
                                         accom_coef_h2so4, r_universal, r_pi, beta_inp, &! in
                                         ntot_amode, nghq,                                &! in
                                         dgn_awet(1:ntot_amode), alnsg_aer(1:ntot_amode), &! in
                                         uptkaer_ref(1:ntot_amode)                        )! inout

      !========================================================================================
      ! Assign uptake rate to each gas species and each mode using the ref. value uptkaer_ref
      ! calculated above and the uptake rate factor specified as constants at the
      ! beginning of the module
      !========================================================================================
        uptkaer(:,:) = 0._wp  !default is no uptake

        do igas = 1,ngas
           where( l_gas_condense_to_mode(igas,1:ntot_amode) )  &
           uptkaer(igas,1:ntot_amode) = uptkaer_ref(1:ntot_amode) * uptk_rate_factor(igas)
        end do

        ! total uptake rate (sum of all aerosol modes) for h2so4.
        ! Diagnosed for calling routine. Not used in this subroutne.
        uptkrate_h2so4 = sum( uptkaer(igas_h2so4,1:ntot_amode) )

      end if ! (l_calc_gas_uptake_coeff)

      !==================================================================================
      ! Solve condensation equation for non-volatile species
      !==================================================================================
      ! Using quasi-analytical solution (with no time sub-stepping)
      !---------------------------------------------------------------------------------
      do igas = 1,ngas

        if (eqn_and_numerics_category(igas).ne.ANAL) cycle

        iaer = idx_gas_to_aer(igas)

        call mam_gasaerexch_1subarea_1gas_nonvolatile( dt, qgas_netprod_otrproc(igas),  &! in, in
                                                       n_mode,  uptkaer(igas,:),        &! in, in
                                                       qgas_cur(igas),  qgas_avg(igas), &! inout, out
                                                       qaer_cur(iaer,:)                 )! inout
      end do

      !---------------------------------------------------------------------------------
      ! Clip condensation rate when nh3 is on the list of non-volatile gases:
      ! limit the condensation of nh3 so that nh4 does not exceed
      ! aer_nh4_so4_molar_ratio_max * so4 (molar basis).
      ! (Hui Wan's comment from Dec. 2020: chose to leave the following block here instead
      ! of moving it to the subroutine mam_gasaerexch_1subarea_nonvolatile_quasi_analytical,
      ! to make it a bit easier to see the assumed relationship between species.)
      !---------------------------------------------------------------------------------
      if ( igas_nh3 > 0 ) then

         igas = igas_nh3
         iaer = idx_gas_to_aer(igas)

         do n = 1, n_mode
            if (uptkaer(igas,n) <= 0.0_wp) cycle

            ! if nh4 exceeds aer_nh4_so4_molar_ratio_max*so4 (molar basis),
            ! put the excessive amount back to gas phase (nh3)

            mass_excess = qaer_cur(iaer,n) - aer_nh4_so4_molar_ratio_max*qaer_cur(iaer_so4,n)

            if (mass_excess > 0.0_wp) then
               qaer_cur(iaer,n) = qaer_cur(iaer,n) - mass_excess
               qgas_cur(igas)   = qgas_cur(igas)   + mass_excess
               qgas_avg(igas)   = qgas_avg(igas)   + mass_excess*0.5_wp
            end if
         end do

      end if ! igas_nh3 > 0

      !============================================
      ! Solve condensation equations for SOA
      !============================================
      call mam_soaexch_1subarea(                                    &
         max_gas, max_aer, iaer_pom,                                &
         nsoa, npoa, npca, pstd, r_universal,                 &
         lund,              dt,                                     &
         temp,              pmid,             aircon,               &
         n_mode,                                                    &
         ntot_amode,                                                &
         max_mode,                                                  &
         qgas_cur,          qgas_avg,                               &
         qaer_cur,                                                  &
         qnum_cur,                                                  &
         uptkaer, mode_aging_optaa, lptr2_soa_a_amode               )

  end subroutine run


  !--------------------------------------------------------------------------------
  function mean_molecular_speed( temp, rmw, r_universal, pi )

    implicit none

    real(wp) :: mean_molecular_speed    ! (m/s)
    real(wp),intent(in) :: temp         ! temperature (K)
    real(wp),intent(in) :: rmw          ! molec. weight (g/mol)
    real(wp),intent(in) :: r_universal  ! universal gas constant
    real(wp),intent(in) :: pi           ! pi

    mean_molecular_speed = sqrt(8._wp*r_universal*temp/(pi*rmw))

  end function mean_molecular_speed


  !--------------------------------------------------------------------------------
  function gas_diffusivity( T_in_K, p_in_atm, mw_gas, mw_air, vd_gas, vd_air )

    implicit none

    real(wp) :: gas_diffusivity       ! (m2/s)
    real(wp),intent(in) :: T_in_K     ! temperature (K)
    real(wp),intent(in) :: p_in_atm   ! pressure (atmospheres)
    real(wp),intent(in) :: mw_gas     ! molec. weight of the condensing gas (g/mol)
    real(wp),intent(in) :: mw_air     ! molec. weight of air (g/mol)
    real(wp),intent(in) :: vd_gas     ! molec. diffusion volume of the condensing gas
    real(wp),intent(in) :: vd_air     ! molec. diffusion volume of air

    real(wp),parameter :: onethird = 1._wp/3._wp

    gas_diffusivity = (1.0e-7_wp * T_in_K**1.75_wp        &
                    * sqrt(1._wp/mw_gas + 1._wp/mw_air))  &
                    / (p_in_atm * (vd_gas**onethird + vd_air**onethird)**2 )

  end function gas_diffusivity


  !--------------------------------------------------------------------------------
  subroutine gas_aer_uptkrates_1box1gas( l_condense_to_mode, &
                                         temp, pmid, pstd, mw_gas, mw_air, vol_molar_gas, vol_molar_air, &
                                         accom, r_universal, pi, beta_inp, &
                                         n_mode, nghq, dgncur_awet, lnsg, &
                                         uptkaer )
!---------------------------------------------------------------------------------
!   Computes   uptake rate parameter uptkaer(1:n_mode) = uptkrate(1:n_mode)
!
!                            /
!   where      uptkrate(i) = | gas_conden_rate(Dp) n_i(lnDp) dlnDp
!                            /
!
!   is the uptake rate for aerosol mode i with size distribution n_i(lnDp)
!   and number concentration of = 1 #/m3; aernum(i) is the actual number
!   mixing ratio of mode i in the unit of  #/kmol-air, and aircon is
!   the air concentration in the unit of kmol/m3.
!
!   gas_conden_rate(D_p) = 2 * pi * gasdiffus * D_p * F(Kn,ac), with
!           gasdiffus = gas diffusivity
!           F(Kn,ac) = Fuchs-Sutugin correction factor
!           Kn = Knudsen number (which is a function of Dp)
!           ac = accomodation coefficient (constant for each gas species)
!---------------------------------------------------------------------------------
!   using Gauss-Hermite quadrature of order nghq=2
!
!       D_p = particle diameter (cm)
!       x = ln(D_p)
!       dN/dx = log-normal particle number density distribution
!---------------------------------------------------------------------------------
      implicit none

      logical,  intent(in) :: l_condense_to_mode(n_mode) ! flags indicating whether gas can condense to aerosol
      real(wp), intent(in) :: temp             ! air temperature (K)
      real(wp), intent(in) :: pmid             ! air pressure at model levels (Pa)
      real(wp), intent(in) :: pstd             ! 101325 Pa
      real(wp), intent(in) :: mw_gas           ! molecular weight of gas
      real(wp), intent(in) :: mw_air           ! molecular weight of air
      real(wp), intent(in) :: vol_molar_gas    ! The molecular weight of aerosol as assumed by the host atm model
      real(wp), intent(in) :: vol_molar_air    ! The molecular weight of air as assumed by the host atm model

      real(wp), intent(in)  :: accom           ! accomodation coefficient (--)

      real(wp), intent(in)  :: r_universal     ! universal gas constant
      real(wp), intent(in)  :: pi              ! pi
      real(wp), intent(in)  :: beta_inp        ! quadrature parameter (--)

      integer,  intent(in)  :: n_mode              ! number of modes
      integer,  intent(in)  :: nghq                ! number of ghq points
      real(wp), intent(in)  :: dgncur_awet(n_mode) ! mode-median wet diameter of number distribution (m)
      real(wp), intent(in)  :: lnsg(n_mode)        ! ln( sigmag )  (--)

      real(wp), intent(out) :: uptkaer(n_mode)  ! gas-to-aerosol mass transfer rates (1/s)

      ! local
      integer :: iq, n

      real(wp), parameter :: tworootpi = 3.5449077018110320_wp
      real(wp), parameter :: root2     = 1.4142135623730950_wp

      real(wp), parameter :: one       = 1.0_wp
      real(wp), parameter :: two       = 2.0_wp

      real(wp) :: accomxp283, accomxp75
      real(wp) :: beta
      real(wp) :: const
      real(wp) :: D_p
      real(wp) :: knudsen
      real(wp) :: lndp, lndpgn
      real(wp) :: sumghq
      real(wp) :: tmpa

      real(wp) :: hh

      real(wp) :: p_in_atm
      real(wp) :: gasdiffus            ! gas diffusivity (m2/s)
      real(wp) :: gasfreepath          ! gas mean free path (m)

      real(wp) :: uptkrate(n_mode) ! gas-to-aerosol mass transfer rates (1/s)
                                   ! for number concentration = 1 #/m3

    ! Dick's old version
    ! integer, parameter :: nghq = 2
    ! real(wp), save :: xghq(nghq), wghq(nghq) ! quadrature abscissae and weights
    ! data xghq / 0.70710678_wp, -0.70710678_wp /
    ! data wghq / 0.88622693_wp,  0.88622693_wp /

      real(wp) :: xghq(nghq)
      real(wp) :: wghq(nghq)
      !------------------choose the integration points and weights ----------------------------------
      if ( nghq == 20) then
         xghq(1:nghq) = (/  &
         -5.3874808900112_wp, -4.6036824495507_wp, -3.9447640401156_wp, -3.3478545673832_wp, &
         -2.7888060584281_wp, -2.2549740020893_wp, -1.7385377121166_wp, -1.2340762153953_wp, &
         -0.73747372854539_wp, -0.2453407083009_wp, 0.2453407083009_wp, 0.73747372854539_wp, &
         1.2340762153953_wp, 1.7385377121166_wp, 2.2549740020893_wp, 2.7888060584281_wp, &
         3.3478545673832_wp, 3.9447640401156_wp, 4.6036824495507_wp, 5.3874808900112_wp /)
         wghq(1:nghq) = (/  &
         2.229393645534E-13_wp, 4.399340992273E-10_wp, 1.086069370769E-7_wp, 7.80255647853E-6_wp, &
         2.283386360164E-4_wp, 0.003243773342238_wp, 0.024810520887464_wp, 0.10901720602002_wp, &
         0.28667550536283_wp, 0.46224366960061_wp, 0.46224366960061_wp, 0.28667550536283_wp, &
         0.10901720602002_wp, 0.024810520887464_wp, 0.003243773342238_wp, 2.283386360164E-4_wp, &
         7.80255647853E-6_wp, 1.086069370769E-7_wp, 4.399340992273E-10_wp, 2.229393645534E-13_wp /)
      else if (nghq == 10) then
         xghq(1:nghq) = (/  &
       -3.436159118837737603327_wp, -2.532731674232789796409_wp, -1.756683649299881773451_wp, -1.036610829789513654178_wp, &
       -0.3429013272237046087892_wp, 0.3429013272237046087892_wp, 1.036610829789513654178_wp, 1.756683649299881773451_wp, &
       2.532731674232789796409_wp, 3.436159118837737603327_wp/)
         wghq(1:nghq) = (/  &
       7.64043285523262062916E-6_wp, 0.001343645746781232692202_wp,0.0338743944554810631362_wp,0.2401386110823146864165_wp, &
       0.6108626337353257987836_wp,0.6108626337353257987836_wp,0.2401386110823146864165_wp,0.03387439445548106313616_wp, &
       0.001343645746781232692202_wp,7.64043285523262062916E-6_wp/)
      else if (nghq == 4) then
         xghq(1:nghq) = (/ -1.6506801238858_wp, -0.52464762327529_wp, 0.52464762327529_wp, 1.6506801238858_wp /)
         wghq(1:nghq) = (/ 0.081312835447245_wp, 0.8049140900055_wp, 0.8049140900055_wp, 0.081312835447245_wp /)
      else if (nghq == 2) then
         xghq(1:nghq) = (/ -7.0710678118654746e-01_wp, 7.0710678118654746e-01_wp /)
         wghq(1:nghq) = (/  8.8622692545275794e-01_wp, 8.8622692545275794e-01_wp /)
      else
      write(*,*) 'Quadrature option is not available. Supported orders: 2, 4, 10 and 20 but not',nghq
      end if


      !----------------------------------------------------------------------------------------------

      p_in_atm  = pmid/pstd
      gasdiffus = gas_diffusivity( temp, p_in_atm, mw_gas, mw_air, vol_molar_gas, vol_molar_air)
      gasfreepath = 3.0_wp * gasdiffus / mean_molecular_speed( temp, mw_gas, r_universal, pi )

      accomxp283 = accom * 0.283_wp
      accomxp75  = accom * 0.75_wp

      ! outermost loop over all modes
      do n = 1, n_mode

         lndpgn = log( dgncur_awet(n) )   ! (m)

         ! beta = dln(uptake_rate)/dln(D_p)
         !      = 2.0 in free molecular regime, 1.0 in continuum regime
         ! if uptake_rate ~= a * (D_p**beta), then the 2 point quadrature is very accurate
         if (abs(beta_inp-1.5_wp) > 0.5_wp) then

           !D_p = dgncur_awet(n) * exp( 1.5_wp*(lnsg(n)**2) )
            D_p = dgncur_awet(n)


            knudsen = two*gasfreepath/D_p
            ! tmpa = dln(fuchs_sutugin)/d(knudsen)
            tmpa = one/(one+knudsen) - (two*knudsen + one + accomxp283) / &
                      ( knudsen*( knudsen + one + accomxp283 ) + accomxp75 )
            beta = one - knudsen*tmpa


            beta = max( one, min( two, beta ) )

         else
            beta = beta_inp
         end if

         const  = tworootpi * exp( beta*lndpgn + 0.5_wp*(beta*lnsg(n))**2 )

         ! sum over gauss-hermite quadrature points
         sumghq = 0.0_wp
         do iq = 1, nghq

            lndp = lndpgn + beta*lnsg(n)**2 + root2*lnsg(n)*xghq(iq)
            D_p = exp(lndp)

            hh = fuchs_sutugin( D_p, gasfreepath, accomxp283, accomxp75)

            sumghq = sumghq + wghq(iq)*   D_p*hh    /(D_p**beta)
         end do

         uptkrate(n) = const * gasdiffus * sumghq

      end do   ! n = 1, n_mode

      !-----------------------------------------------------------------------------------
      !  Unit of uptkrate is for number = 1 #/m3.
      !-----------------------------------------------------------------------------------
      where ( l_condense_to_mode(:) )
        uptkaer(:) = uptkrate(:)
      elsewhere
        uptkaer(:) = 0.0_wp  ! zero means no uptake
      end where

  contains

  !------
  function fuchs_sutugin( D_p, gasfreepath, accomxp283, accomxp75)

    real(wp) :: fuchs_sutugin

    real(wp),intent(in) :: D_p                    ! diameter (m) of a single particle
    real(wp),intent(in) :: gasfreepath            ! gas mean free path (m)
    real(wp),intent(in) :: accomxp283, accomxp75  ! accommodation coefficient * 0.283 or 0.75

    real(wp) :: knudsen

    knudsen = 2._wp * gasfreepath /D_p

    !fkn = ( 0.75*accomcoef*(1. + xkn) ) / &
    !      ( xkn**2 + xkn + 0.283*xkn*accomcoef + 0.75*accomcoef )

    fuchs_sutugin = ( accomxp75*(1._wp + knudsen) ) / &
                    ( knudsen*( knudsen + 1._wp + accomxp283 ) + accomxp75 )

  end function fuchs_sutugin

  end subroutine gas_aer_uptkrates_1box1gas

  !----------------------------------------------------------------------
  subroutine mam_gasaerexch_1subarea_1gas_nonvolatile( dt, qgas_netprod_otrproc, &
                                                       n_mode,  uptkaer,         &
                                                       qgas_cur,  qgas_avg,      &
                                                       qaer_cur )

      integer,  intent(in) :: n_mode    ! current number of modes (including temporary)
      real(wp), intent(in) :: dt        ! current integration timestep (s)

      real(wp), intent(in) :: qgas_netprod_otrproc
                ! qgas_netprod_otrproc = gas net production rate from other processes
                !    such as gas-phase chemistry and emissions (mol/mol/s)
                ! this allows the condensation (gasaerexch) routine to apply production and condensation loss
                !    together, which is more accurate numerically
                ! NOTE - must be >= zero, as numerical method can fail when it is negative
                ! NOTE - currently only the values for h2so4 and nh3 should be non-zero

      real(wp), intent(in)    :: uptkaer(1:max_mode)
      real(wp), intent(inout) :: qgas_cur
      real(wp), intent(out)   :: qgas_avg
      real(wp), intent(inout) :: qaer_cur(1:max_mode)

      !local

      integer ::  n

      real(wp) :: qgas_prv
      real(wp) :: qaer_prv(1:max_mode)

      real(wp) :: total_uptake_coeff_of_all_modes  !the uptake coefficient of a single gas species
                                                   !summed over all aerosol modes
      real(wp) :: tmp_kxt, tmp_kxt2, tmp_pxt

      real(wp) :: zqgas_init, zqgas_end, zqgas_avg, zqgas_equil
      real(wp) :: tmp_qdel_cond

         ! Save current values as the "previous" value for new time step (is this needed?)

         qgas_prv           = qgas_cur
         qaer_prv(1:n_mode) = qaer_cur(1:n_mode)

         ! misc. quantities that will occur repeatedly later equations

         total_uptake_coeff_of_all_modes = sum( uptkaer(1:n_mode) )

         tmp_kxt = total_uptake_coeff_of_all_modes*dt    ! uptake rate * dt
         tmp_pxt = qgas_netprod_otrproc*dt               ! prodc. rate * dt

         ! zqgas_init (was tmp_q1) = mix-rat at t=tcur
         zqgas_init = qgas_prv

         ! zqgas_end (was tmp_q3) = mix-rat at t=tcur+dt
         ! zqgas_avg (was tmp_q4) = avg mix-rat between t=tcur and t=tcur+dt

         if (tmp_kxt < 1.0e-20_wp) then
         ! Consider uptake to aerosols = 0.0, hence there is no change in aerosol mass;
         ! gas concentration is updated using the passed-in production rate
         ! which is considered step-wise constant

            qgas_cur = zqgas_init + tmp_pxt
            qgas_avg = zqgas_init + tmp_pxt*0.5_wp

         else
         ! tmp_kxt >= 1.0e-20_wp: there is non-negligible condensation;
         ! calculate amount of gas condensation and update gas concentration.

            if (tmp_kxt > 0.001_wp) then
            ! Stronger condesation, no worry about division by zero or singularity;
            ! calculate analytical solution assuming step-wise constant condensation coeff.

               zqgas_equil = tmp_pxt/tmp_kxt
               zqgas_end   = (zqgas_init - zqgas_equil)*exp(-tmp_kxt) + zqgas_equil
               zqgas_avg   = (zqgas_init - zqgas_equil)*(1.0_wp - exp(-tmp_kxt))/tmp_kxt + zqgas_equil

            else
            ! Weaker condensation, use Taylor expansion to avoid large error
            ! resulting from small denominator

               tmp_kxt2  = tmp_kxt*tmp_kxt
               zqgas_end = zqgas_init *(1.0_wp - tmp_kxt        + tmp_kxt2*0.5_wp) &
                         + tmp_pxt    *(1.0_wp - tmp_kxt*0.5_wp + tmp_kxt2/6.0_wp)
               zqgas_avg = zqgas_init *(1.0_wp - tmp_kxt*0.5_wp + tmp_kxt2/6.0_wp) &
                         + tmp_pxt    *(0.5_wp - tmp_kxt/6.0_wp + tmp_kxt2/24.0_wp)
            end if
            qgas_cur      = zqgas_end
            qgas_avg      = zqgas_avg
            tmp_qdel_cond = (zqgas_init + tmp_pxt) - zqgas_end   ! amount of condensed mass

            ! Distribute the condensed mass to different aerosol modes

            do n = 1, n_mode
               if (uptkaer(n) <= 0.0_wp) cycle
               qaer_cur(n) = qaer_prv(n) &
                           + tmp_qdel_cond*(uptkaer(n)/total_uptake_coeff_of_all_modes)
            end do

         end if

  end subroutine mam_gasaerexch_1subarea_1gas_nonvolatile
  !----------------------------------------------------------------------
end module mam_gasaerexch


