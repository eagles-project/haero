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
   use haero,             only: model_t, aerosol_species_t, &
                                prognostics_t, atmosphere_t, diagnostics_t, tendencies_t
   use haero_constants,   only: pi, pstd => pressure_stp, &
                                r_universal => r_gas, &
                                mw_air => molec_weight_dry_air, &  ! kg/mole 
                                vol_molar_air => molec_diffusion_dry_air

   implicit none
   private

   public :: init           ! init subroutine
   public :: run            ! main subroutine for condensation
   public :: finalize       ! end  subroutine for condensation

   public :: mam_soaexch_1subarea                      ! for unit testing only, 
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

   !> save off model%num_modes for use in subroutines
   integer ::  max_mode  

   !> Save off model%num_gases
   integer :: max_gas
  
   !> Save off maxval(model%num_mode_species)
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
   integer :: igas_soag_bgn
   integer :: igas_soag_end

   integer, dimension(:), allocatable :: idx_gas_to_aer                ! dimension: (gas species)

   integer, dimension(:), allocatable ::  mode_aging_optaa
  !----------------------------------------------------------------------

contains

  subroutine init( model)
   type(model_t), intent(in) :: model

   ! local variables
   integer :: igas, iaer

   type(aerosol_species_t) h2so4

   max_mode= model%num_modes
   max_gas = model%num_gases
   ngas    = max_gas
   max_aer = maxval(model%num_mode_species)

   allocate(l_mode_can_contain_species(maxval(model%num_mode_species), model%num_modes))
   allocate(l_gas_condense_to_mode(model%num_gases, model%num_modes))
   allocate(l_mode_can_age(maxval(model%num_mode_species)))
   allocate(idx_gas_to_aer(model%num_gases))
   allocate(eqn_and_numerics_category(model%num_gases))
   allocate(uptk_rate_factor(model%num_gases))

   allocate(qgas_cur(model%num_gases))
   allocate(qgas_avg(model%num_gases))
   allocate(qgas_netprod_otrproc(model%num_gases))

   allocate(uptkaer(model%num_gases,model%num_modes))      ! gas to aerosol mass transfer rate (1/s)

   allocate(alnsg_aer(model%num_modes))
   allocate(mode_aging_optaa(model%num_modes))

   !------------------------------------------------------------------
   ! MAM currently assumes that the uptake rate of other gases 
   ! are proportional to the uptake rate of sulfuric acid gas (H2SO4).
   ! Here the array uptk_rate_factor contains the uptake rate ratio 
   ! w.r.t. that of H2SO4.
   !------------------------------------------------------------------
   ! Set default to 0, meaning no uptake
   uptk_rate_factor(:) = 0._wp

   ! H2SO4 is the ref species, so the ratio is 1
   uptk_rate_factor(igas_h2so4) = 1._wp

   ! For SOAG. (igas_soag_bgn and igas_soag_end are the start- and end-indices)
   uptk_rate_factor(igas_soag_bgn:igas_soag_end) = soag_h2so4_uptake_coeff_ratio

   ! For NH3 (if implemented)
   if (igas_nh3 > 0) then
      uptk_rate_factor(igas_nh3) = nh3_h2so4_uptake_coeff_ratio
   end if

   ! Record the index for h2so4 aerosol within the Aitken mode and fetch some
   ! properties.
   iaer_h2so4 = model%aerosol_index(nait, "H2SO4")
   if (iaer_h2so4 > 0) then
     h2so4 = model%aero_species(nait, iaer_h2so4)
     mw_h2so4 = h2so4%molecular_wt
   end if

   nsoa = model%aerosol_index(nait, "SOA")
   npoa = model%aerosol_index(nait, "POA")

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
   eqn_and_numerics_category(igas_soag_bgn:igas_soag_end) = IMPL

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

  subroutine finalize(model)
    implicit none

    ! Arguments
    type(model_t), intent(in) :: model

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
  subroutine run(model, t, dt, prognostics, atmosphere, diagnostics, tendencies)

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
      type(model_t), intent(in)         :: model
      real(wp), value, intent(in)       :: t
      real(wp), value, intent(in)       :: dt    ! integration timestep (s)
      type(prognostics_t), intent(in)   :: prognostics
      type(atmosphere_t), intent(in)    :: atmosphere
      type(diagnostics_t), intent(in)   :: diagnostics
      type(tendencies_t), intent(inout) :: tendencies

      integer   :: igas_h2so4, igas_nh3, iaer_so4

      logical   :: l_calc_gas_uptake_coeff 
      integer   :: lund        ! logical unit for diagnostic output
      !integer   :: ngas        ! number of potentially condensable gases 
      integer   :: n_mode      ! current number of modes (including temporary) !Q: what are temporary modes?
      integer   :: ntot_amode  ! total of modes (not including temporary?) 

      real(wp)  :: temp        ! air temperature (K)
      real(wp)  :: pmid        ! air pressure at model levels (Pa)
      real(wp)  :: aircon      ! air molar concentration (kmol/m3)


     ! qgas/aer values are updated during the dt integration
      real(wp) :: qaer_cur(max_aer,model%num_modes)     ! current aerosol mass mix ratios (mol/mol)
      real(wp) :: qnum_cur(model%num_modes)             ! current aerosol number mix ratios (#/kmol)

      real(wp) :: dgn_awet(model%num_modes)            ! wet geo. mean dia. (m) of number distrib.

      real(wp) :: uptkrate_h2so4
                      ! h2so4(g) to aerosol mass transfer rate, summed over all modes (1/s)
                      ! this is needed by the nucleation routine (mam_newnuc_1subarea)

      integer  ::  lptr2_soa_a_amode(model%num_modes, nsoa)
      ! local variables

      integer :: iaer, igas
      integer :: n

      real(wp) :: mass_excess  !amount of nh4 mass exceeding aer_nh4_so4_molar_ratio_max*so4 on 
                               !molar basis after condensation of nh3. This 
                               !amount is then reduced to 0 by effectively evaporating nh4 back to nh3

      logical :: l_condense_to_mode(model%num_modes)

      real(wp) :: uptkaer_ref(1:model%num_modes) 

      real(wp) :: r_pi = pi

      n_mode = model%num_modes
      ntot_amode = model%num_modes
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
                                         ntot_amode,                                      &! in
                                         dgn_awet(1:ntot_amode), alnsg_aer(1:ntot_amode), &! in
                                         qnum_cur(1:ntot_amode), aircon,                  &! in
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
                                         n_mode, dgncur_awet, lnsg, &
                                         aernum, aircon, &
                                         uptkaer )
!---------------------------------------------------------------------------------
!   Computes   uptkaer(1:n_mode) = uptkrate(1:n_mode) * aernum(1:n_mode) * aircon
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
      real(wp), intent(in)  :: dgncur_awet(n_mode) ! mode-median wet diameter of number distribution (m)
      real(wp), intent(in)  :: lnsg(n_mode)        ! ln( sigmag )  (--)

      real(wp), intent(in)  :: aernum(n_mode)   ! aerosol number mixing ratio 
      real(wp), intent(in)  :: aircon           ! air molar concentration (kmol/m3)

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

      integer, parameter :: nghq = 2
      real(wp),parameter :: xghq(nghq) = (/ -7.0710678118654746e-01_wp, 7.0710678118654746e-01_wp /)
      real(wp),parameter :: wghq(nghq) = (/  8.8622692545275794e-01_wp, 8.8622692545275794e-01_wp /)
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
      ! Unit conversion: uptkrate is for number = 1 #/m3, so mult. by number conc. (#/m3)
      !-----------------------------------------------------------------------------------
      where ( l_condense_to_mode(:) ) 
        uptkaer(:) = uptkrate(:) * (aernum(:) * aircon)
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


!----------------------------------------------------------------------
  subroutine mam_soaexch_1subarea(                              &
         lund,                 &
         dt,                   &
         temp,                 &
         pmid,                 &
         aircon,               &
         n_mode,               &
         ntot_amode,           &
         max_mode,             &
         qgas_cur,             &
         qgas_avg,             &
         qaer_cur,             &
         qnum_cur,             &
         uptkaer,              &
         mode_aging_optaa,     &
         lptr2_soa_a_amode     )
!
! calculate soa condensation/evaporation over time dt
!

      implicit none

! arguments
      integer,  intent(in) :: lund                  ! logical unit for diagnostic output
      integer,  intent(in) :: n_mode                ! current number of modes (including temporary)
      integer,  intent(in) :: ntot_amode 
      integer,  intent(in) :: max_mode

      real(wp), intent(in) :: dt        ! current integration timestep (s)
      real(wp), intent(in) :: temp             ! temperature (K)
      real(wp), intent(in) :: pmid             ! pressure at model levels (Pa)
      real(wp), intent(in) :: aircon           ! air molar concentration (kmol/m3)

      real(wp), intent(inout), dimension( 1:max_gas ) :: qgas_cur, qgas_avg
      real(wp), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: qaer_cur
      real(wp), intent(inout), dimension( 1:max_mode ) :: qnum_cur
      real(wp), intent(in   ), dimension( 1:max_gas, 1:max_mode ) :: uptkaer
      integer,  intent(in   ), dimension( 1:max_mode) ::  mode_aging_optaa
      integer,  intent(in   ), dimension( 1:max_mode, 1:nsoa) ::  lptr2_soa_a_amode

! local
      integer :: ntot_poaspec 
      integer :: ntot_soaspec 

      integer :: iaer, igas
      integer :: ll
      integer :: n, niter, niter_max
      integer :: ntot_soamode

      logical, parameter :: flag_pcarbon_opoa_frac_zero   = .true.

      logical :: skip_soamode(max_mode)   ! true if this mode does not have soa

      real(wp), dimension( 1:max_gas ) :: &
         qgas_prv

      real(wp), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_prv

      real(wp) :: uptkaer_soag_tmp(nsoa,max_mode)

      real(wp), parameter :: a_min1 = 1.0e-20_wp
      real(wp), parameter :: g_min1 = 1.0e-20_wp
#if ( defined( MAM_STANDALONE ) )
      real(wp) :: alpha_astem
#else
      real(wp), parameter :: alpha_astem = 0.05_wp ! parameter used in calc of time step
#endif
      real(wp), parameter :: dtsub_fixed = -1.0_wp ! fixed sub-step for time integration (s)
!     real(wp), parameter :: dtsub_fixed = 10.0    ! fixed sub-step for time integration (s)
      real(wp), parameter :: rgas = 8.3144_wp      ! gas constant in J/K/mol

      real(wp) :: a_ooa_sum_tmp(max_mode)          ! total ooa (=soa+opoa) in a mode
      real(wp) :: a_opoa(max_mode)                 ! oxidized-poa aerosol mixrat (mol/mol at actual mw)
      real(wp) :: a_soa(nsoa,max_mode)     ! soa aerosol mixrat (mol/mol at actual mw)
      real(wp) :: a_soa_tmp(nsoa,max_mode) ! temporary soa aerosol mixrat (mol/mol)
      real(wp) :: beta(nsoa,max_mode)      ! dtcur*xferrate
      real(wp) :: delh_vap_soa(nsoa)       ! delh_vap_soa = heat of vaporization for gas soa (J/mol)
      real(wp) :: del_g_soa_tmp(nsoa)
      real(wp) :: dtcur                            ! current time step (s)
      real(wp) :: dtfull                           ! full time step (s)
      real(wp) :: dtmax                            ! = (dtfull-tcur)
      real(wp) :: dtsum_qgas_avg
      real(wp) :: g0_soa(nsoa)             ! ambient soa gas equilib mixrat (mol/mol at actual mw)
      real(wp) :: g_soa(nsoa)              ! soa gas mixrat (mol/mol at actual mw)
      real(wp) :: g_star(nsoa,max_mode)    ! soa gas mixrat that is in equilib
                                                   ! with each aerosol mode (mol/mol)
!     real(wp) :: mw_poa(npoa)             ! actual molec wght of poa
!     real(wp) :: mw_soa(nsoa)             ! actual molec wght of soa
      real(wp) :: opoa_frac(npoa,max_mode) ! fraction of poa that is opoa
      real(wp) :: phi(nsoa,max_mode)       ! "relative driving force"
      real(wp) :: p0_soa(nsoa)             ! soa gas equilib vapor presssure (atm)
      real(wp) :: p0_soa_298(nsoa)         ! p0_soa_298 = soa gas equilib vapor presssure (atm) at 298 k
      real(wp) :: sat(nsoa,max_mode)       ! sat(m,ll) = g0_soa(ll)/a_ooa_sum_tmp(m) = g_star(m,ll)/a_soa(m,ll)
                                                   !    used by the numerical integration scheme -- it is not a saturation rato!
      real(wp) :: tcur                             ! current integration time (from 0 s)

      real(wp) :: tmpa, tmpb, tmpc

      real(wp) :: tot_soa(nsoa)            ! g_soa + sum( a_soa(:) )

      ntot_poaspec = npoa
      ntot_soaspec = nsoa
! calc ntot_soamode = "last" mode on which soa is allowed to condense
      ntot_soamode = 0
      do n = 1, ntot_amode
!        if (n == nufi) cycle
         if (mode_aging_optaa(n) > 0) ntot_soamode = n
         if (lptr2_soa_a_amode(n,1) > 0) ntot_soamode = n
      end do
#ifdef l_debug_print
#if ( defined( MAM_STANDALONE ) )
      if ( i*k == top_lev .and. ldiagd1 ) write(lund,'(/a,5i5)') &
         'ntot_amode, ntot_amode_extd, n_mode, ntot_soamode', &
          ntot_amode, ntot_amode_extd, n_mode, ntot_soamode
#endif
#endif

      opoa_frac = 0.1_wp
! for primary carbon mode, set opoa_frac=0 for consistency with older code
! (this could be changed)
!     if ( flag_pcarbon_opoa_frac_zero ) then
!        if (npca > 0) opoa_frac(:,npca) = 0.0_wp
!     end if

      delh_vap_soa = 156.0e3_wp
!     delh_vap_soa =  30.0e3  ! 11-jun-2012
      p0_soa_298 = 1.0e-10_wp

! calc ambient equilibrium soa gas
      do ll = 1, ntot_soaspec
         p0_soa(ll) = p0_soa_298(ll) * &
! JS changes on 08-28-2019
!                  exp( -(delh_vap_soa(ll)/rgas)*((1._wp/temp)-(1._wp/298._wp)) )
                  exp( -(delh_vap_soa(ll)/(r_universal/1.e3))*((1._wp/temp)-(1._wp/298._wp)) )
         g0_soa(ll) = pstd*p0_soa(ll)/pmid
      end do

      niter_max = 1000
      niter = 0
      dtfull = dt
      tcur = 0._wp
      dtcur = 0._wp
      phi(:,:) = 0._wp
      g_star(:,:) = 0._wp
      g_soa(:) = 0._wp
      a_opoa(:) = 0._wp
      a_soa(:,:) = 0._wp

#if ( defined( MAM_STANDALONE ) )
      alpha_astem = alpha_astem_soa_boxtest
      niter_max = niter_max_soa_boxtest
#endif

!
! main integration loop -- does multiple substeps to reach dtfull
!
      qgas_avg(1:nsoa) = 0.0_wp
      dtsum_qgas_avg = 0.0_wp

time_loop: &
      do while (tcur < dtfull-1.0e-3_wp )

      niter = niter + 1
      if (niter > niter_max) exit


! set qxxx_prv to be current value
      qgas_prv(1:nsoa) = qgas_cur(1:nsoa)
      qaer_prv = qaer_cur
!     qaer_num = qnum_cur


! determine which modes have non-zero transfer rates
!    and are involved in the soa gas-aerosol transfer
! for diameter = 1 nm and number = 1 #/cm3, xferrate ~= 1e-9 s-1
      do n = 1, ntot_soamode
         skip_soamode(n) = .true.
         do ll = 1, ntot_soaspec
            if (uptkaer(ll,n) > 1.0e-15_wp) then
               uptkaer_soag_tmp(ll,n) = uptkaer(ll,n)
               skip_soamode(n) = .false.
            else
               uptkaer_soag_tmp(ll,n) = 0.0_wp
            end if
         end do
      end do

! load incoming soag and soaa into temporary arrays
! force things to be non-negative
! calc tot_soa(ll)
! calc a_opoa (always slightly >0)
!
! *** questions ***
! > why not use qgas and qaer instead of g_soa and a_soa
! > why not calc the following on every substep because
!      nuc and coag may change things:
!      skip)soamode, uptkaer_soag_tmp, tot_soa, a_opoa
! > include gasprod for soa ??
! > create qxxx_bgn = qxxx_cur at the very beginning (is it needed)
!
      do ll = 1, ntot_soaspec
         g_soa(ll) = max( qgas_prv(ll), 0.0_wp )
         tot_soa(ll) = g_soa(ll)
         do n = 1, ntot_soamode
            if ( skip_soamode(n) ) cycle
            a_soa(ll,n) = max( qaer_prv(ll,n), 0.0_wp )
            tot_soa(ll) = tot_soa(ll) + a_soa(ll,n)
         end do
      end do

      do n = 1, ntot_soamode
         if ( skip_soamode(n) ) cycle
         a_opoa(n) = 0.0_wp
         do ll = 1, ntot_poaspec
            a_opoa(n) = a_opoa(n) + opoa_frac(ll,n) * max( qaer_prv(ll,n), 0.0_wp )
         end do
      end do


! determine time step
      tmpa = 0.0_wp  ! time integration parameter for all soa species
      do n = 1, ntot_soamode
         if ( skip_soamode(n) ) cycle
         a_ooa_sum_tmp(n) = a_opoa(n) + sum( a_soa(1:ntot_soaspec,n) )
      end do
      do ll = 1, ntot_soaspec
         tmpb = 0.0  ! time integration parameter for a single soa species
         do n = 1, ntot_soamode
            if ( skip_soamode(n) ) cycle
            sat(ll,n) = g0_soa(ll)/max( a_ooa_sum_tmp(n), a_min1 )
            g_star(ll,n) = sat(ll,n)*a_soa(ll,n)
            phi(ll,n) = (g_soa(ll) - g_star(ll,n))/max( g_soa(ll), g_star(ll,n), g_min1 )
            tmpb = tmpb + uptkaer_soag_tmp(ll,n)*abs(phi(ll,n))
         end do
         tmpa = max( tmpa, tmpb )
      end do

      if (dtsub_fixed > 0.0_wp) then
         dtcur = dtsub_fixed
         tcur = tcur + dtcur
      else
         dtmax = dtfull-tcur
         if (dtmax*tmpa <= alpha_astem) then
! here alpha_astem/tmpa >= dtmax, so this is final substep
            dtcur = dtmax
            tcur = dtfull
         else
            dtcur = alpha_astem/tmpa
            tcur = tcur + dtcur
         end if
      end if


! step 1 - for modes where soa is condensing, estimate "new" a_soa(ll,n)
!    using an explicit calculation with "old" g_soa
!    and g_star(ll,n) calculated using "old" a_soa(ll,n)
! do this to get better estimate of "new" a_soa(ll,n) and sat(ll,n)
      do n = 1, ntot_soamode
         if ( skip_soamode(n) ) cycle
         do ll = 1, ntot_soaspec
            ! first ll loop calcs a_soa_tmp(ll,n) & a_ooa_sum_tmp
            a_soa_tmp(ll,n) = a_soa(ll,n)
            beta(ll,n) = dtcur*uptkaer_soag_tmp(ll,n)
            del_g_soa_tmp(ll) = g_soa(ll) - g_star(ll,n)
            if (del_g_soa_tmp(ll) > 0.0_wp) then
               a_soa_tmp(ll,n) = a_soa(ll,n) + beta(ll,n)*del_g_soa_tmp(ll)
            end if
         end do
         a_ooa_sum_tmp(n) = a_opoa(n) + sum( a_soa_tmp(1:ntot_soaspec,n) )
         do ll = 1, ntot_soaspec
            ! second ll loop calcs sat & g_star
            if (del_g_soa_tmp(ll) > 0.0_wp) then
               sat(ll,n) = g0_soa(ll)/max( a_ooa_sum_tmp(n), a_min1 )
               g_star(ll,n) = sat(ll,n)*a_soa_tmp(ll,n)   ! this just needed for diagnostics
            end if
         end do
      end do


! step 2 - implicit in g_soa and semi-implicit in a_soa,
!    with g_star(ll,n) calculated semi-implicitly
      do ll = 1, ntot_soaspec
         tmpa = 0.0_wp
         tmpb = 0.0_wp
         do n = 1, ntot_soamode
            if ( skip_soamode(n) ) cycle
            tmpa = tmpa + a_soa(ll,n)/(1.0_wp + beta(ll,n)*sat(ll,n))
            tmpb = tmpb + beta(ll,n)/(1.0_wp + beta(ll,n)*sat(ll,n))
         end do

         g_soa(ll) = (tot_soa(ll) - tmpa)/(1.0_wp + tmpb)
         g_soa(ll) = max( 0.0_wp, g_soa(ll) )
         do n = 1, ntot_soamode
            if ( skip_soamode(n) ) cycle
            a_soa(ll,n) = (a_soa(ll,n) + beta(ll,n)*g_soa(ll))/   &
                       (1.0_wp + beta(ll,n)*sat(ll,n))
         end do
      end do


! update mix ratios for soa species
      do igas = 1, nsoa
         iaer = igas
         qgas_cur(igas) = g_soa(igas)
         tmpc = qgas_cur(igas) - qgas_prv(igas)
         qgas_avg(igas) = qgas_avg(igas) + dtcur*(qgas_prv(igas) + 0.5_wp*tmpc)
         do n = 1, ntot_soamode
            qaer_cur(iaer,n) = a_soa(iaer,n)
            tmpc = qaer_cur(iaer,n) - qaer_prv(iaer,n)
         end do
      end do


      dtsum_qgas_avg = dtsum_qgas_avg + dtcur

      end do time_loop

! convert qgas_avg from sum_over[ qgas*dt_cut ] to an average
      do igas = 1, nsoa
         qgas_avg(igas) = max( 0.0_wp, qgas_avg(igas)/dtsum_qgas_avg )
      end do

      end subroutine mam_soaexch_1subarea

end module mam_gasaerexch


