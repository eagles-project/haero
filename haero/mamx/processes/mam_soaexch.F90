module mam_soaexch

   use haero_precision,   only: wp

   implicit none
   private

   public :: mam_soaexch_1subarea

contains

!----------------------------------------------------------------------
  subroutine mam_soaexch_1subarea(                              &
         max_gas,              &
         max_aer,              &
         iaer_pom,             &
         nsoa,                 &
         npoa,                 &
         npca,                 &
         pstd,                 &
         r_universal,          &
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
      integer,  intent(in) :: max_gas
      integer,  intent(in) :: max_aer
      integer,  intent(in) :: iaer_pom
      integer,  intent(in) :: nsoa
      integer,  intent(in) :: npoa
      integer,  intent(in) :: npca
      real(wp), intent(in) :: pstd
      real(wp), intent(in) :: r_universal
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

      integer :: ll
      integer :: n, niter_max
      integer :: ntot_soamode

      logical, parameter :: flag_pcarbon_opoa_frac_zero   = .true.

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

      real(wp) :: delh_vap_soa(nsoa)       ! delh_vap_soa = heat of vaporization for gas soa (J/mol)
      real(wp) :: dtfull                           ! full time step (s)
      real(wp) :: g0_soa(nsoa)             ! ambient soa gas equilib mixrat (mol/mol at actual mw)
!     real(wp) :: mw_poa(npoa)             ! actual molec wght of poa
!     real(wp) :: mw_soa(nsoa)             ! actual molec wght of soa
      real(wp) :: opoa_frac(npoa,max_mode) ! fraction of poa that is opoa
      real(wp) :: p0_soa(nsoa)             ! soa gas equilib vapor presssure (atm)
      real(wp) :: p0_soa_298(nsoa)         ! p0_soa_298 = soa gas equilib vapor presssure (atm) at 298 k

      ntot_poaspec = npoa
      ntot_soaspec = nsoa
! calc ntot_soamode = "last" mode on which soa is allowed to condense
      ntot_soamode = 0
      do n = 1, ntot_amode
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
      if ( flag_pcarbon_opoa_frac_zero ) then
         if (npca > 0) opoa_frac(:,npca) = 0.0_wp
      end if

      delh_vap_soa = 156.0e3_wp
!     delh_vap_soa =  30.0e3  ! 11-jun-2012
      p0_soa_298 = 1.0e-10_wp

! calc ambient equilibrium soa gas
      do ll = 1, ntot_soaspec
         p0_soa(ll) = p0_soa_298(ll) * &
! JS changes on 08-28-2019
!                  exp( -(delh_vap_soa(ll)/rgas)*((1._wp/temp)-(1._wp/298._wp)) )
                  exp( -(delh_vap_soa(ll)/(r_universal/1.e3_wp))*((1._wp/temp)-(1._wp/298._wp)) )
         g0_soa(ll) = pstd*p0_soa(ll)/pmid
      end do

      niter_max = 1000
      dtfull = dt

#if ( defined( MAM_STANDALONE ) )
      alpha_astem = alpha_astem_soa_boxtest
      niter_max = niter_max_soa_boxtest
#endif

!
! main integration loop -- does multiple substeps to reach dtfull
!
      qgas_avg(1:nsoa) = 0.0_wp

! calculate internal substep tendency
      call mam_soaexch_tend(max_gas, max_aer, max_mode, iaer_pom, nsoa, npoa, &
                            dtfull, niter_max, alpha_astem,    &
                            ntot_soamode, uptkaer,             &
                            g0_soa, opoa_frac,qgas_cur, qaer_cur, qgas_avg)

   end subroutine mam_soaexch_1subarea


subroutine mam_soaexch_tend( &
         max_gas,              &
         max_aer,              &
         max_mode,             &
         iaer_pom,             &
         nsoa,                 &
         npoa,                 &
         dtfull,               &
         niter_max,            &
         alpha_astem,          &
         ntot_soamode,         &
         uptkaer,              &
         g0_soa,               &
         opoa_frac,            &
         qgas_cur,             &
         qaer_cur,             &
         qgas_avg)

      implicit none

      integer,  intent(in) :: max_gas
      integer,  intent(in) :: max_aer
      integer,  intent(in) :: max_mode
      integer,  intent(in) :: iaer_pom
      integer,  intent(in) :: nsoa
      integer,  intent(in) :: npoa
      real(wp), intent(in) :: dtfull                ! host model dt (s)
      integer, intent(in)  :: niter_max             ! max interations
      real(wp), intent(in) :: alpha_astem           ! error control parameter for sub-timesteps
      integer, intent(in)  :: ntot_soamode          !"last" mode on which soa is allowed to condense

      real(wp), intent(inout), dimension( 1:max_gas ) :: qgas_cur, qgas_avg
      real(wp), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: qaer_cur
      real(wp), intent(in   ), dimension( 1:max_gas, 1:max_mode ) :: uptkaer
      real(wp), intent(in   ), dimension( 1:npoa, 1:max_mode ) :: opoa_frac ! fraction of poa that is opoa
      real(wp), intent(in   ), dimension( 1:nsoa) :: g0_soa    ! ambient soa gas equilib mixrat (mol/mol at actual mw)

      logical :: skip_soamode(max_mode)  ! true if this mode does not have soa
      real(wp) :: uptkaer_soag_tmp(nsoa,max_mode)
      real(wp), dimension( 1:max_gas ) :: &
         qgas_prv

      real(wp), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_prv
      integer :: ntot_poaspec
      integer :: ntot_soaspec
      real(wp), parameter :: a_min1 = 1.0e-20_wp
      real(wp), parameter :: g_min1 = 1.0e-20_wp
      real(wp), parameter :: dtsub_fixed = -1.0_wp ! fixed sub-step for time integration (s)

      real(wp) :: a_ooa_sum_tmp(max_mode)          ! total ooa (=soa+opoa) in a mode
      real(wp) :: a_opoa(max_mode)                 ! oxidized-poa aerosol mixrat (mol/mol at actual mw)
      real(wp) :: a_soa(nsoa,max_mode)             ! soa aerosol mixrat (mol/mol at actual mw)
      real(wp) :: a_soa_tmp(nsoa,max_mode) ! temporary soa aerosol mixrat (mol/mol)
      real(wp) :: beta(nsoa,max_mode)      ! dtcur*xferrate
      real(wp) :: del_g_soa_tmp(nsoa)
      real(wp) :: g_soa(nsoa)              ! soa gas mixrat (mol/mol at actual mw)
      real(wp) :: g_star(nsoa,max_mode), g_star_test(nsoa,max_mode)    ! soa gas mixrat that is in equilib
                                                   ! with each aerosol mode (mol/mol)
      real(wp) :: sat(nsoa,max_mode)       ! sat(m,ll) = g0_soa(ll)/a_ooa_sum_tmp(m) = g_star(m,ll)/a_soa(m,ll)
                                                   !    used by the numerical integration scheme -- it is not a saturation rato!
      real(wp) :: dtsum_qgas_avg
      real(wp) :: tmpa, tmpb, tmpc
      real(wp) :: dtcur, tcur           ! current timestep
      integer :: iaer, igas
      integer :: ll
      integer :: n, niter
      real(wp) :: phi(nsoa,max_mode)       ! "relative driving force"
      real(wp) :: tot_soa(nsoa)            ! g_soa + sum( a_soa(:) )

      ntot_poaspec = npoa
      ntot_soaspec = nsoa
      niter = 0
      tcur = 0._wp
      dtcur = 0._wp
      g_star(:,:) = 0._wp
      g_star_test(:,:) = 0._wp
      g_soa(:) = 0._wp
      a_opoa(:) = 0._wp
      a_soa(:,:) = 0._wp
      phi(:,:) = 0._wp
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
            a_opoa(n) = a_opoa(n) + opoa_frac(ll,n) * max( qaer_prv(iaer_pom+ll-1,n), 0.0_wp )
         end do
      end do

       call mam_soaexch_substep(max_mode, nsoa, npoa,              &
                                ntot_soamode, dtsub_fixed, dtfull, &
                                skip_soamode, uptkaer_soag_tmp,    &
                                a_soa, a_opoa, g_soa,              &
                                a_ooa_sum_tmp, sat ,               &
                                g0_soa, alpha_astem,               &
                                g_star, dtcur, tcur)

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
   end subroutine mam_soaexch_tend

 subroutine mam_soaexch_substep( &
         max_mode,             &
         nsoa,                 &
         npoa,                 &
         ntot_soamode,         &
         dtsub_fixed,          &
         dtfull,               &
         skip_soamode,         &
         uptkaer_soag_tmp,     &
         a_soa,                &
         a_opoa,               &
         g_soa,                &
         a_ooa_sum_tmp,        &
         sat ,                 &
         g0_soa,               &
         alpha_astem,          &
         g_star,               &
         dtcur,                &
         tcur)

      implicit none

      integer,  intent(in) :: max_mode
      integer,  intent(in) :: nsoa
      integer,  intent(in) :: npoa
      real(wp), intent(in) :: dtfull          ! host model dt (s)
      real(wp), intent(in) :: dtsub_fixed     ! fixed integration timestep (s)
      real(wp), intent(in) :: alpha_astem     ! error control parameter for sub-timesteps
      integer, intent(in)  :: ntot_soamode    !"last" mode on which soa is allowed to condense

      real(wp), intent(in   ), dimension( 1:nsoa, 1:max_mode ) :: uptkaer_soag_tmp
      real(wp), intent(in   ), dimension( 1:nsoa, 1:max_mode ) :: a_soa     ! soa aerosol mixrat (mol/mol at actual mw)
      logical , intent(in   ), dimension( 1:nsoa) :: skip_soamode
      real(wp), intent(in   ), dimension( 1:nsoa) :: g0_soa    ! ambient soa gas equilib mixrat (mol/mol at actual mw)
      real(wp), intent(in   ), dimension( 1:nsoa) :: g_soa     ! soa gas mixrat (mol/mol at actual mw)

      real(wp), intent(in   ), dimension( 1:max_mode) :: a_opoa  ! oxidized-poa aerosol mixrat (mol/mol at actual mw))
      real(wp), intent(inout) :: tcur         ! current time (s)
      real(wp), intent(inout  ), dimension( 1:nsoa, 1:max_mode ) :: g_star     ! soa gas mixrat (mol/mol at actual mw)
      real(wp), intent(inout  ), dimension(1:max_mode) :: a_ooa_sum_tmp  ! total ooa (=soa+opoa) in a mode
      real(wp), intent(inout  ), dimension( 1:nsoa, 1:max_mode ) :: sat  ! sat(m,ll) = g0_soa(ll)/a_ooa_sum_tmp(m) = g_star(m,ll)/a_soa(m,ll)
                                                   !    used by the numerical integration scheme -- it is not a saturation rato!
      real(wp), intent(out) :: dtcur        ! current integration timestep (s)

      real(wp), parameter :: a_min1 = 1.0e-20_wp
      real(wp), parameter :: g_min1 = 1.0e-20_wp
      integer :: ntot_poaspec
      integer :: ntot_soaspec
      real(wp) :: phi       ! "relative driving force"
      real(wp) :: dtmax, tmpa, tmpb
      integer :: ll, n


      ntot_poaspec = npoa
      ntot_soaspec = nsoa
      dtcur = 0.0_wp
      phi = 0.0_wp
      dtmax = 0.0_wp

! determine time step
      tmpa = 0.0_wp  ! time integration parameter for all soa species
      do n = 1, ntot_soamode
         if ( skip_soamode(n) ) cycle
         a_ooa_sum_tmp(n) = a_opoa(n) + sum( a_soa(1:ntot_soaspec,n) )
      end do
      do ll = 1, ntot_soaspec
         tmpb = 0.0_wp  ! time integration parameter for a single soa species
         do n = 1, ntot_soamode
            if ( skip_soamode(n) ) cycle
            sat(ll,n) = g0_soa(ll)/max( a_ooa_sum_tmp(n), a_min1 )
            g_star(ll,n) = sat(ll,n)*a_soa(ll,n)
            phi = (g_soa(ll) - g_star(ll,n))/max( g_soa(ll), g_star(ll,n), g_min1 )
            tmpb = tmpb + uptkaer_soag_tmp(ll,n)*abs(phi)
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

   end subroutine mam_soaexch_substep

end module mam_soaexch

