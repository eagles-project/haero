module mam_rename_subarea

  use haero_precision, only: wp
  use haero, only: model_t, aerosol_species_t, gas_species_t, &
       prognostics_t, atmosphere_t, diagnostics_t, tendencies_t

  implicit none
  public :: init, &
            run, &
            finalize, &
            set_integer_param, &
            set_logical_param, &
            set_real_param, &
            find_renaming_pairs

  ! Dummy variables which will eventually be fields extracted from the model
  ! ---
  real(wp) :: dgnumlo_aer(1), &
              dgnumhi_aer(1), &
              alnsg_aer(1), &
              dgnum_aer(1), &
              fac_m2v_aer(1)

  integer :: max_mode, &
             rename_method_optaa, &
             ntot_amode, &
             naer, &
             max_aer
  ! ---

contains
  subroutine init(model)

    implicit none

    ! Arguments
    type(model_t), intent(in) :: model
  end subroutine init

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

    print *, 'in f routine'

  end subroutine run

  subroutine finalize(model)
    implicit none

    ! Arguments
    type(model_t), intent(in) :: model
  end subroutine finalize


      subroutine mam_rename_1subarea(                               &
         nstep,             lchnk,                                  &
         i,                 k,                jsub,                 &
         latndx,            lonndx,           lund,                 &
         iscldy_subarea,                                            &
         mtoo_renamexf,                                             &
         n_mode,                                                    &
         qnum_cur,                                                  &
         qaer_cur,          qaer_del_grow4rnam,                     &
         qwtr_cur,                                                  &
         qnumcw_cur,                                                &
         qaercw_cur,        qaercw_del_grow4rnam                    )

       ! Note that the original code used error functions from shr_spfn_mod
       ! or error_function modules. These calls were replaced with the GNU erf
       ! function.

      logical,  intent(in)    :: iscldy_subarea        ! true if sub-area is cloudy
      integer,  intent(in)    :: nstep                 ! model time-step number
      integer,  intent(in)    :: lchnk                 ! chunk identifier
      integer,  intent(in)    :: i, k                  ! column and level indices
      integer,  intent(in)    :: jsub                  ! sub-area index
      integer,  intent(in)    :: latndx, lonndx        ! lat and lon indices
      integer,  intent(in)    :: lund                  ! logical unit for diagnostic output
      integer,  intent(in)    :: mtoo_renamexf(max_mode)
      integer,  intent(in)    :: n_mode                ! current number of modes (including temporary)

      real(wp), intent(inout), dimension( 1:max_mode ) :: &
         qnum_cur
      real(wp), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_cur
      real(wp), intent(in   ), dimension( 1:max_aer, 1:max_mode ) :: &
         qaer_del_grow4rnam
      real(wp), intent(inout), dimension( 1:max_mode ) :: &
         qwtr_cur

      real(wp), intent(inout), optional, dimension( 1:max_mode ) :: &
         qnumcw_cur
      real(wp), intent(inout), optional, dimension( 1:max_aer, 1:max_mode ) :: &
         qaercw_cur
      real(wp), intent(in   ), optional, dimension( 1:max_aer, 1:max_mode ) :: &
         qaercw_del_grow4rnam


! !DESCRIPTION:
! computes TMR (tracer mixing ratio) tendencies for "mode renaming"
!    during a continuous growth process
! currently this transfers number and mass (and surface) from the aitken
!    to accumulation mode after gas condensation or stratiform-cloud
!    aqueous chemistry
! (convective cloud aqueous chemistry not yet implemented)
!
! !REVISION HISTORY:
!

! local variables
      integer :: iaer
      integer :: mfrm, mtoo
      integer :: n, npair

      integer, parameter :: ldiag1 = 0

      real(wp), parameter :: onethird = 1.0_wp/3.0_wp

      real(wp) :: deldryvol_a(ntot_amode)
      real(wp) :: deldryvol_c(ntot_amode)
      real(wp) :: dp_belowcut(max_mode)
      real(wp) :: dp_cut(max_mode)
      real(wp) :: dgn_aftr, dgn_xfer
      real(wp) :: dgn_t_new, dgn_t_old, dgn_t_oldaa
      real(wp) :: dryvol_t_del, dryvol_t_new
      real(wp) :: dryvol_t_old, dryvol_t_oldaa, dryvol_t_oldbnd
      real(wp) :: dryvol_a(ntot_amode)
      real(wp) :: dryvol_c(ntot_amode)
      real(wp) :: dryvol_smallest(ntot_amode)
      real(wp) :: factoraa(ntot_amode)
      real(wp) :: factoryy(ntot_amode)
      real(wp) :: lndp_cut(max_mode)
      real(wp) :: lndgn_new, lndgn_old
      real(wp) :: lndgv_new, lndgv_old
      real(wp) :: num_t_old, num_t_oldbnd
      real(wp) :: tailfr_volnew, tailfr_volold
      real(wp) :: tailfr_numnew, tailfr_numold
      real(wp) :: tmpa, tmpb, tmpd
      real(wp) :: tmp_alnsg2(max_mode)
      real(wp) :: v2nhirlx(ntot_amode), v2nlorlx(ntot_amode)
      real(wp) :: xfercoef, xfertend
      real(wp) :: xferfrac_vol, xferfrac_num, xferfrac_max
      real(wp) :: yn_tail, yv_tail

      xferfrac_max = 1.0_wp - 10.0_wp*epsilon(1.0_wp)   ! 1-eps !FIXME: This can be a parameter???

      !------------------------------------------------------------------------
      !Find mapping between different modes, so that we can move aerosol
      !particles from one mode to another
      !------------------------------------------------------------------------

      !FIXME: All the arrays in find_renaming_pairs subroutine call should be
      !initialized to HUGE or NaNs as they are partially populated

      ! Find (from->to) pairs of modes which can participate in inter-mode particle transfer
      call find_renaming_pairs (ntot_amode, mtoo_renamexf, & !input
           npair, factoraa, factoryy, v2nlorlx, &            !output
           v2nhirlx, tmp_alnsg2, dp_cut, &                   !output
           lndp_cut, dp_belowcut, dryvol_smallest)           !output

      if (npair <= 0) return ! if no transfer required, return

      !^^^^^^^^^^^^^^^^^^ BSINGH - ENDS REFACTOR ^^^^^^^^^^^^^^^^^^^^^^^



      !BSINGH- original comments to be adjusted
! calculate variable used in the renamingm mode" of each renaming pair
! also compute dry-volume change during the continuous growth process

! dryvol_smallest is a very small volume mixing ratio (m3-AP/kmol-air)
! used for avoiding overflow.  it corresponds to dp = 1 nm
! and number = 1e-5 #/mg-air ~= 1e-5 #/cm3-air

      !Original Comments ENDS

! compute aerosol dry-volume for the "from mode" of each renaming pair
! also compute dry-volume change during the continuous growth process
      do n = 1, ntot_amode
         mtoo = mtoo_renamexf(n)
         if (mtoo <= 0) cycle

         tmpa = 0.0_wp ; tmpb = 0.0_wp
         do iaer = 1, naer
!   fac_m2v_aer converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
            tmpa = tmpa + qaer_cur(iaer,n)*fac_m2v_aer(iaer)
            tmpb = tmpb + qaer_del_grow4rnam(iaer,n)*fac_m2v_aer(iaer)
         end do
         dryvol_a(n) = tmpa-tmpb ! dry volume before growth
         deldryvol_a(n) = tmpb   ! change to dry volume due to growth

         if ( iscldy_subarea ) then
         tmpa = 0.0_wp ; tmpb = 0.0_wp
         do iaer = 1, naer
!   fac_m2v_aer converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
            tmpa = tmpa + qaercw_cur(iaer,n)*fac_m2v_aer(iaer)
            tmpb = tmpb + qaercw_del_grow4rnam(iaer,n)*fac_m2v_aer(iaer)
         end do
         dryvol_c(n) = tmpa-tmpb
         deldryvol_c(n) = tmpb
         end if ! ( iscldy_subarea ) then

      end do


!
!   loop over renaming pairs
!
mainloop1_ipair:  do n = 1, ntot_amode

      mfrm = n
      mtoo = mtoo_renamexf(n)
      if (mtoo <= 0) cycle mainloop1_ipair

!   dryvol_t_old is the old total (a+c) dry-volume for the "from" mode
!      in m^3-AP/kmol-air
!   dryvol_t_new is the new total dry-volume
!      (old/new = before/after the continuous growth)
!   num_t_old is total number in particles/kmol-air
      if ( iscldy_subarea ) then
         dryvol_t_old = dryvol_a(mfrm) + dryvol_c(mfrm)
         dryvol_t_del = deldryvol_a(mfrm) + deldryvol_c(mfrm)
         num_t_old = (qnum_cur(mfrm) + qnumcw_cur(mfrm))
      else
         dryvol_t_old = dryvol_a(mfrm)
         dryvol_t_del = deldryvol_a(mfrm)
         num_t_old = qnum_cur(mfrm)
      end if
      dryvol_t_new = dryvol_t_old + dryvol_t_del

!   no renaming if dryvol_t_new ~ 0 or dryvol_t_del ~ 0
      if (dryvol_t_new .le. dryvol_smallest(mfrm)) cycle mainloop1_ipair
      dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )
      if (rename_method_optaa .ne. 40) then
         if (dryvol_t_del .le. 1.0e-6*dryvol_t_oldbnd) cycle mainloop1_ipair
      end if

      num_t_old = max( 0.0_wp, num_t_old )
      dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )
      num_t_oldbnd = min( dryvol_t_oldbnd*v2nlorlx(mfrm), num_t_old )
      num_t_oldbnd = max( dryvol_t_oldbnd*v2nhirlx(mfrm), num_t_oldbnd )

!   no renaming if dgnum < "base" dgnum,
      dgn_t_new = (dryvol_t_new/(num_t_oldbnd*factoraa(mfrm)))**onethird
      if (dgn_t_new .le. dgnum_aer(mfrm)) cycle mainloop1_ipair

!   compute new fraction of number and mass in the tail (dp > dp_cut)
      lndgn_new = log( dgn_t_new )
      lndgv_new = lndgn_new + tmp_alnsg2(mfrm)
      yn_tail = (lndp_cut(mfrm) - lndgn_new)*factoryy(mfrm)
      yv_tail = (lndp_cut(mfrm) - lndgv_new)*factoryy(mfrm)
      tailfr_numnew = 0.5_wp*erf( yn_tail )
      tailfr_volnew = 0.5_wp*erf( yv_tail )

!   compute old fraction of number and mass in the tail (dp > dp_cut)
      dgn_t_old =   &
            (dryvol_t_oldbnd/(num_t_oldbnd*factoraa(mfrm)))**onethird
      dgn_t_oldaa = dgn_t_old
      dryvol_t_oldaa = dryvol_t_old

      if (rename_method_optaa .eq. 40) then
         if (dgn_t_old .gt. dp_belowcut(mfrm)) then
            ! this revised volume corresponds to dgn_t_old == dp_belowcut, and same number conc
            dryvol_t_old = dryvol_t_old * (dp_belowcut(mfrm)/dgn_t_old)**3
            dgn_t_old = dp_belowcut(mfrm)
         end if
         if ((dryvol_t_new-dryvol_t_old) .le. 1.0e-6_wp*dryvol_t_oldbnd) cycle mainloop1_ipair
      else if (dgn_t_new .ge. dp_cut(mfrm)) then
!         if dgn_t_new exceeds dp_cut, use the minimum of dgn_t_old and
!         dp_belowcut to guarantee some transfer
          dgn_t_old = min( dgn_t_old, dp_belowcut(mfrm) )
      end if
      lndgn_old = log( dgn_t_old )
      lndgv_old = lndgn_old + tmp_alnsg2(mfrm)
      yn_tail = (lndp_cut(mfrm) - lndgn_old)*factoryy(mfrm)
      yv_tail = (lndp_cut(mfrm) - lndgv_old)*factoryy(mfrm)
      tailfr_numold = 0.5_wp*erf( yn_tail )
      tailfr_volold = 0.5_wp*erf( yv_tail )

!   transfer fraction is difference between new and old tail-fractions
!   transfer fraction for number cannot exceed that of mass
      tmpa = tailfr_volnew*dryvol_t_new - tailfr_volold*dryvol_t_old
      if (tmpa .le. 0.0_wp) cycle mainloop1_ipair

      xferfrac_vol = min( tmpa, dryvol_t_new )/dryvol_t_new
      xferfrac_vol = min( xferfrac_vol, xferfrac_max )
      xferfrac_num = tailfr_numnew - tailfr_numold
      xferfrac_num = max( 0.0_wp, min( xferfrac_num, xferfrac_vol ) )
#if ( defined( CAMBOX_ACTIVATE_THIS ) )
      if ( ldiag98 ) write(lun98,'(/a,2i3,1p,10e11.3)') &
         'rename i,k, xf n/v', i, k, xferfrac_num, xferfrac_vol
#endif

#if ( defined( CAMBOX_NEVER_ACTIVATE_THIS ) )
!   diagnostic output start ----------------------------------------
       if (ldiag1 > 0) then
       icol_diag = -1
       if ((lonndx(i) == 37) .and. (latndx(i) == 23)) icol_diag = i
       if ((i == icol_diag) .and. (mod(k-1,5) == 0)) then
 !      write(lund,97010) fromwhere, nstep, lchnk, i, k, ipair
       write(lund,97010) fromwhere, nstep, latndx(i), lonndx(i), k, ipair
       write(lund,97020) 'drv olda/oldbnd/old/new/del',   &
             dryvol_t_oldaa, dryvol_t_oldbnd, dryvol_t_old, dryvol_t_new, dryvol_t_del
       write(lund,97020) 'num old/oldbnd, dgnold/new ',   &
             num_t_old, num_t_oldbnd, dgn_t_old, dgn_t_new
       write(lund,97020) 'tailfr v_old/new, n_old/new',   &
             tailfr_volold, tailfr_volnew, tailfr_numold, tailfr_numnew
       tmpa = max(1.0d-10,xferfrac_vol) / max(1.0d-10,xferfrac_num)
       dgn_xfer = dgn_t_new * tmpa**onethird
       tmpa = max(1.0d-10,(1.0d0-xferfrac_vol)) /   &
               max(1.0d-10,(1.0d0-xferfrac_num))
       dgn_aftr = dgn_t_new * tmpa**onethird
       write(lund,97020) 'xferfrac_v/n; dgn_xfer/aftr',   &
             xferfrac_vol, xferfrac_num, dgn_xfer, dgn_aftr
 !97010      format( / 'RENAME ', a, '  nx,lc,i,k,ip', i8, 4i4 )
 97010      format( / 'RENAME ', a, '  nx,lat,lon,k,ip', i8, 4i4 )
 97020      format( a, 6(1pe15.7) )
       end if
       end if ! (ldiag1 > 0)
!   diagnostic output end   ------------------------------------------
#endif


!
!   compute changes to number and species masses
!
      tmpd = qnum_cur(mfrm)*xferfrac_num
      qnum_cur(mfrm) = qnum_cur(mfrm) - tmpd
      qnum_cur(mtoo) = qnum_cur(mtoo) + tmpd
      do iaer = 1, naer
         tmpd = qaer_cur(iaer,mfrm)*xferfrac_vol
         qaer_cur(iaer,mfrm) = qaer_cur(iaer,mfrm) - tmpd
         qaer_cur(iaer,mtoo) = qaer_cur(iaer,mtoo) + tmpd
      end do ! iaer

      if ( iscldy_subarea ) then
      tmpd = qnumcw_cur(mfrm)*xferfrac_num
      qnumcw_cur(mfrm) = qnumcw_cur(mfrm) - tmpd
      qnumcw_cur(mtoo) = qnumcw_cur(mtoo) + tmpd
      do iaer = 1, naer
         tmpd = qaercw_cur(iaer,mfrm)*xferfrac_vol
         qaercw_cur(iaer,mfrm) = qaercw_cur(iaer,mfrm) - tmpd
         qaercw_cur(iaer,mtoo) = qaercw_cur(iaer,mtoo) + tmpd
      end do ! iaer
      end if ! ( iscldy_subarea ) then


#if ( defined( CAMBOX_NEVER_ACTIVATE_THIS ) )
!   diagnostic output start ----------------------------------------
                if (ldiag1 > 0) then
                if ((i == icol_diag) .and. (mod(k-1,5) == 0)) then
                  if (lstooa .gt. 0) then
                    write(lund,'(a,i4,2(2x,a),1p,10e14.6)') 'RENAME qdels', iq,   &
                        cnst_name(lsfrma+loffset), cnst_name(lstooa+loffset),   &
                        deltat*dqdt(i,k,lsfrma), deltat*(dqdt(i,k,lsfrma) - xfertend),   &
                        deltat*dqdt(i,k,lstooa), deltat*(dqdt(i,k,lstooa) + xfertend)
                  else
                    write(lund,'(a,i4,2(2x,a),1p,10e14.6)') 'RENAME qdels', iq,   &
                        cnst_name(lsfrma+loffset), cnst_name(lstooa+loffset),   &
                        deltat*dqdt(i,k,lsfrma), deltat*(dqdt(i,k,lsfrma) - xfertend)
                  end if
                end if
                end if
!   diagnostic output end   ------------------------------------------
#endif


      end do mainloop1_ipair


      return
      end subroutine mam_rename_1subarea


  subroutine find_renaming_pairs (nmodes, to_mode_of_mode, &    !input
       num_pairs, sz_factor, fmode_dist_tail_fac, v2n_lo_rlx, & !output
       v2n_hi_rlx, ln_diameter_tail_fac, diameter_cutoff, &     !output
       ln_dia_cutoff, diameter_belowcutoff, dryvol_smallest)    !output

    !arguments (intent-ins)
    integer, intent(in) :: nmodes              !total number of modes
    integer, intent(in) :: to_mode_of_mode(:)  !array carry info about the "to" mode of a particular mode

    !intent-outs
    integer,  intent(out) :: num_pairs         ! total number of pairs to be found
    real(wp), intent(out) :: sz_factor(:), fmode_dist_tail_fac(:) !precomputed factors to be used later
    real(wp), intent(out) :: v2n_lo_rlx(:), v2n_hi_rlx(:)         !relaxed volume to num high and low ratio limits
    real(wp), intent(out) :: ln_diameter_tail_fac(:)              !log of diameter factor for distribution tail
    real(wp), intent(out) :: diameter_cutoff(:), ln_dia_cutoff(:) !cutoff (threshold) for deciding the  do inter-mode transfer
    real(wp), intent(out) :: diameter_belowcutoff(:), dryvol_smallest(:) ! some limiters/factors

    !local variables
    integer :: to_mode, from_mode, imode

    !some parameters
    real(wp), parameter :: sqrt_half = sqrt(0.5)
    real(wp), parameter :: frelax = 27.0_wp !(3^3)
    real(wp), parameter :: smallest_dryvol_value = 1.0e-25

  !number of pairs allowed to do inter-mode particle transfer
  ! (e.g. if we have a pair "mode_1<-->mode_2", mode_1 and mode_2 can participate in
  ! inter-mode aerosol particle transfer where like particles in mode_1 can be
  ! transferred to mode_2 and vice-versa)
  num_pairs = 0 ! Let us assume there are none to start with

  !if there can be no possible pairs, just return
  if (all(to_mode_of_mode(:)<=0)) return

  !Go through all the modes to find if we have atleast one or more than one pairs
  do imode = 1, nmodes
     to_mode   = to_mode_of_mode(imode) ! transfer "to" mode for mode "imode"

     !if to_mode is <=0, transfer is not possible for this mode, cycle the loop for the next mode
     if(to_mode <= 0)cycle

     from_mode = imode                  ! transfer "from" mode is the current mode (i.e. imode)

     !^^At this point, we know that particles can be tranfered from the
     ! "from_mode" to "to_mode". "from_mode" is the current mode (i.e. imode)

     !update number of pairs found so far
     num_pairs = num_pairs + 1    !increment npair

     !-------------------------------------------------------
     !now precompute some common factors to be used later
     !-------------------------------------------------------

     call compute_size_factor (from_mode, sz_factor) !size factor for "from mode"
     call compute_size_factor (to_mode,   sz_factor) !size factor for "to mode"

     !---------------------------------------------------------------------------------------------------------
     ! We compute few factors below for the "from_mode", which will be used for inter-mode particle transfer
     !---------------------------------------------------------------------------------------------------------

     fmode_dist_tail_fac(from_mode) = sqrt_half/alnsg_aer(from_mode) !factor for computing distribution tails of the  "from mode"

     dryvol_smallest(from_mode) = smallest_dryvol_value
     !compute volume to number high and low limits with relaxation coefficients (watch out for repeated calculations)
     v2n_lo_rlx(from_mode) = compute_vol_to_num_ratio(from_mode, dgnumlo_aer) * frelax
     v2n_hi_rlx(from_mode) = compute_vol_to_num_ratio(from_mode, dgnumhi_aer) / frelax

     !A factor for computing diameter at the tails of the distribution
     ln_diameter_tail_fac(from_mode) = 3.0 * (alnsg_aer(from_mode)**2)

     !Cut-off (based on geometric mean) for making decision to do inter-mode transfers
     diameter_cutoff(from_mode) = sqrt(   &
        dgnum_aer(from_mode)*exp(1.5*(alnsg_aer(from_mode)**2)) *   &
        dgnum_aer(to_mode)*exp(1.5*(alnsg_aer(to_mode)**2)) )

     ln_dia_cutoff(from_mode) = log(diameter_cutoff(from_mode)) !log of cutt-off
     diameter_belowcutoff(from_mode) = 0.99*diameter_cutoff(from_mode) !99% of the cutoff

  enddo

end subroutine find_renaming_pairs


  subroutine compute_size_factor(imode, size_factor)
    ! Compute size factor for a mode
    use haero_constants, only: pi_sixth
    implicit none

    integer,  intent(in) :: imode     !mode number
    real(wp), intent(inout) :: size_factor(:) !size factor

    size_factor(imode) = (pi_sixth/6.)*exp(4.5*(alnsg_aer(imode)**2))

  end subroutine compute_size_factor


  pure function compute_vol_to_num_ratio(imode, diameter) result(v2n)
    !compute volume to number ratio for a mode
    use haero_constants, only: pi_sixth
    implicit none
    integer,  intent(in) :: imode
    real(wp), intent(in) :: diameter(:) ![m]

    real(wp) :: v2n !return value

    v2n = ( 1._wp / ( (pi_sixth/6._wp)* &
         (diameter(imode)**3._wp)*exp(4.5_wp*alnsg_aer(imode)**2._wp) ) )

  end function compute_vol_to_num_ratio


  subroutine set_integer_param(name, val)
    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: val

    ! No integers to set!
  end subroutine set_integer_param


  subroutine set_logical_param(name, val)
    implicit none

    character(len=*), intent(in) :: name
    logical, intent(in) :: val

    ! No logicals to set!
  end subroutine set_logical_param


  subroutine set_real_param(name, val)
    implicit none

    character(len=*), intent(in) :: name
    real(wp), intent(in) :: val

    ! No reals to set!
  end subroutine set_real_param


end module mam_rename_subarea
