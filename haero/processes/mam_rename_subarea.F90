
! TODO: _alnsg_aer_ means "a natural log of standard deviation". Find more
! descriptive names. Are these already in haero? used in calcsize, refer to that
! PR. Extract standard dev, pass to _find_renaming_pairs_, take log of it in
! the routine.

! TODO: _dgnum_ means "diameter of particles". Not available in Haero atm.

! TODO: use dummy variables for parameters I do not know.

! TODO: send _alnsg_ to _compute_vol_to_num_ratio_. Precompute, pass to each
! routine.
! alnsg(imode) = log(model%modes(imode)%mean_std_dev)

#define rename_subarea_log write (*,*) 'mam_rename_subarea.F90', __LINE__,

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
            set_real_param

  ! TODO: Dummy variables which will eventually be fields extracted from the model
  ! TODO: Refer to TODOs at top of file
  ! ---
  real(wp) :: dgnumlo_aer(1), & ! TODO: get dim from model (model%modes)
              dgnumhi_aer(1), & ! TODO: get dim from model (model%modes)
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

  ! TODO: Move all parts of this subroutine that are only called once into the
  ! init routine
  subroutine run(model, t, dt, prognostics, atmosphere, diagnostics, tendencies)
    implicit none

    ! Parameters
    type(model_t), intent(in)         :: model
    real(wp), value, intent(in)       :: t
    real(wp), value, intent(in)       :: dt
    type(prognostics_t), intent(in)   :: prognostics
    type(atmosphere_t), intent(in)    :: atmosphere
    type(diagnostics_t), intent(in)   :: diagnostics
    type(tendencies_t), intent(inout) :: tendencies

    integer  :: nmodes
    integer  :: to_mode_of_mode(model%num_modes)

    ! total number of pairs to be found
    integer  :: num_pairs

    ! precomputed factors to be used later
    real(wp) :: sz_factor(model%num_modes), fmode_dist_tail_fac(model%num_modes)

    ! relaxed volume to num high and low ratio limits
    real(wp) :: v2n_lo_rlx(model%num_modes), v2n_hi_rlx(model%num_modes)

    ! log of diameter factor for distribution tail
    real(wp) :: ln_diameter_tail_fac(model%num_modes)

    ! cutoff (threshold) for deciding the do inter-mode transfer
    real(wp) :: diameter_cutoff(model%num_modes), &
                ln_dia_cutoff(model%num_modes)

    ! some limiters/factors
    real(wp) :: diameter_belowcutoff(model%num_modes), &
                dryvol_smallest(model%num_modes)

    !  alnsg(imode) = log(model%modes(imode)%mean_std_dev)
    real(wp) :: alnsg(model%num_modes)

    nmodes = model%num_modes
    to_mode_of_mode(:) = [0, 2, 0, 0]

    rename_subarea_log 'Calling find_renaming_pairs'
    rename_subarea_log 'Got num_modes = ', nmodes

    call populate_ln_of_std_dev(alnsg, model)

    ! TODO: new parameters needed for _find_renaming_pairs_; extract data from
    ! _model_ and pass to _find_renaming_pairs_. DONT pass model to subroutines.
    call find_renaming_pairs(nmodes, to_mode_of_mode, alnsg, &         ! input
       num_pairs, sz_factor, fmode_dist_tail_fac, v2n_lo_rlx, & ! output
       v2n_hi_rlx, ln_diameter_tail_fac, diameter_cutoff, &     ! output
       ln_dia_cutoff, diameter_belowcutoff, dryvol_smallest)    ! output

  end subroutine run

  ! alnsg is used throughout the original code. The model stores an array of
  ! the mean std deviations, which we can use to create a comperable array.
  ! This is not performant, but is being used temporarily during the porting
  ! process.
  !
  ! TODO: identify a better way to perform the same calculations as the original
  ! code, or remove this comment if this subroutine suffices.
  subroutine populate_ln_of_std_dev(alnsg, model)

    ! Parameters
    type(model_t), intent(in)     :: model
    real(wp),      intent(inout)  :: alnsg(:)

    ! Local variables
    integer :: imode

    rename_subarea_log 'Populating log of stddev array'

    do imode = 1, model%num_modes
      alnsg(imode) = log(model%modes(imode)%mean_std_dev)
    end do

  end subroutine populate_ln_of_std_dev

  subroutine finalize(model)
    implicit none

    ! Arguments
    type(model_t), intent(in) :: model
  end subroutine finalize

  subroutine find_renaming_pairs (nmodes, to_mode_of_mode, alnsg_aer, &    !input
       num_pairs, sz_factor, fmode_dist_tail_fac, v2n_lo_rlx, & !output
       v2n_hi_rlx, ln_diameter_tail_fac, diameter_cutoff, &     !output
       ln_dia_cutoff, diameter_belowcutoff, dryvol_smallest)    !output

    !arguments (intent-ins)

    real(wp), intent(in) :: alnsg_aer(:)
    integer,  intent(in) :: nmodes              !total number of modes
    integer,  intent(in) :: to_mode_of_mode(:)  !array carry info about the "to" mode of a particular mode

    !intent-outs
    integer,  intent(out) :: num_pairs         ! total number of pairs to be found
    real(wp), intent(out) :: sz_factor(:), fmode_dist_tail_fac(:) !precomputed factors to be used later
    real(wp), intent(out) :: v2n_lo_rlx(:), v2n_hi_rlx(:)         !relaxed volume to num high and low ratio limits
    real(wp), intent(out) :: ln_diameter_tail_fac(:)              !log of diameter factor for distribution tail
    real(wp), intent(out) :: diameter_cutoff(:), ln_dia_cutoff(:) !cutoff (threshold) for deciding the  do inter-mode transfer
    real(wp), intent(out) :: diameter_belowcutoff(:), dryvol_smallest(:) ! some limiters/factors

    !local variables
    integer  :: to_mode, from_mode, imode
    real(wp) :: alnsg_for_current_mode

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
    if (all(to_mode_of_mode(:)<=0)) then
      rename_subarea_log 'Got all(to_mode_of_mode(:)<=0), returning early'
      return
    endif

    ! Find >=1 pair
    do imode = 1, nmodes
       rename_subarea_log 'Searching for pairs with imode = ', imode

       to_mode   = to_mode_of_mode(imode) ! transfer "to" mode for mode "imode"

       ! if to_mode is <=0, transfer is not possible for this mode, cycle the
       ! loop for the next mode
       if(to_mode <= 0) then
         rename_subarea_log 'Got to_mode <= 0, skipping this mode.'
         cycle
       end if

       ! transfer "from" mode is the current mode (i.e. imode)
       from_mode = imode

       ! log of stddev for current mode
       alnsg_for_current_mode = alnsg_aer(from_mode)

       !^^At this point, we know that particles can be tranfered from the
       ! "from_mode" to "to_mode". "from_mode" is the current mode (i.e. imode)

       ! update number of pairs found so far
       num_pairs = num_pairs + 1

       !-------------------------------------------------------
       ! precompute common factors to be used later
       !-------------------------------------------------------

       ! size factor for "from mode"
       call compute_size_factor (from_mode, alnsg_for_current_mode, sz_factor)

       ! size factor for "to mode"
       call compute_size_factor (to_mode, alnsg_for_current_mode, sz_factor)

       !------------------------------------------------------------------------
       ! We compute few factors below for the "from_mode", which will be used
       ! for inter-mode particle transfer
       !------------------------------------------------------------------------

       ! factor for computing distribution tails of the "from mode"
       fmode_dist_tail_fac(from_mode) = sqrt_half/alnsg_for_current_mode

       dryvol_smallest(from_mode) = smallest_dryvol_value
       ! compute volume to number high and low limits with relaxation
       ! coefficients (watch out for repeated calculations)

       ! TODO: see calcsize pr for these values, how to extract from model
       v2n_lo_rlx(from_mode) = &
         compute_vol_to_num_ratio(from_mode, alnsg_for_current_mode, dgnumlo_aer) * frelax

       v2n_hi_rlx(from_mode) = &
         compute_vol_to_num_ratio(from_mode, alnsg_for_current_mode, dgnumhi_aer) / frelax

       ! A factor for computing diameter at the tails of the distribution
       ln_diameter_tail_fac(from_mode) = 3.0 * (alnsg_for_current_mode**2)

       ! Cut-off (based on geometric mean) for making decision to do inter-mode transfers

       ! TODO: use dummy values for _dgnum_aer_, or assign to dgnum_low for the
       ! moment. Have to figure out how to compute this. We will extract from
       ! model at some point.
       diameter_cutoff(from_mode) = sqrt(   &
          dgnum_aer(from_mode)*exp(1.5*(alnsg_for_current_mode**2)) *   &
          dgnum_aer(to_mode)*exp(1.5*(alnsg_aer(to_mode)**2)) )

       ln_dia_cutoff(from_mode) = log(diameter_cutoff(from_mode)) !log of cutt-off
       diameter_belowcutoff(from_mode) = 0.99*diameter_cutoff(from_mode) !99% of the cutoff

    enddo

  end subroutine find_renaming_pairs


  subroutine compute_size_factor(imode, alnsg, size_factor)
    ! Compute size factor for a mode
    use haero_constants, only: pi_sixth
    implicit none

    integer,  intent(in) :: imode     !mode number
    real(wp), intent(in) :: alnsg
    real(wp), intent(inout) :: size_factor(:) !size factor

    size_factor(imode) = (pi_sixth)*exp(4.5*(alnsg**2))

  end subroutine compute_size_factor


  pure function compute_vol_to_num_ratio(imode, alnsg, diameter) result(v2n)
    !compute volume to number ratio for a mode
    use haero_constants, only: pi_sixth
    implicit none
    integer,  intent(in) :: imode
    real(wp), intent(in) :: alnsg
    real(wp), intent(in) :: diameter(:) ![m]

    real(wp) :: v2n !return value

    v2n = ( 1._wp / ( (pi_sixth)* &
         (diameter(imode)**3._wp)*exp(4.5_wp*alnsg**2._wp) ) )

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
