
! TODO Evaluate whether the _initialize_ methods are helpful or just extra
! overhead.

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

    ! --- Parameters
    type(model_t), intent(in)         :: model
    real(wp), value, intent(in)       :: t
    real(wp), value, intent(in)       :: dt
    type(prognostics_t), intent(in)   :: prognostics
    type(atmosphere_t), intent(in)    :: atmosphere
    type(diagnostics_t), intent(in)   :: diagnostics
    type(tendencies_t), intent(inout) :: tendencies

    ! --- Local variables

    ! Total number of modes
    integer  :: nmodes

    ! Contains information about the destination mode for a given mode
    integer  :: dest_mode_of_mode(model%num_modes)

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

    ! Upper and lower limits for diameters
    real(wp) :: dgnumlo_aer(model%num_modes), &
                dgnumhi_aer(model%num_modes), &
                dgnum_aer(model%num_modes)

    ! Natural log of the mean std dev for each mode
    real(wp) :: alnsg(model%num_modes)

    nmodes = model%num_modes
    dest_mode_of_mode(:) = [0, 2, 0, 0]

    call initialize_diameters(dgnumlo_aer, dgnumhi_aer, dgnum_aer, model)

    call initialize_ln_of_std_dev(alnsg, model)

    ! TODO: new parameters needed for _find_renaming_pairs_; extract data from
    ! _model_ and pass to _find_renaming_pairs_. DONT pass model to subroutines.
    call find_renaming_pairs(nmodes, dest_mode_of_mode, alnsg, &    ! input
        dgnumlo_aer, dgnumhi_aer, dgnum_aer,                   &    ! input
        num_pairs, sz_factor, fmode_dist_tail_fac, v2n_lo_rlx, &    ! output
        v2n_hi_rlx, ln_diameter_tail_fac, diameter_cutoff,     &    ! output
        ln_dia_cutoff, diameter_belowcutoff, dryvol_smallest)       ! output

    rename_subarea_log 'Final results:'
    rename_subarea_log 'nmodes=', nmodes
    rename_subarea_log 'dest_mode_of_mode=', dest_mode_of_mode
    rename_subarea_log 'alnsg=', alnsg
    rename_subarea_log 'dgnumlo_aer=', dgnumlo_aer
    rename_subarea_log 'dgnumhi_aer=', dgnumhi_aer
    rename_subarea_log 'dgnum_aer=', dgnum_aer
    rename_subarea_log 'num_pairs=', num_pairs
    rename_subarea_log 'sz_factor=', sz_factor
    rename_subarea_log 'fmode_dist_tail_fac=', fmode_dist_tail_fac
    rename_subarea_log 'v2n_lo_rlx=', v2n_lo_rlx
    rename_subarea_log 'v2n_hi_rlx=', v2n_hi_rlx

  end subroutine run

  subroutine initialize_diameters(dgnumlo_aer, dgnumhi_aer, dgnum_aer, model)
    
    type(model_t), intent(in)    :: model
    real(wp),      intent(inout) :: dgnumlo_aer(:), &
                                    dgnumhi_aer(:), &
                                    dgnum_aer(:)

    ! Initialize min and max diameters
    dgnumlo_aer(:) = model%modes(:)%min_diameter
    dgnumhi_aer(:) = model%modes(:)%max_diameter

    ! Initialize this to the minimum diameter for now.
    ! TODO: this will be updated with the correct calculation later.
    dgnum_aer(:) = dgnumlo_aer(:)

  end subroutine initialize_diameters

  ! alnsg is used throughout the original code. The model stores an array of
  ! the mean std deviations, which we can use to create a comperable array.
  ! This is not performant, but is being used temporarily during the porting
  ! process.
  !
  ! TODO: identify a better way to perform the same calculations as the original
  ! code, or remove this comment if this subroutine suffices.
  subroutine initialize_ln_of_std_dev(alnsg, model)

    ! Parameters
    type(model_t), intent(in)     :: model
    real(wp),      intent(inout)  :: alnsg(:)

    rename_subarea_log 'Populating log of stddev array:'

    alnsg(:) = log(model%modes(:)%mean_std_dev)

  end subroutine initialize_ln_of_std_dev

  subroutine finalize(model)
    implicit none

    ! Arguments
    type(model_t), intent(in) :: model
  end subroutine finalize

  subroutine find_renaming_pairs (nmodes, dest_mode_of_mode,  & ! input
       alnsg_aer, dgnumlo_aer, dgnumhi_aer, dgnum_aer,        & ! input
       num_pairs, sz_factor, fmode_dist_tail_fac, v2n_lo_rlx, & ! output
       v2n_hi_rlx, ln_diameter_tail_fac, diameter_cutoff,     & ! output
       ln_dia_cutoff, diameter_belowcutoff, dryvol_smallest)    ! output

    ! --- arguments (intent-ins)

    ! Natural log of std deviation
    real(wp), intent(in) :: alnsg_aer(:)

    ! Diameters
    real(wp), intent(in) :: dgnumlo_aer(:)
    real(wp), intent(in) :: dgnumhi_aer(:)
    real(wp), intent(in) :: dgnum_aer(:)

    ! total number of modes
    integer,  intent(in) :: nmodes

    ! Array for information about the destination mode of a particular mode
    integer,  intent(in) :: dest_mode_of_mode(:)

    ! --- intent-outs

    ! total number of pairs to be found
    integer,  intent(out) :: num_pairs

    ! precomputed factors to be used later
    real(wp), intent(out) :: sz_factor(:), fmode_dist_tail_fac(:)

    ! relaxed volume to num high and low ratio limits
    real(wp), intent(out) :: v2n_lo_rlx(:), v2n_hi_rlx(:)

    ! log of diameter factor for distribution tail
    real(wp), intent(out) :: ln_diameter_tail_fac(:)

    ! cutoff (threshold) for deciding whether or not to do inter-mode transfer
    real(wp), intent(out) :: diameter_cutoff(:), ln_dia_cutoff(:)

    ! limiters/factors
    real(wp), intent(out) :: diameter_belowcutoff(:), dryvol_smallest(:)

    ! --- local variables

    ! Destination mode
    integer  :: dest_mode

    ! Source mode
    integer  :: src_mode

    ! Iterator for looping over modes
    integer  :: imode

    ! Ln of stddev for source mode
    real(wp) :: alnsg_for_current_mode

    ! --- Parameters
    real(wp), parameter :: sqrt_half = sqrt(0.5)
    real(wp), parameter :: frelax = 27.0_wp !(3^3)
    real(wp), parameter :: smallest_dryvol_value = 1.0e-25

    ! number of pairs allowed to do inter-mode particle transfer
    ! (e.g. if we have a pair "mode_1<-->mode_2", mode_1 and mode_2 can participate in
    ! inter-mode aerosol particle transfer where like particles in mode_1 can be
    ! transferred to mode_2 and vice-versa)
    !
    ! Let us assume there are none to start with.
    num_pairs = 0

    ! if there can be no possible pairs, just return
    if (all(dest_mode_of_mode(:)<=0)) then
      rename_subarea_log 'Found no possible mode pairs, returning early'
      return
    endif

    ! Find >=1 pair
    do imode = 1, nmodes
      rename_subarea_log 'Searching for pairs with imode=', imode

      ! Destination mode for mode _imode_
      dest_mode = dest_mode_of_mode(imode)

      ! if dest_mode is <=0, transfer is not possible for this mode. cycle the
      ! loop for the next mode
      if(dest_mode <= 0) then
        rename_subarea_log 'Got dest_mode <= 0, skipping this mode.'
        cycle
      end if

      ! Source mode is the current mode (i.e. _imode_)
      src_mode = imode

      ! log of stddev for current mode
      alnsg_for_current_mode = alnsg_aer(src_mode)

      !^^At this point, we know that particles can be tranfered from the
      ! _src_mode_ to _dest_mode_. _src_mode_ is the current mode (i.e. imode)

      ! update number of pairs found so far
      num_pairs = num_pairs + 1

      !-------------------------------------------------------
      ! precompute common factors to be used later
      !-------------------------------------------------------

      ! size factor for _src_mode_
      call compute_size_factor (src_mode, alnsg_for_current_mode, sz_factor)

      ! size factor for _dest_mode_
      call compute_size_factor (dest_mode, alnsg_for_current_mode, sz_factor)

      !------------------------------------------------------------------------
      ! We compute few factors below for the _src_mode_, which will be used
      ! for inter-mode particle transfer
      !------------------------------------------------------------------------

      ! factor for computing distribution tails of the _src_mode_
      fmode_dist_tail_fac(src_mode) = sqrt_half/alnsg_for_current_mode

      dryvol_smallest(src_mode) = smallest_dryvol_value
      ! compute volume to number high and low limits with relaxation
      ! coefficients (watch out for repeated calculations)

      ! TODO: see calcsize pr for these values, how to extract from model
      v2n_lo_rlx(src_mode) = compute_vol_to_num_ratio(src_mode, &
                               alnsg_for_current_mode, &
                               dgnumlo_aer) &
                               * frelax

      v2n_hi_rlx(src_mode) = compute_vol_to_num_ratio(src_mode, &
                               alnsg_for_current_mode, &
                               dgnumhi_aer) &
                               / frelax

      ! A factor for computing diameter at the tails of the distribution
      ln_diameter_tail_fac(src_mode) = 3.0 * (alnsg_for_current_mode**2)

      ! Cut-off (based on geometric mean) for making decision to do inter-mode
      ! transfers

      ! TODO: use dummy values for _dgnum_aer_, or assign to dgnum_low for the
      ! moment. Have to figure out how to compute this. We will extract from
      ! model at some point.
      diameter_cutoff(src_mode) = sqrt(   &
         dgnum_aer(src_mode)*exp(1.5*(alnsg_for_current_mode**2)) *   &
         dgnum_aer(dest_mode)*exp(1.5*(alnsg_aer(dest_mode)**2)) )

      ln_dia_cutoff(src_mode) = log(diameter_cutoff(src_mode)) !log of cutt-off
      diameter_belowcutoff(src_mode) = 0.99*diameter_cutoff(src_mode) !99% of the cutoff

    enddo

  end subroutine find_renaming_pairs


  ! Compute size factor for a mode
  subroutine compute_size_factor(imode, alnsg, size_factor)

    use haero_constants, only: pi_sixth
    implicit none

    integer,  intent(in) :: imode     !mode number
    real(wp), intent(in) :: alnsg
    real(wp), intent(inout) :: size_factor(:) !size factor

    size_factor(imode) = (pi_sixth)*exp(4.5*(alnsg**2))

    rename_subarea_log 'size_factor=', size_factor(imode), 'with:'
    rename_subarea_log 'imode=', imode
    rename_subarea_log 'alnsg=', alnsg

  end subroutine compute_size_factor

  ! compute volume to number ratio for a mode
  pure function compute_vol_to_num_ratio(imode, alnsg, diameter) result(v2n)

    use haero_constants, only: pi_sixth
    implicit none

    ! Parameters
    integer,  intent(in) :: imode
    real(wp), intent(in) :: alnsg
    real(wp), intent(in) :: diameter(:) ![m]

    ! Return value
    real(wp) :: v2n

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
