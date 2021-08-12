
module mam_rename

  use haero_precision, only: wp
  use haero, only: model_t, aerosol_species_t, gas_species_t, &
       prognostics_t, atmosphere_t, diagnostics_t, tendencies_t

  implicit none
  private

  ! Upper and lower limits for diameters
  real(wp), save, allocatable :: dgnumlo_aer(:), &
                                 dgnumhi_aer(:), &
                                 dgnum_aer(:),   &
                                 alnsg(:)

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

    call initialize_diameters(model)

    call initialize_ln_of_std_dev(model)

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

    ! Start new variables from original fortran routine
    ! TODO: Find better names for these variables ported from the fortran routine
    real(wp) :: dryvol_t_del, dryvol_t_new
    real(wp) :: dryvol_t_old, dryvol_t_oldaa, dryvol_t_oldbnd
    logical  :: iscldy_subarea        ! true if sub-area is cloudy
    real(wp) :: qnum_cur(model%num_modes)
    integer  :: mtoo_renamexf(model%num_modes)
    integer  :: naer, max_aer, max_mode, ntot_amode

    ! TODO: how is max_aer defined? In the original fortran routine, it's often
    ! defined like: 
    ! integer, parameter :: max_aer = nsoa + npoa + nbc + 8

    ! For the following variables, the original declaration is in the comment
    ! preceding the declaration.

    ! real(wp) :: qaercw_del_grow4rnam(1:max_aer, 1:max_mode)
    real(wp) :: qaercw_del_grow4rnam(1:model%num_modes, 1:model%num_modes)

    ! real(wp), intent(in   ), dimension( 1:max_aer, 1:max_mode ) :: qaer_del_grow4rnam
    real(wp) :: qaer_del_grow4rnam(1:model%num_modes, 1:model%num_modes)

    ! real(wp), intent(inout), dimension( 1:max_aer, 1:max_mode ) :: qaer_cur
    real(wp) :: qaer_cur(1:model%num_modes, 1:model%num_modes)

    ! real(wp), intent(inout), optional, dimension( 1:max_aer, 1:max_mode ) :: qaercw_cur
    real(wp) :: qaercw_cur(1:model%num_modes, 1:model%num_modes)

    ! These lengths were originally ntot_amode
    real(wp) :: dryvol_a(model%num_modes)
    real(wp) :: dryvol_c(model%num_modes)
    real(wp) :: deldryvol_a(model%num_modes)
    real(wp) :: deldryvol_c(model%num_modes)
    ! End newly ported variables

    ! TODO: How should ntot_amode be initialized? It seems to come from the global
    ! aero config in the original code
    ntot_amode = model%num_modes
    ! ---

    nmodes = model%num_modes

    iscldy_subarea = .false.

    ! TODO: Initialize newly ported arrays to 0. How should all these really be
    ! initialized?
    dryvol_t_del    = 0._wp
    dryvol_t_new    = 0._wp
    dryvol_t_old    = 0._wp
    dryvol_t_oldaa  = 0._wp
    dryvol_t_oldbnd = 0._wp

    qnum_cur(:)      = 0
    mtoo_renamexf(:) = 0

    qaercw_del_grow4rnam(:, :) = 0
    qaer_del_grow4rnam(:, :)   = 0
    qaer_cur(:, :)             = 0
    qaercw_cur(:, :)           = 0

    dryvol_c(:)    = 0
    dryvol_a(:)    = 0
    deldryvol_a(:) = 0
    deldryvol_c(:) = 0
    ! ---

    ! TODO: This should not be hardwired here but should be either part of the
    ! metadata or otherwise populated.
    dest_mode_of_mode(:) = [0, 1, 0, 0]

    ! TODO: new parameters needed for _find_renaming_pairs_; extract data from
    ! _model_ and pass to _find_renaming_pairs_. DONT pass model to subroutines.
    call find_renaming_pairs(nmodes, dest_mode_of_mode,        &    ! input
        num_pairs, sz_factor, fmode_dist_tail_fac, v2n_lo_rlx, &    ! output
        v2n_hi_rlx, ln_diameter_tail_fac, diameter_cutoff,     &    ! output
        ln_dia_cutoff, diameter_belowcutoff, dryvol_smallest)       ! output

    ! Compute initial (before growth) aerosol dry volume and also the growth in
    ! dryvolume for both interstitial and cloud-borne (if iscldy_subaera is
    ! true) aerosols of the "src" mode
    call compute_dryvol_change_in_src_mode(ntot_amode, naer, mtoo_renamexf, &              !input
        iscldy_subarea, qaer_cur, qaer_del_grow4rnam, qaercw_cur, qaercw_del_grow4rnam, & !input
        dryvol_a, deldryvol_a, dryvol_c, deldryvol_c)                                     !output

  end subroutine run

subroutine compute_dryvol_change_in_src_mode(nmode, nspec, dest_mode_of_mode, &
    iscldy, qi_vmr, qi_del_growth, qcld_vmr, qcld_del_growth, &
    dryvol_a, deldryvol_a, dryvol_c, deldryvol_c)

  integer,  intent(in):: nmode ! total number of modes
  integer,  intent(in):: nspec !total number of species in a mode
  integer,  intent(in):: dest_mode_of_mode(:) ! destination mode for a mode

  logical,  intent(in) :: iscldy ! true if it is a cloudy cell

  real(wp), intent(in) :: qi_vmr(:,:)           ! mass mixing ratios (mmr) [kmol/kmol]
  real(wp), intent(in) :: qi_del_growth(:,:) !growth in mmr [kmol/kmol]

  real(wp), intent(in), optional :: qcld_vmr(:,:)
  real(wp), intent(in), optional :: qcld_del_growth(:,:)

  !intent-outs
  real(wp), intent(out) :: dryvol_a(:), dryvol_c(:)       !dry volumes (before growth) [m3/kmol-air]
  real(wp), intent(out) :: deldryvol_a(:), deldryvol_c(:) !change in dry volumes [m3/kmol-air]

  integer :: imode
  integer :: dest_mode

  !For each mode, compute the initial (before growth) dryvolume and the growth in dryvolume
  do imode = 1, nmode
    !compute dry volume only for modes participating in inter-modal transfer
    dest_mode = dest_mode_of_mode(imode)
    if (dest_mode <= 0) cycle

    !compute dry volumes (before growth) and its change for interstitial aerosols
    call dryvolume_change(imode, nspec, qi_vmr, qi_del_growth, & !input
      dryvol_a(imode), deldryvol_a(imode)) !output

    if ( iscldy ) then ! if this grid cell has cloud
      !if a grid cell is cloudy, clloud borne quantities has to be present
      if(.not. present(qcld_vmr) .or. .not. present(qcld_del_growth)) then
        print *, 'If a grid cell is cloudy, dryvol_c and deldryvol_c should be present'
        stop 1
      endif
      !compute dry volume (before growth) and its change for cloudborne aerosols
      call dryvolume_change(imode, nspec, qcld_vmr, qcld_del_growth, &!input
            dryvol_c(imode), deldryvol_c(imode)) !output
    end if !iscldy then
  end do

  end subroutine compute_dryvol_change_in_src_mode

  subroutine dryvolume_change (imode, nspec, q_vmr, q_del_growth, &!input
       dryvol, deldryvol) !output

    !intent-ins
    integer,  intent(in) :: imode           !current mode number
    integer,  intent(in) :: nspec           !number of species in the current mode
    real(wp), intent(in) :: q_vmr(:,:)        !volume mixing ratio [kmol/kmol] FIXME: units needs to be reverified
    real(wp), intent(in) :: q_del_growth(:,:) !change (delta) in volume mixing ratio [kmol/kmol]

    !intent-outs
    real(wp), intent(out) :: dryvol, deldryvol !dry volume (before growth) and its grwoth [m3/kmol]

    !local variables
    integer  :: ispec, s_spec_ind, e_spec_ind
    real(wp) :: mass_2_vol(nspec) ! converts specie mass to dry volume !DO NOT PORT, we will construct it during "init"
    real(wp) :: tmp_dryvol, tmp_del_dryvol

    ! Original comment: "Temporary variable name change- do not port"
    ! I've therefore commented out this line. It seems this is copying to/from
    ! a module-level variable that we shouldn't touch in hearo?
    !
    ! mass_2_vol(:) = fac_m2v_aer(:)

    !For each mode, we compute a dry volume by combining (accumulating) mass/density for each specie in that mode.
    !conversion from mass to volume is accomplished by multiplying with precomputed "mass_2_vol" factor

    s_spec_ind = 1     !start specie index for this mode [These will be subroutine args]
    e_spec_ind = nspec !end specie index for this mode

    !initialize tmp accumulators
    tmp_dryvol     = 0.0_wp !dry volume accumulator
    tmp_del_dryvol = 0.0_wp !dry volume growth(change) accumulator

    !Notes on mass_2_vol factor:Units:[m3/kmol-s]; where kmol-s is the amount of a specie "s"
    ! This factor is obtained by  (molecular_weight/density) of a specie. That is,
    ! [ (kg/kmol-s) / (kg/m3) ]; where molecular_weight has units [kg/kmol-s] and density units are [kg/m3]
    ! which results in the units of m3/kmol-s

    do ispec = s_spec_ind, e_spec_ind
      !Multiply by mass_2_vol[m3/kmol-s] to convert q_vmr[kmol-s/kmol-air]) to volume units[m3/kmol-air]
      tmp_dryvol     = tmp_dryvol     + q_vmr(ispec,imode)*mass_2_vol(ispec)        !compute current dryvolume
      !accumulate the "grwoth" in volume units as well
      tmp_del_dryvol = tmp_del_dryvol + q_del_growth(ispec,imode)*mass_2_vol(ispec) !compute dryvolume growth
    end do

    dryvol    = tmp_dryvol-tmp_del_dryvol ! This is dry volume before the growth
    deldryvol = tmp_del_dryvol          ! change in dry volume due to growth

  end subroutine dryvolume_change

  subroutine initialize_diameters(model)

    implicit none

    type(model_t), intent(in)    :: model

    ! Holds most recent error code
    integer :: ierr

    ierr = 0

    allocate(dgnumlo_aer(model%num_modes), stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not allocate dgnumlo_aer with length ', model%num_modes
      stop ierr
    endif

    allocate(dgnumhi_aer(model%num_modes), stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not allocate dgnumhi_aer with length ', model%num_modes
      stop ierr
    endif

    allocate(dgnum_aer(model%num_modes), stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not allocate dgnum_aer with length ', model%num_modes
      stop ierr
    endif

    ! Initialize min and max diameters
    dgnumlo_aer(:) = model%modes(:)%min_diameter
    dgnumhi_aer(:) = model%modes(:)%max_diameter

    ! Initialize this to the minimum diameter for now.
    ! TODO: this will be updated with the correct calculation later.
    dgnum_aer(:) = dgnumlo_aer(:)

  end subroutine initialize_diameters

  subroutine finalize_diameters()
    implicit none

    ! Holds most recent error code
    integer :: ierr

    ierr = 0

    deallocate(dgnumlo_aer, stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not deallocate dgnumlo_aer'
      stop ierr
    endif

    deallocate(dgnumhi_aer, stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not deallocate dgnumhi_aer'
      stop ierr
    endif

    deallocate(dgnum_aer, stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not deallocate dgnum_aer'
      stop ierr
    endif

  end subroutine finalize_diameters

  subroutine initialize_ln_of_std_dev(model)
    implicit none

    ! Parameters
    type(model_t), intent(in)     :: model

    integer :: ierr

    ierr = 0

    allocate(alnsg(model%num_modes), stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not allocate alnsg with length ', model%num_modes
      stop ierr
    endif

    alnsg(:) = log(model%modes(:)%mean_std_dev)

  end subroutine initialize_ln_of_std_dev

  subroutine finalize_ln_of_std_dev()

    implicit none

    ! Local variables
    integer :: ierr

    ierr = 0

    deallocate(alnsg, stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not deallocate alnsg'
      stop ierr
    endif

  end subroutine finalize_ln_of_std_dev

  subroutine finalize(model)
    implicit none

    ! Arguments
    type(model_t), intent(in) :: model

    call finalize_diameters()

    call finalize_ln_of_std_dev()

  end subroutine finalize

  subroutine find_renaming_pairs (nmodes, dest_mode_of_mode,  & ! input
       num_pairs, sz_factor, fmode_dist_tail_fac, v2n_lo_rlx, & ! output
       v2n_hi_rlx, ln_diameter_tail_fac, diameter_cutoff,     & ! output
       ln_dia_cutoff, diameter_belowcutoff, dryvol_smallest)    ! output

    ! --- arguments (intent-ins)

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
    real(wp), parameter :: sqrt_half = sqrt(0.5_wp)
    real(wp), parameter :: frelax = 27.0_wp !(3^3)
    real(wp), parameter :: smallest_dryvol_value = 1.0e-25_wp

    ! number of pairs allowed to do inter-mode particle transfer
    ! (e.g. if we have a pair "mode_1<-->mode_2", mode_1 and mode_2 can participate in
    ! inter-mode aerosol particle transfer where like particles in mode_1 can be
    ! transferred to mode_2 and vice-versa)
    !
    ! Let us assume there are none to start with.
    num_pairs = 0

    ! if there can be no possible pairs, just return
    if (all(dest_mode_of_mode(:)<=0)) then
      return
    endif

    ! Find >=1 pair
    do imode = 1, nmodes

      ! Destination mode for mode _imode_
      dest_mode = dest_mode_of_mode(imode)

      ! if dest_mode is <=0, transfer is not possible for this mode. cycle the
      ! loop for the next mode
      if(dest_mode <= 0) then
        cycle
      end if

      ! Source mode is the current mode (i.e. _imode_)
      src_mode = imode

      ! log of stddev for current mode
      alnsg_for_current_mode = alnsg(src_mode)

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
      ln_diameter_tail_fac(src_mode) = 3.0_wp * (alnsg_for_current_mode**2)

      ! Cut-off (based on geometric mean) for making decision to do inter-mode
      ! transfers

      ! TODO: use dummy values for _dgnum_aer_, or assign to dgnum_low for the
      ! moment. Have to figure out how to compute this. We will extract from
      ! model at some point.
      diameter_cutoff(src_mode) = sqrt(   &
         dgnum_aer(src_mode)*exp(1.5_wp*(alnsg_for_current_mode**2)) *   &
         dgnum_aer(dest_mode)*exp(1.5_wp*(alnsg(dest_mode)**2)) )

      ln_dia_cutoff(src_mode) = log(diameter_cutoff(src_mode)) !log of cutt-off
      diameter_belowcutoff(src_mode) = 0.99_wp*diameter_cutoff(src_mode) !99% of the cutoff

    enddo

  end subroutine find_renaming_pairs


  ! Compute size factor for a mode
  subroutine compute_size_factor(imode, alnsg, size_factor)

    use haero_constants, only: pi_sixth
    implicit none

    integer,  intent(in) :: imode     !mode number
    real(wp), intent(in) :: alnsg
    real(wp), intent(inout) :: size_factor(:) !size factor

    size_factor(imode) = (pi_sixth)*exp(4.5_wp*(alnsg**2))

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


end module mam_rename
