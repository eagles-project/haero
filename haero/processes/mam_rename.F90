
module mam_rename

  use haero_precision, only: wp
  use haero, only: model_t, aerosol_species_t, gas_species_t, &
       prognostics_t, atmosphere_t, diagnostics_t, tendencies_t

  implicit none
  private

  integer, save :: max_aer
  integer, save :: num_populations
  integer, save :: naer
  integer, save :: num_modes
  integer, save, allocatable :: population_offsets(:)

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

    ! Locals
    integer :: ierr

    num_modes = model%num_modes
    num_populations = model%num_populations

    allocate(population_offsets(num_modes), stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not allocate population_offsets with length ', num_modes
      stop ierr
    endif
    population_offsets(:) = model%population_offsets(:)

    ! FIXME: naer is the total number of species in a mode, it is hardwired here
    ! FIXME: but it should be computed based on the mode number in a mode loop
    naer = population_offsets(2) - population_offsets(1)

    ! FIXME: max_aer is number of species in the mode with most species, it
    ! FIXME: should be computed dynamically using population_offsets
    max_aer = model%num_modes

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

    ! Contains information about the destination mode for a given mode
    integer  :: dest_mode_of_mode(num_modes)

    ! total number of pairs to be found
    integer  :: num_pairs

    ! precomputed factors to be used later
    real(wp) :: sz_factor(num_modes), fmode_dist_tail_fac(num_modes)

    ! relaxed volume to num high and low ratio limits
    real(wp) :: v2n_lo_rlx(num_modes), v2n_hi_rlx(num_modes)

    ! log of diameter factor for distribution tail
    real(wp) :: ln_diameter_tail_fac(num_modes)

    ! cutoff (threshold) for deciding the do inter-mode transfer
    real(wp) :: diameter_cutoff(num_modes), &
                ln_dia_cutoff(num_modes)

    ! some limiters/factors
    real(wp) :: diameter_belowcutoff(num_modes), &
                dryvol_smallest(num_modes)

    ! true if sub-area is cloudy
    logical  :: iscldy_subarea

    integer  :: max_mode

    real(wp) :: qaercw_del_grow4rnam(max_aer, num_modes)
    real(wp) :: qaer_del_grow4rnam(max_aer, num_modes)

    real(wp), pointer :: q_interstitial(:, :)
    real(wp), pointer :: q_cloudborne(:, :)

    real(wp) :: dryvol_a(num_modes)
    real(wp) :: dryvol_c(num_modes)
    real(wp) :: deldryvol_a(num_modes)
    real(wp) :: deldryvol_c(num_modes)

    iscldy_subarea = .true.

    ! FIXME: diagnostics should be updated
    qaercw_del_grow4rnam(:, :) = 0.0_wp
    qaer_del_grow4rnam(:, :)   = 0.0_wp
    ! ---

    q_interstitial => prognostics%interstitial_aerosols()
    q_cloudborne => prognostics%cloud_aerosols()

    !------------------------------------------------------------------------
    ! Find mapping between different modes, so that we can move aerosol
    ! particles from one mode to another
    !------------------------------------------------------------------------

    ! FIXME: All the arrays in find_renaming_pairs subroutine call should be
    ! initialized to HUGE or NaNs as they are partially populated
      
    ! Find (src->destination) pairs of modes which can participate in inter-mode particle transfer
      
    ! NOTE:dryvol_smallest is a very small volume mixing ratio [m3-spc/kmol-air] (where m3-spc
    ! is meter cubed volume of a specie "spc") used for avoiding overflow.  it corresponds to dp = 1 nm
    ! and number = 1e-5 #/mg-air ~= 1e-5 #/cm3-air

    ! TODO: This should not be hardwired here but should be either part of the
    ! metadata or otherwise populated.
    dest_mode_of_mode(:) = [0, 1, 0, 0]

    call find_renaming_pairs(num_modes, dest_mode_of_mode,        &    ! input
        num_pairs, sz_factor, fmode_dist_tail_fac, v2n_lo_rlx, &    ! output
        v2n_hi_rlx, ln_diameter_tail_fac, diameter_cutoff,     &    ! output
        ln_dia_cutoff, diameter_belowcutoff, dryvol_smallest)       ! output

    ! Compute initial (before growth) aerosol dry volume and also the growth in
    ! dryvolume for both interstitial and cloud-borne (if iscldy_subaera is
    ! true) aerosols of the "src" mode
    call compute_dryvol_change_in_src_mode(num_modes, dest_mode_of_mode, &              !input
        iscldy_subarea, q_interstitial, qaer_del_grow4rnam, q_cloudborne, qaercw_del_grow4rnam, & !input
        dryvol_a, deldryvol_a, dryvol_c, deldryvol_c)                                     !output

  end subroutine run

subroutine compute_dryvol_change_in_src_mode(num_modes, dest_mode_of_mode, &
    iscldy, qi_vmr, qi_del_growth, qcld_vmr, qcld_del_growth, &
    dryvol_a, deldryvol_a, dryvol_c, deldryvol_c)

  integer,  intent(in):: num_modes ! total number of modes
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
  integer :: start_species_index, end_species_index

  !For each mode, compute the initial (before growth) dryvolume and the growth in dryvolume
  do imode = 1, num_modes
    !compute dry volume only for modes participating in inter-modal transfer
    dest_mode = dest_mode_of_mode(imode)
    if (dest_mode <= 0) cycle

    ! find start and end index of species in this mode in the "population" array
    ! The indices are same for interstitial and cloudborne species
    start_species_index = population_offsets(imode)
    end_species_index = population_offsets(imode+1) - 1

    if (imode .eq. num_modes) then
      end_species_index = num_populations
    endif

    !compute dry volumes (before growth) and its change for interstitial aerosols
    call dryvolume_change(imode, qi_vmr, qi_del_growth, & !input
      start_species_index, end_species_index, & ! input
      dryvol_a(imode), deldryvol_a(imode)) !output

    if ( iscldy ) then ! if this grid cell has cloud
      !if a grid cell is cloudy, clloud borne quantities has to be present
      if(.not. present(qcld_vmr) .or. .not. present(qcld_del_growth)) then
        print *, 'If a grid cell is cloudy, dryvol_c and deldryvol_c should be present'
        stop 1
      endif

      !compute dry volume (before growth) and its change for cloudborne aerosols
      call dryvolume_change(imode, qcld_vmr, qcld_del_growth, &!input
            start_species_index, end_species_index, & ! input
            dryvol_c(imode), deldryvol_c(imode)) !output

    end if !iscldy then
  end do

  end subroutine compute_dryvol_change_in_src_mode

  subroutine dryvolume_change (imode, q_vmr, q_del_growth, & ! input
                               start_species_index, end_species_index, & ! input
                               dryvol, deldryvol) ! output

    !intent-ins
    integer,  intent(in) :: imode           !current mode number
    integer,  intent(in) :: start_species_index, end_species_index
    real(wp), intent(in) :: q_vmr(:,:)        !volume mixing ratio [kmol/kmol] FIXME: units needs to be reverified
    real(wp), intent(in) :: q_del_growth(:,:) !change (delta) in volume mixing ratio [kmol/kmol]

    !intent-outs
    real(wp), intent(out) :: dryvol, deldryvol !dry volume (before growth) and its grwoth [m3/kmol]

    !local variables
    integer  :: ispec

    !FIXME: This factor needs to be precomputed
    ! mass_2_vol(:) = fac_m2v_aer(:)
    real(wp) :: mass_2_vol(naer)

    real(wp) :: tmp_dryvol, tmp_del_dryvol

    !For each mode, we compute a dry volume by combining (accumulating) mass/density for each specie in that mode.
    !conversion from mass to volume is accomplished by multiplying with precomputed "mass_2_vol" factor

    !initialize tmp accumulators
    tmp_dryvol     = 0.0_wp !dry volume accumulator
    tmp_del_dryvol = 0.0_wp !dry volume growth(change) accumulator

    !Notes on mass_2_vol factor:Units:[m3/kmol-s]; where kmol-s is the amount of a specie "s"
    ! This factor is obtained by  (molecular_weight/density) of a specie. That is,
    ! [ (kg/kmol-s) / (kg/m3) ]; where molecular_weight has units [kg/kmol-s] and density units are [kg/m3]
    ! which results in the units of m3/kmol-s

    do ispec = start_species_index, end_species_index
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

    allocate(dgnumlo_aer(num_modes), stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not allocate dgnumlo_aer with length ', num_modes
      stop ierr
    endif

    allocate(dgnumhi_aer(num_modes), stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not allocate dgnumhi_aer with length ', num_modes
      stop ierr
    endif

    allocate(dgnum_aer(num_modes), stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not allocate dgnum_aer with length ', num_modes
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

    deallocate(population_offsets, stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not deallocate population_offsets'
      stop ierr
    endif

  end subroutine finalize_diameters

  subroutine initialize_ln_of_std_dev(model)
    implicit none

    ! Parameters
    type(model_t), intent(in)     :: model

    integer :: ierr

    ierr = 0

    allocate(alnsg(num_modes), stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not allocate alnsg with length ', num_modes
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

  subroutine find_renaming_pairs (num_modes, dest_mode_of_mode,  & ! input
       num_pairs, sz_factor, fmode_dist_tail_fac, v2n_lo_rlx, & ! output
       v2n_hi_rlx, ln_diameter_tail_fac, diameter_cutoff,     & ! output
       ln_dia_cutoff, diameter_belowcutoff, dryvol_smallest)    ! output

    ! --- arguments (intent-ins)

    ! total number of modes
    integer,  intent(in) :: num_modes

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
    do imode = 1, num_modes

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
