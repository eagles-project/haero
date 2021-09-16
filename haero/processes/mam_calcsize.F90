module mam_calcsize


  !TODO:
  !1.  do_adjust and do_aitacc_transfer should be set somewhere else and should be an input to this process
  !2. We might need to revist "close_to_one" variable for single precision???
  use haero_precision, only: wp
  use haero, only: model_t, aerosol_species_t, gas_species_t, &
       prognostics_t, atmosphere_t, diagnostics_t, tendencies_t

  implicit none
  private
  ! Module functions
  public :: init, &
            run, &
            finalize, &
            set_integer_param, &
            set_logical_param, &
            set_real_param

  !FIXME: The following parameters should set somewhere else
  logical, parameter :: do_adjust = .true.
  logical, parameter :: do_aitacc_transfer = .true.

  !global parameters
  real(wp), parameter :: third   = 1.0_wp/3.0_wp

  !module-level private variables used accross this module (initialized in "init" routine)
  integer, save :: nlevs      !number of levels
  integer, save :: nmodes     !number of modes
  integer, save :: max_nspec  !number of species in the mode with the max species
  integer, save :: num_populations !total number of species
  integer, save :: aitken_idx !index of aitken mode
  integer, save :: accum_idx  !index of accumulation mode

  integer, save, allocatable :: population_offsets(:) ! start and end indices of species in a mode innpopulation arrays
  integer, save, allocatable :: num_mode_species(:)   ! number of species in a mode
  integer, save, allocatable :: spec_density(:,:)     ! density of species

  real(wp), save, allocatable :: v2nmin_nmodes(:)
  real(wp), save, allocatable :: v2nmax_nmodes(:)
  real(wp), save, allocatable :: dgnmin_nmodes(:)
  real(wp), save, allocatable :: dgnmax_nmodes(:)

contains
  subroutine init(model)

    implicit none

    ! Arguments
    type(model_t), intent(in) :: model

    !local
    integer :: ierr !error code
    integer :: imode

    !initialize module-level variables
    nlevs     = model%num_levels !number of levels
    nmodes    = model%num_modes  !number of modes
    num_populations = model%num_populations !total number of species

    allocate(population_offsets(nmodes), stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not allocate population_offsets with length ', nmodes
      stop
    endif

    !An array that stores species' start and end indices in a specie spopulation index
    population_offsets(1:nmodes) = model%population_offsets(1:nmodes)

    allocate(num_mode_species(nmodes), stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not allocate num_mode_species with length ', nmodes
      stop
    endif
    num_mode_species(1:nmodes) = model%num_mode_species(1:nmodes) !total number of species in a mode

    max_nspec = maxval(num_mode_species(1:nmodes))!number of species in the mode with the max species


    allocate(spec_density(nmodes,max_nspec), stat=ierr)
    if (ierr .ne. 0) then
       print *, 'Could not allocate spec_density with length ', nmodes
       stop
    endif

    spec_density(1:nmodes,1:max_nspec) = model%aero_species(1:nmodes,1:max_nspec)%density !specie density TODO: units

    allocate(v2nmin_nmodes(nmodes), v2nmax_nmodes(nmodes), dgnmin_nmodes(nmodes), dgnmax_nmodes(nmodes), stat=ierr)
    if (ierr .ne. 0) then
       print *, 'Could not allocate v2nmin_nmodes, v2nmax_nmodes, dgnmin_nmodes, dgnmax_nmodes with length ', nmodes
       stop
    endif

    do imode = 1, nmodes
       v2nmin_nmodes(imode) = model%modes(imode)%min_vol_to_num_ratio()
       v2nmax_nmodes(imode) = model%modes(imode)%max_vol_to_num_ratio()
       dgnmin_nmodes(imode) = model%modes(imode)%min_diameter
       dgnmax_nmodes(imode) = model%modes(imode)%max_diameter
    enddo

    aitken_idx = model%mode_index("aitken")
    accum_idx  = model%mode_index("accum")

  end subroutine init

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  subroutine run(model, t, dt, prognostics, atmosphere, diagnostics, tendencies)

    use haero_constants, only: pi_sixth

    implicit none

    ! Arguments
    type(model_t), intent(in)         :: model
    real(wp), value, intent(in)       :: t
    real(wp), value, intent(in)       :: dt
    type(prognostics_t), intent(in)   :: prognostics
    type(atmosphere_t), intent(in)    :: atmosphere
    type(diagnostics_t), intent(in)   :: diagnostics
    type(tendencies_t), intent(inout) :: tendencies

    !local variables
    integer  :: imode, nspec

    integer  :: s_spec_ind, e_spec_ind ! species starting and ending index in the population (q_i and q_c) arrays for a mode
    integer  :: ierr, klev

    real(wp) :: v2nmin, v2nmax, dgnmin, dgnmax, cmn_factor !max and mins for diameter and volume to number ratios
    real(wp) :: v2nminrl, v2nmaxrl !relaxed counterparts of max and mins for volume to number ratios

    real(wp) :: drv_a, drv_c !dry volume [FIXME: unit????]
    real(wp) :: num_a, num_c ![#/kg(of air)]

    real(wp) :: dgncur_a(nlevs, nmodes) !interstitial particle diameter[m](FIXME: This should be diagnostic variable)
    real(wp) :: dgncur_c(nlevs, nmodes) !cldborne particle diameter [m](FIXME: This should be diagnostic variable)

    real(wp) :: v2ncur_a(nlevs, nmodes) !interstitial particle diameter[FIXME: units????]
    real(wp) :: v2ncur_c(nlevs, nmodes) !cldborne vol2num ratio [FIXME:units???]

    real(wp) :: dryvol_a(nlevs), dryvol_c(nlevs) !dry volume of a particle[m3/kg(of air)]
    real(wp) :: init_num_a, init_num_c !initial number mixing ratios entering this process [#/kg(of air)]
    real(wp) :: interstitial_tend, cloudborne_tend ![#/kg(of air)/s]
    real(wp):: density(max_nspec) !specie density array for each mode [kg/m3]

    real(wp), pointer, dimension(:,:) :: q_i    ! interstitial aerosol mix ratios [kg/kg(of air)]]
    real(wp), pointer, dimension(:,:) :: q_c    ! cldborne aerosol mix ratios [kg/kg(of air)]

    real(wp), pointer, dimension(:,:) :: n_i    ! interstitial aerosol number mixing ratios [#/kg(of air)]
    real(wp), pointer, dimension(:,:) :: n_c    ! cldborne aerosol number mixing ratios [#/kg(of air)]

    real(wp), pointer, dimension(:,:) :: dnidt    ! interstitial aerosol number mixing ratios tendency [#/kg(of air)/s]
    real(wp), pointer, dimension(:,:) :: dncdt    ! cldborne     aerosol number mixing ratios tendency [#/kg(of air)/s]

    !parameters
    integer, parameter :: top_lev = 1 ![FIXME: This may not be always true]

    !------------------------------
    !Extract relevant data
    !------------------------------

    !interstitial mass and number mixing ratios
    q_i => prognostics%interstitial_aerosols()
    n_i => prognostics%interstitial_num_mix_ratios()

    !cloud-borne mass and number mixing ratios
    q_c => prognostics%cloud_aerosols()
    n_c => prognostics%cloud_num_mix_ratios()

    !tendencies for interstitial number mixing ratios
    dnidt => tendencies%interstitial_num_concs()

    !tendencies for cloud-borne number mixing ratios
    dncdt => tendencies%cloudborne_num_concs()

    !Loop through each mode and find particle diameter
    do imode = 1, nmodes

       !Initialize diameter(dgnum), volume to number ratios(v2ncur) and dry volume (dryvol) for both
       !interstitial and cloudborne aerosols

       !NOTE: In Haero we do not carry default dgnum, so we initialize dgncur_* and v2ncur_* to zero
       !That is why, we do not need to send "list_idx" as an argument. This call can be removed
       !as there is no need to initialize these fields to zero, they will be updated eventually
       !with the valid values
       call set_initial_sz_and_volumes (imode, top_lev, nlevs, dgncur_a, v2ncur_a, dryvol_a) !for interstitial aerosols
       call set_initial_sz_and_volumes (imode, top_lev, nlevs, dgncur_c, v2ncur_c, dryvol_c) !for cloud-borne aerosols

       !----------------------------------------------------------------------
       !Algorithm to compute dry aerosol diameter:
       !calculate aerosol diameter volume, volume is computed from mass and density
       !----------------------------------------------------------------------

       !find start and end index of species in this mode in the "population" array
       !The indices are same for interstitial and cloudborne species
       s_spec_ind = population_offsets(imode)       !start index
       e_spec_ind = population_offsets(imode+1) - 1 !end index of species for all (modes expect the last mode)

       if(imode.eq.nmodes) then ! for last mode
          e_spec_ind = num_populations !if imode==nmodes, end index is the total number of species
       endif

       nspec = num_mode_species(imode) !total number of species in mode "imode"

       !capture densities for each specie in this mode
       density(1:max_nspec) = huge(density) !initialize the whole array to a huge value [FIXME: NaN would be better than huge]
       density(1:nspec) = spec_density(imode, 1:nspec) !assign density till nspec (as nspec can be different for each mode)

       !FIXME:Density needs to be fixed, the values are not right as they are all for the coarse mode currently, probably for simplification!

       call compute_dry_volume(imode, top_lev, nlevs, s_spec_ind, e_spec_ind, density, q_i, q_c, dryvol_a, dryvol_c)

       !compute upper and lower limits for volume to num (v2n) ratios and diameters (dgn)
       v2nmin = v2nmin_nmodes(imode)
       v2nmax = v2nmax_nmodes(imode)
       dgnmin = dgnmin_nmodes(imode)
       dgnmax = dgnmax_nmodes(imode)

       !Get relaxed limits for volume_to_num
       !(we use relaxed limits for aerosol number "adjustment" calculations via "adjust_num_sizes" subroutine.
       !Note: The relaxed limits will be artifically inflated (or deflated) for the aitken and accumulation modes
       !if "do_aitacc_transfer" flag is true to effectively shut-off aerosol number "adjustment" calculations
       !for these modes because we do the explicit transfer (via "aitken_accum_exchange" subroutine) from one
       !mode to another instead of adjustments for these modes)
       call get_relaxed_v2n_limits(do_aitacc_transfer, & !inputs
            imode == aitken_idx, imode == accum_idx, & !inputs
            v2nmin, v2nmax, v2nminrl, v2nmaxrl)!outputs (NOTE: v2nmin and v2nmax are only updated for aitken and accumulation modes)

       do klev = top_lev, nlevs

          !interstital aerosols
          drv_a = dryvol_a(klev) !dry volume
          init_num_a = n_i(klev,imode) !inital value of num_a for this level and mode
          num_a = max( 0.0_wp, init_num_a) ! Make it non-negative

          !cloudborne aerosols
          drv_c = dryvol_c(klev) !dry volume
          init_num_c = n_c(klev,imode) !inital value of num_c for this level and mode
          num_c = max( 0.0_wp, init_num_c) ! Make it non-negative

          !compute a common factor
          cmn_factor = exp(4.5_wp*log(model%modes(imode)%mean_std_dev)**2.0_wp)*pi_sixth

          !FIXME: size adjustment is done here based on volume to num ratios
          if (do_adjust) then


             !-----------------------------------------------------------------
             ! Do number adjustment for interstitial and activated particles
             !-----------------------------------------------------------------
             !Adjustments that are applied over time-scale deltat (model time step in seconds):
             !(1) make numbers non-negative or
             !(2) make numbers zero when volume is zero
             !
             !
             !Adjustments that are applied over time-scale of a day (in seconds)
             !(3) bring numbers to within specified bounds
             ![Adjustment details are explained in the process]
             !-----------------------------------------------------------------


             !number tendencies to be updated by adjust_num_sizes subroutine
             interstitial_tend = dnidt(klev,imode)
             cloudborne_tend   = dncdt(klev,imode)

             call adjust_num_sizes (drv_a, drv_c, init_num_a, init_num_c, dt, v2nmin, v2nmax, v2nminrl, v2nmaxrl, & !input
                  num_a, num_c, interstitial_tend, cloudborne_tend )!output



          endif !do_adjust



          !FIXME: in (or better done after) the following update_diameter_and_vol2num calls, we need to update mmr as well
          !but we are currently skipping that update. That update will require additional arguments

          !update diameters and volume to num ratios for interstitial aerosols
          call update_diameter_and_vol2num(klev, imode, drv_a, num_a, &
               v2nmin, v2nmax, dgnmin, dgnmax, cmn_factor, &
               dgncur_a, v2ncur_a)

          !update diameters and volume to num ratios for cloudborne aerosols
          call update_diameter_and_vol2num(klev, imode, drv_c, num_c, &
               v2nmin, v2nmax, dgnmin, dgnmax, cmn_factor, &
               dgncur_c, v2ncur_c)

       enddo! klev

    enddo! imodes

  end subroutine run

  !-----------------------------------------------------------------------------
  !Set initial defaults for the dry diameter, volume to num
  ! and dry volume
  !-----------------------------------------------------------------------------
  subroutine set_initial_sz_and_volumes(imode, top_lev, nlevs, dgncur, v2ncur, dryvol)

    implicit none

    !inputs
    integer, intent(in) :: top_lev, nlevs !for model level loop
    integer, intent(in) :: imode   !mode index

    !outputs
    real(wp), intent(out) :: dgncur(:,:) !diameter [m]
    real(wp), intent(out) :: v2ncur(:,:) !volume to number [FIXME:units???]
    real(wp), intent(out) :: dryvol(:)   !dry volume [m3/kg(of air)]

    !local variables
    integer  :: klev

    do klev = top_lev, nlevs
       dgncur(klev,imode) = 0.0_wp !diameter
       v2ncur(klev,imode) = 0.0_wp !volume to number
       dryvol(klev)       = 0.0_wp !initialize dry vol
    end do

    return

  end subroutine set_initial_sz_and_volumes


  !-----------------------------------------------------------------------------
  !Compute initial dry volume based on bulk mass mixing ratio (mmr) and specie density
  ! volume = mmr/density
  !-----------------------------------------------------------------------------
  subroutine compute_dry_volume(imode, top_lev, nlevs, s_spec_ind, e_spec_ind, density, &
       q_i, q_c, dryvol_a, dryvol_c)

    implicit none

    !inputs
    integer,  intent(in) :: top_lev, nlevs         !for model level loop
    integer,  intent(in) :: imode                  !mode index
    integer,  intent(in) :: s_spec_ind, e_spec_ind !start and end indices of population array for this mode

    real(wp), intent(in) :: density(:)             !species density [kg/m3]
    real(wp), intent(in) :: q_i(:,:), q_c(:,:)     !interstitial and cldborne mix ratios [kg/kg(of air)]

    !in-outs
    real(wp), intent(inout) :: dryvol_a(:)         ! interstital aerosol dry volume [m3/kg(of air)]
    real(wp), intent(inout) :: dryvol_c(:)         ! cloud borne aerosol dry volume [m3/kg(of air)]

    !local vars
    integer  :: ispec, klev, density_ind
    real(wp) :: inv_density !density inverse [m3/kg]

    do ispec = s_spec_ind, e_spec_ind
       density_ind = ispec - s_spec_ind + 1 !density array index goes from 1 to nspec

       inv_density = 1.0_wp / density(density_ind) !inverse of density

       !compute dry volume as a function of space (i,k)
       do klev = top_lev, nlevs
          ! volume is mass*inv_density = [kg/kg(of air)] * [1/(kg/m3)] = [m3/kg(of air)]
          dryvol_a(klev) = dryvol_a(klev) + max(0.0_wp,q_i(klev,ispec))*inv_density
          dryvol_c(klev) = dryvol_c(klev) + max(0.0_wp,q_c(klev,ispec))*inv_density
       end do
    end do

  end subroutine compute_dry_volume

  !----------------------------------------------------------------------------------------
  ! Compute particle diameter and volume to number ratios using dry bulk volume (drv)
  !----------------------------------------------------------------------------------------

  subroutine update_diameter_and_vol2num(klev, imode, drv, num, &
       v2nmin, v2nmax, dgnmin, dgnmax, cmn_factor, &
       dgncur, v2ncur)

    implicit none

    !inputs
    integer,  intent(in) :: klev, imode
    real(wp), intent(in) :: drv ![m3/kg(of air)]
    real(wp), intent(in) :: num ![#/kg(of air)]
    real(wp), intent(in) :: v2nmin, v2nmax, dgnmin, dgnmax, cmn_factor

    !output
    real(wp), intent(inout) :: dgncur(:,:), v2ncur(:,:)

    if (drv > 0.0_wp) then
       if (num <= drv*v2nmin) then
          dgncur(klev,imode) = dgnmin !set to minimum diameter for this mode
          v2ncur(klev,imode) = v2nmin !set to minimum vol2num ratio for this mode
       else if (num >= drv*v2nmax) then
          dgncur(klev,imode) = dgnmax !set to maximum diameter for this mode
          v2ncur(klev,imode) = v2nmax !set to maximum vol2num ratio for this mode
       else
          dgncur(klev,imode) = (drv/(cmn_factor*num))**third !compute diameter based on dry volume (drv)
          v2ncur(klev,imode) = num/drv
       end if
    end if

    return
  end subroutine update_diameter_and_vol2num

  subroutine get_relaxed_v2n_limits( do_aitacc_transfer, & !inputs
       is_aitken_mode, is_accum_mode, & ! inputs
       v2nmin, v2nmax, v2nminrl, v2nmaxrl ) !outputs

    implicit none

    !intent-ins
    logical,  intent(in) :: do_aitacc_transfer !flag to control whether to transfer aerosols from one mode to another
    logical,  intent(in) :: is_aitken_mode     !true if this mode is aitken mode
    logical,  intent(in) :: is_accum_mode      !true if this mode is accumulation mode

    !intent-(in)outs
    real(wp), intent(inout) :: v2nmin, v2nmax     !volume_to_num min/max ratios
    real(wp), intent(out)   :: v2nminrl, v2nmaxrl ! relaxed counterparts of volume_to_num min/max ratios

    !local
    !(relaxation factor is currently assumed to be a factor of 3 in diameter
    !which makes it 3**3=27 for volume)
    !i.e. dgnumlo_relaxed = dgnumlo/3 and dgnumhi_relaxed = dgnumhi*3; therefore we use
    !3**3=27 as a relaxation factor for volume
    real(wp), parameter :: relax_factor = 27.0_wp !relax_factor=3**3=27

    !factor to artifically inflate or deflate v2nmin and v2nmax
    real(wp), parameter :: szadj_block_fac = 1.0e6_wp

    !default relaxation:
    v2nminrl = v2nmin/relax_factor
    v2nmaxrl = v2nmax*relax_factor

    !if do_aitacc_transfer is turned on, we will do the ait<->acc tranfer separately in
    !aitken_accum_exchange subroutine, so we are effectively turning OFF the size adjustment for these
    !two modes here by artifically inflating (or deflating) v2min and v2nmax using "szadj_block_fac"
    !and then computing v2minrl and v2nmaxrl based on newly computed v2min and v2nmax.
    if ( do_aitacc_transfer ) then
       !for aitken mode, divide v2nmin by 1.0e6 to effectively turn off the
       !         adjustment when number is too small (size is too big)
       if (is_aitken_mode) v2nmin = v2nmin/szadj_block_fac

       !for accumulation, multiply v2nmax by 1.0e6 to effectively turn off the
       !         adjustment when number is too big (size is too small)
       if (is_accum_mode) v2nmax = v2nmax*szadj_block_fac

       !Also change the v2nmaxrl/v2nminrl so that
       !the interstitial<-->activated number adjustment is effectively turned off
       v2nminrl = v2nmin/relax_factor
       v2nmaxrl = v2nmax*relax_factor
    end if

    return
  end subroutine get_relaxed_v2n_limits

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------
  subroutine adjust_num_sizes (drv_a, drv_c, init_num_a, init_num_c, dt, v2nmin, v2nmax, v2nminrl, v2nmaxrl, & !input
       num_a, num_c, dqdt, dqqcwdt )!output

    implicit none

    !intent-ins
    real(wp), intent(in) :: drv_a, drv_c      !dry volumes [TODO:units]
    real(wp), intent(in) :: init_num_a, init_num_c    !initial number mixing ratios [TODO:units]
    real(wp), intent(in) :: dt                !time step [s]
    real(wp), intent(in) :: v2nmin, v2nmax    !volume to number min and max[TODO:units]
    real(wp), intent(in) :: v2nminrl, v2nmaxrl!volume to number "relaxed" min and max[TODO:units]

    !intent-outs
    real(wp), intent(out):: num_a, num_c  !final number  mixing ratios after size adjument
    real(wp), intent(out):: dqdt, dqqcwdt ! number mixing ratio tendencies

    !local variables
    real(wp) :: numbnd, num_a_stp1, num_c_stp1, num_a_stp2, num_c_stp2, num_a_stp3,num_c_stp3
    real(wp) :: delnum_a_stp2, delnum_c_stp2, delnum_a_stp3,delnum_c_stp3
    real(wp) :: total_drv, total_num ![TODO: units???]
    real(wp) :: min_number_bound, max_number_bound
    real(wp) :: dtinv, adj_tscale, adj_tscale_inv, frac_adj_in_dt, delnum_t3

    real(wp), parameter :: close_to_one = 1.0_wp + 1.0e-15_wp
    real(wp), parameter :: seconds_in_a_day = 86400.0_wp

    !-----------------------------------------------------------------------------------------------
    !The logic behind the number adjustment is described in detail in the "else"
    !section of the following "if" condition. In short, here is the logic:
    !
    !We accomplish number adjustments in 3 stesp:
    !1. Ensure that number mixing ratios are either zero or positive to begin with. If both of them
    !   are zero (or less), we make them zero and update tendencies accordingly (logic in the first
    !   "if" block")
    !2. In this step, we use "relaxed" bounds for bringing number mixing ratios in their bounds. This
    !   is accomplished in three sub-steps [(a), (b) and (c)] described in "Step 2" below.
    !3. In this step, we use the actual bounds for bringing number mixing ratios in their bounds. This
    !   is also accomplished in three sub-steps [(a), (b) and (c)] described in "Step 3" below.
    !
    !TODO: I have a strong feeling that all this logic can be encapsulated in a subroutine
    !which can be called for all parts ("elseif" and else" blocks) of the following "if" condition
    !-----------------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------------
    !If the number mixing ratio in a mode is out of mode's min/max range, we re-balance interstitial
    !and cloud borne aerosols such that the number mixing ratio comes within the range. Time step for
    !such an operation is assumed to be one day (in seconds). That is, it is assumed that number
    !mixing ratios will be within range in a day. "adj_tscale" represents that time scale
    !-----------------------------------------------------------------------------------------------

    adj_tscale     = max( seconds_in_a_day, dt ) !time scale for number adjustment
    adj_tscale_inv = 1.0_wp/(adj_tscale*close_to_one)  !inverse of the adjustment time scale
    frac_adj_in_dt = max( 0.0_wp, min( 1.0_wp, dt*adj_tscale_inv ) ) !fraction of adj_tscale covered in the current time step "dt"

    !inverse of time step
    dtinv = 1.0_wp/(dt*close_to_one)

    ! If both interstitial (drv_a) and cloud borne (drv_c) dry volumes are zero (or less)
    ! adjust numbers(num_a and num_c respectively) for both of them to be zero for this mode and level
    if ((drv_a <= 0.0_wp) .and. (drv_c <= 0.0_wp)) then

       num_a   = 0.0_wp
       num_c   = 0.0_wp
       dqdt    = update_num_tends(num_a, init_num_a, dtinv)
       dqqcwdt = update_num_tends(num_c, init_num_c, dtinv)

    else if (drv_c <= 0.0_wp) then
       ! if cloud borne dry volume (drv_c) is zero(or less), the interstitial number/volume == total/combined
       ! apply step 1 and 3, but skip the relaxed adjustment (step 2, see below)
       num_c = 0.0_wp
       numbnd = min_max_bounded(drv_a, v2nmin, v2nmax, num_a)
       num_a  = num_a + (numbnd - num_a)*frac_adj_in_dt

    else if (drv_a <= 0.0_wp) then
       ! interstitial volume is zero, treat similar to above
       num_a = 0.0_wp
       numbnd = min_max_bounded(drv_c, v2nmin, v2nmax, num_c)
       num_c  = num_c + (numbnd - num_c)*frac_adj_in_dt

    else
       !The number adjustment is done in 3 steps:
       !Step 1: assumes that num_a and num_c are non-negative (nothing to be done here)
       !------
       num_a_stp1 = num_a
       num_c_stp1 = num_c

       !Step 2 [Apply relaxed bounds] has 3 parts (a), (b) and (c)
       !Step 2: (a)Apply relaxed bounds to bound num_a and num_c within "relaxed" bounds.
       !------
       numbnd = min_max_bounded(drv_a, v2nminrl, v2nmaxrl, num_a_stp1) !bounded to relaxed min and max

       !------:(b)Ideally, num_* should be in range. If they are not, we assume that
       !        they will reach their maximum (or minimum)for this mode within a day (time scale).
       !        We then compute how much num_* will change in a time step by multiplying the difference
       !        between num_* and its maximum(or minimum) with "frac_adj_in_dt".
       delnum_a_stp2 = (numbnd - num_a_stp1)*frac_adj_in_dt
       num_a_stp2 = num_a_stp1 + delnum_a_stp2 !change in num_a in one time step

       numbnd = min_max_bounded(drv_c, v2nminrl, v2nmaxrl, num_c_stp1) !bounded to relaxed min and max
       delnum_c_stp2 = (numbnd - num_c_stp1)*frac_adj_in_dt
       num_c_stp2 = num_c_stp1 + delnum_c_stp2 !change in num_a in one time step


       !------:(c)We now also need to balance num_* incase only one among the interstitial or cloud-
       !        borne is changing. If interstitial stayed the same (i.e. it is within range)
       !        but cloud-borne is predicted to reach its maximum(or minimum), we modify
       !        interstitial number (num_a), so as to accomodate change in the cloud-borne aerosols
       !        (and vice-versa). We try to balance these by moving the num_* in the opposite
       !        direction as much as possible to conserve num_a + num_c (such that num_a+num_c
       !        stays close to its original value)

       if ((delnum_a_stp2 == 0.0_wp) .and. (delnum_c_stp2 /= 0.0_wp)) then
          num_a_stp2 = min_max_bounded(drv_a, v2nminrl, v2nmaxrl, num_a_stp1-delnum_c_stp2)
       else if ((delnum_a_stp2 /= 0.0_wp) .and. (delnum_c_stp2 == 0.0_wp)) then
          num_c_stp2 = min_max_bounded(drv_c, v2nminrl, v2nmaxrl, num_c_stp1-delnum_a_stp2)
       end if


       !Step3[apply stricter bounds] has 3 parts (a), (b) and (c)
       !Step 3:(a) compute combined total of num_a and num_c
       total_drv = drv_a + drv_c
       total_num = num_a_stp2 + num_c_stp2

       !-----:(b) We now compute amount of num_* to change if total_num
       !          is out of range. If total_num is within range, we don't do anything (i.e.
       !          delnuma3 and delnum_c_stp3 remain zero)
       delnum_a_stp3 = 0.0_wp
       delnum_c_stp3 = 0.0_wp

       min_number_bound = total_drv*v2nmin !"total_drv*v2nmin" represents minimum number for this mode
       max_number_bound = total_drv*v2nmax !"total_drv*v2nmxn" represents maximum number for this mode

       if (total_num < min_number_bound) then
          delnum_t3 = (min_number_bound - total_num)*frac_adj_in_dt!change in total_num in one time step

          !Now we need to decide how to distribute "delnum"(change in number) for num_a and num_c
          if ((num_a_stp2 < drv_a*v2nmin) .and. (num_c_stp2 < drv_c*v2nmin)) then
             !if both num_a and num_c are less than the lower bound
             !distribute "delnum" using weighted ratios
             delnum_a_stp3 = delnum_t3*(num_a_stp2/total_num)
             delnum_c_stp3 = delnum_t3*(num_c_stp2/total_num)

          else if (num_c_stp2 < drv_c*v2nmin) then
             !if only num_c is less than lower bound, assign total change to num_c
             delnum_c_stp3 = delnum_t3
          else if (num_a_stp2 < drv_a*v2nmin) then
             !if only num_a is less than lower bound, assign total change to num_a
             delnum_a_stp3 = delnum_t3
          end if

       else if (total_num > max_number_bound) then
          delnum_t3 = (max_number_bound - total_num)*frac_adj_in_dt !change in total_num in one time step

          !Now we need to decide how to distribute "delnum"(change in number) for num_a and num_c
          if ((num_a_stp2 > drv_a*v2nmax) .and. (num_c_stp2 > drv_c*v2nmax)) then
             !if both num_a and num_c are more than the upper bound
             !distribute "delnum" using weighted ratios
             delnum_a_stp3 = delnum_t3*(num_a_stp2/total_num)
             delnum_c_stp3 = delnum_t3*(num_c_stp2/total_num)
          else if (num_c_stp2 > drv_c*v2nmax) then
             !if only num_c is more than the upper bound, assign total change to num_c
             delnum_c_stp3 = delnum_t3
          else if (num_a_stp2 > drv_a*v2nmax) then
             !if only num_a is more than the upper bound, assign total change to num_a
             delnum_a_stp3 = delnum_t3
          end if
       end if

       !update num_a and num_c
       num_a = num_a_stp2 + delnum_a_stp3
       num_c = num_c_stp2 + delnum_c_stp3
    end if

    !update tendencies
    dqdt    = update_num_tends(num_a, init_num_a, dtinv)
    dqqcwdt = update_num_tends(num_c, init_num_c, dtinv)

  end subroutine adjust_num_sizes

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  pure function update_num_tends(num, num0, dtinv) result (function_return)

    real(wp) , intent(in) :: num, num0, dtinv

    real(wp) :: function_return

    function_return = (num - num0)*dtinv

    !return function_return

  end function update_num_tends

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  pure function min_max_bounded(drv, v2nmin, v2nmax, num) result (function_return)

    real(wp), intent(in) :: drv, v2nmin, v2nmax, num

    real(wp) :: function_return

    function_return = max( drv*v2nmin, min( drv*v2nmax, num ) )
    !return function_return

  end function min_max_bounded

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  subroutine finalize(model)
    implicit none

    ! Arguments
    type(model_t), intent(in) :: model

    !local
    integer :: ierr

    !deallocate arrays

    deallocate(population_offsets, stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not deallocate population_offsets'
      stop
    endif


    !deallocate()
  end subroutine finalize

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  subroutine set_integer_param(name, val)
    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: val

    ! No integers to set!
  end subroutine set_integer_param

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  subroutine set_logical_param(name, val)
    implicit none

    character(len=*), intent(in) :: name
    logical, intent(in) :: val

    ! No logicals to set!
  end subroutine set_logical_param

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  subroutine set_real_param(name, val)
    implicit none

    character(len=*), intent(in) :: name
    real(wp), intent(in) :: val

    ! No reals to set!
  end subroutine set_real_param


end module mam_calcsize
