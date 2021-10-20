module mam_calcsize
  !TODO:
  !1. do_adjust, do_aitacc_transfer, N_DIAG, list_idx and update_mmr should be set somewhere else and should be an input to this process
  !2. We might need to revist "close_to_one" variable for single precision???
  !3. "top_lev" is asumed to be 1 in this code. It can be different from 1
  !4. Form arrays to pair " to and from modes" in the init routine
  !5. E3SM has "dotend" logicals which are assigned in the interface, do we need those?
  !   If yes, we need to add that capability.
  !6. We need nominal dgnumsso that the code can compute geometric mean for ait<-->accum transfer
  !7. We need state%pdel for computing surface fluxes(qsrflx) but we are skipping it for now.
  !8. Define noxf_acc2ait and do we need its length to be equalt to population array?
  !9  Remove unused variables
  !10. improve var names....remove dum, add units and dimensions
  !11. deallocate everything which ia allocated
  !12. NDIAG should be set to max number of radiation diagnostics the code is allowed to handle
  !13. Check if do_adjust and do_adjust_allowed are in sync, otherwise die
  !14. Check if do_aitacc_transfer and do_aitacc_transfer_allowed are in sync, otherwise die



  use haero_precision, only: wp
  use haero, only: model_t, aerosol_species_t, gas_species_t, &
       prognostics_t, atmosphere_t, diagnostics_t, tendencies_t, &
       get_strt_end_spec_ind

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
  !Switch to perform number adjustment after dry diameter calculations
  logical, parameter :: do_adjust = .true.

  !Switch to perform aitken<->accumulation species transfer after number adjustment calculations
  logical, parameter :: do_aitacc_transfer = .true.

  !Switch to know whether to compute tendencies for mass mixing ratios (mmr) or not
  logical, parameter :: update_mmr = .true.

  !Max number of radiation diagnostics FIXME: it should be set somewhere in the "init" codes for the model
  integer, parameter :: N_DIAG = 0

  !"list_idx=0" is reservered for the prognostic call
  !FIXME: We are currently supporting only list_idx=0, generalize and test the code with other list_idx values
  integer, parameter :: list_idx = 0


  !Mimic enumerators for aerosol types
  integer, parameter:: inter_aero   = 1 !interstitial aerosols
  integer, parameter:: cld_brn_aero = 2 !cloud borne species

  !Maximum number of aitken-accumulation pairs we can have (one for each diagnostic list).
  !NOTE: "0" is reserved for the prognostic call, so maxpair_csizxf represents just
  !the diagnostic pairs.
  integer, parameter :: maxpair_csizxf = N_DIAG

  !Total number of pairs of aitken-accumulation modes
  !---------------------------------------------------------------------------------
  !Note: For diagnostic calls, users can ask for any diagnostics like rad_diag_1
  !      and rad_diag_3 (e.g. skipping rad_diag_2). Therefore the arrays should
  !      have a length of N_DIAG (unless we define another array which maps info
  !      such as to include info about the missing diagnostics)
  !---------------------------------------------------------------------------------
  integer, parameter :: npair_csizxf = N_DIAG !total number of possible diagnostic calls

  !global parameters
  real(wp), parameter :: third   = 1.0_wp/3.0_wp
  real(wp), parameter :: close_to_one = 1.0_wp + 1.0e-15_wp

  !------------------------------------------------------------------------------------------
  !module-level private variables used accross this module (initialized in "init" routine)
  !------------------------------------------------------------------------------------------

  !Switch which decides whether we are allowed to perform number adjustment after dry diameter calculations
  !we issue an error message when there is a conflict between "do_adjust_allowed" and "do_allowed"
  logical, save :: do_adjust_allowed

  !Switch which decides whether we are allowed to perform aitken<->accumulation species transfer after number adjustment calculations
  !we issue an error message when there is a conflict between "do_aitacc_transfer_allowed" and "do_aitacc_transfer"
  logical, save :: do_aitacc_transfer_allowed(0:npair_csizxf) ! This is an array as it can be different for each radiatio diagnostic call

  integer, save :: nlevs      !number of levels
  integer, save :: nmodes     !number of modes
  integer, save :: max_nspec  !number of species in the mode with the max species
  integer, save :: num_populations !total number of species
  integer, save :: aitken_idx !index of aitken mode
  integer, save :: accum_idx  !index of accumulation mode

  !"modefrm_csizxf" stores mode number "from" which species will be moved
  !to a mode stored in "modetoo_csizxf".
  ![E.g. if modefrm_csizxf(3)=2 and modetoo_csizxf(3)=1, for rad_diag_3 (notice that
  !these arrays are indexed "3" as rad_diag is rad_diag_3), species will be moved
  !from 2nd mode to the 1st mode.]
  integer, save :: modefrm_csizxf(0:maxpair_csizxf)
  integer, save :: modetoo_csizxf(0:maxpair_csizxf)

  integer, save, allocatable :: lspecfrma_csizxf(:), lspectooa_csizxf(:), lspecfrmc_csizxf(:), lspectooc_csizxf(:)

  integer, save, allocatable :: population_offsets(:) ! start and end indices of species in a mode innpopulation arrays
  integer, save, allocatable :: num_mode_species(:)   ! number of species in a mode
  integer, save, allocatable :: spec_density(:,:)     ! density of species

  real(wp), save, allocatable :: v2nmin_nmodes(:) !Minimum value of volume to number for each mode
  real(wp), save, allocatable :: v2nnom_nmodes(:) !Nominal value of volume to number for each mode
  real(wp), save, allocatable :: v2nmax_nmodes(:) !Maximum value of volume to number for each mode
  real(wp), save, allocatable :: dgnmin_nmodes(:) !Minimum diameter value for each mode
  real(wp), save, allocatable :: dgnnom_nmodes(:) !Maximum diameter value for each mode
  real(wp), save, allocatable :: dgnmax_nmodes(:) !Maximum diameter value for each mode
  real(wp), save, allocatable :: cmn_factor_nmodes(:)    !A common factor used in size calculation

contains
  subroutine init(model)

    use haero_constants, only: pi_sixth

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
      stop 1
    endif

    !An array that stores species' start and end indices in a specie spopulation index
    population_offsets(1:nmodes) = model%population_offsets(1:nmodes)

    allocate(num_mode_species(nmodes), stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not allocate num_mode_species with length ', nmodes
      stop 1
    endif
    num_mode_species(1:nmodes) = model%num_mode_species(1:nmodes) !total number of species in a mode

    max_nspec = maxval(num_mode_species(1:nmodes))!number of species in the mode with the max species

    allocate(spec_density(nmodes,max_nspec), stat=ierr)
    if (ierr .ne. 0) then
       print *, 'Could not allocate spec_density with length ', nmodes
       stop 1
    endif

    spec_density(1:nmodes,1:max_nspec) = model%aero_species(1:nmodes,1:max_nspec)%density !specie density TODO: units

    allocate(v2nmin_nmodes(nmodes), v2nnom_nmodes(nmodes), v2nmax_nmodes(nmodes), dgnmin_nmodes(nmodes), dgnnom_nmodes(nmodes), &
         dgnmax_nmodes(nmodes), cmn_factor_nmodes(nmodes), stat=ierr)
    if (ierr .ne. 0) then
       print *, 'Could not allocate v2nmin_nmodes, v2nnom_nmodes, v2nmax_nmodes, dgnmin_nmodes, dgnmax_nmodes, ', &
            'dgnnom_nmodes, cmn_factor_nmodes with length ', nmodes
       stop 1
    endif

    do imode = 1, nmodes
       v2nmin_nmodes(imode)     = model%modes(imode)%min_vol_to_num_ratio()!Minimum value of volume to number for each mode
       v2nnom_nmodes(imode)     = model%modes(imode)%nom_vol_to_num_ratio()!Maximum value of volume to number for each mode
       v2nmax_nmodes(imode)     = model%modes(imode)%max_vol_to_num_ratio()!Maximum value of volume to number for each mode
       dgnmin_nmodes(imode)     = model%modes(imode)%min_diameter          !Minimum diameter value for each mode
       dgnnom_nmodes(imode)     = model%modes(imode)%nom_diameter          !Nominal diameter value for each mode
       dgnmax_nmodes(imode)     = model%modes(imode)%max_diameter          !Maximum diameter value for each mode
       cmn_factor_nmodes(imode) = exp(4.5_wp*log(model%modes(imode)%mean_std_dev)**2.0_wp)*pi_sixth !A common factor
    enddo

    aitken_idx = model%mode_index("aitken")
    accum_idx  = model%mode_index("accum")

    call find_species_mapping(aitken_idx, accum_idx)

  end subroutine init

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  subroutine find_species_mapping(mode_a, mode_b)

    !intent-ins
    integer, intent(in) :: mode_a, mode_b

    !local
    integer :: ilist
    integer :: nait, nacc, iait, iacc, icnt
    logical :: accum_exists, aitken_exists

    !Initialize by assigning a non-existent mode (-1) to
    !mode "from" and mode "to" for each diagnostic list
    !NOTE: "0" is reserved for the prognostic call
    do ilist = 0, npair_csizxf
      modefrm_csizxf(ilist) = -1
      modetoo_csizxf(ilist) = -1
   enddo

   !do_adjust_allowed allows aerosol size  adjustment to be turned on/off
   !set it to true by default
   do_adjust_allowed = .true.

   !------------------------------------------------------------------------------------------
   !Original comment in the code:
   !"do_aitacc_transfer_allowed" allows aitken <--> accum mode transfer to be turned on/off
   !NOTE: it can only be true when aitken & accum modes are both present
   !      and have prognosed number and diagnosed surface/sigmag
   !------------------------------------------------------------------------------------------

   !find out mode index for accum and aitken modes in the prognostic radiation list (rad_climate)
   !These are used to get strings stored for these modes in modename_amode array
   nait = aitken_idx !mode number of aitken mode
   nacc = accum_idx  !mode number of accumulation mode

   !Go through the radiation list and decide on do_aitacc_transfer_allowed value(true/false)
   !(Note:"0" index is reserved for the prognostic call)
   do ilist = 0, npair_csizxf

      !find out accumulation and aitken modes in the radiation list
      !FIXME: This should get aitken and accumulation mode indices for diagnostic lists
      iacc = nacc !FIXME:Hardwired!! we should get this index frm the diagnostic lists
      iait = nait !FIXME:Hardwired!! we should get this index frm the diagnostic lists

      !find out if aitken or accumulation modes exist in the radiation list
      !(a positive value means that the mode exists)
      accum_exists  = ( iacc > 0)
      aitken_exists = ( iait > 0)

      do_aitacc_transfer_allowed(ilist)=.false. !False by default
      !if both aitken and accumulation modes exist, make it True and assign mode "from" and mode "to" variables
      if(accum_exists .and. aitken_exists .and. iacc .ne. iait) then
         do_aitacc_transfer_allowed(ilist)=.true.
         modefrm_csizxf(ilist) = iait
         modetoo_csizxf(ilist) = iacc
      endif
   enddo


   !-------------------------------------------------------------------------------
   !Find mapping between the corresponding species of aitken and accumulation nodes
   !
   ! For example, soa may exist in accumulation and aitken modes. Here were are
   ! finding mapping between their indices so that if we need to move soa from
   ! accumulation to aitken, we can do that using this mapping
   !-------------------------------------------------------------------------------

#if 0

   ! Only find the mapping if atleast one of the do_aitacc_transfer_allowed is true
   if (any(do_aitacc_transfer_allowed(:))) then

      !Go through the radiation list and find "aitken<-->accumulation" mapping for each list
      do ilist = 0, npair_csizxf
         if(do_aitacc_transfer_allowed(ilist)) then
            icnt = 0
            imode_ait = modefrm_csizxf(ilist) !aitken  mode of this list
            imode_acc = modetoo_csizxf(ilist) !accumulation mode of this list

            !--------------------------------------------------------------------------------------
            !find aerosol *number* indices mapping between aitken and accumulation modes in this list
            !--------------------------------------------------------------------------------------

            !For aitken mode
            call rad_cnst_get_info(ilist, imode_ait, num_name=num_name_ait)
            call cnst_get_ind(num_name_ait, ind_ait)

            !For accumulation mode
            call rad_cnst_get_info(ilist, imode_acc, num_name=num_name_acc)
            call cnst_get_ind(num_name_acc, ind_acc)

#endif


  end subroutine find_species_mapping


  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  subroutine run(t, dt, prognostics, atmosphere, diagnostics, tendencies)

    implicit none

    ! Arguments
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
    real(wp) :: inter_num_tend, cld_num_tend   ![#/kg(of air)/s]
    real(wp) :: density(max_nspec) !specie density array for each mode [kg/m3]
    real(wp) :: adj_tscale, adj_tscale_inv !adjustment time scales in [s] and [1/s] respectively

    !Work variables for aitken<-->accumulation transfer sub process
    real(wp) :: drv_a_aitsv(nlevs), drv_a_accsv(nlevs) !saves aitken and accumulation interstitial modes dryvolume
    real(wp) :: drv_c_aitsv(nlevs), drv_c_accsv(nlevs) !saves aitken and accumulation cloudborne modes dryvolume
    real(wp) :: num_a_aitsv(nlevs), num_a_accsv(nlevs) !saves aitken and accumulation interstitial modes num concentrations
    real(wp) :: num_c_aitsv(nlevs), num_c_accsv(nlevs) !saves aitken and accumulation cloudborne modes num concentrations

    real(wp) :: drv_a_sv(nlevs,nmodes), drv_c_sv(nlevs,nmodes)!saves dryvolume for each mode and level
    real(wp) :: num_a_sv(nlevs,nmodes), num_c_sv(nlevs,nmodes)!saves num conc. for each mode and level


    real(wp), pointer, dimension(:,:) :: q_i    ! interstitial aerosol mix ratios [kg/kg(of air)]]
    real(wp), pointer, dimension(:,:) :: q_c    ! cldborne aerosol mix ratios [kg/kg(of air)]

    real(wp), pointer, dimension(:,:) :: n_i    ! interstitial aerosol number mixing ratios [#/kg(of air)]
    real(wp), pointer, dimension(:,:) :: n_c    ! cldborne aerosol number mixing ratios [#/kg(of air)]

    real(wp), pointer, dimension(:,:) :: didt    ! interstitial aerosol mass mixing ratios tendency [kg/kg(of air)/s]
    real(wp), pointer, dimension(:,:) :: dcdt    ! cldborne     aerosol mass mixing ratios tendency [kg/kg(of air)/s]

    real(wp), pointer, dimension(:,:) :: dnidt    ! interstitial aerosol number mixing ratios tendency [#/kg(of air)/s]
    real(wp), pointer, dimension(:,:) :: dncdt    ! cldborne     aerosol number mixing ratios tendency [#/kg(of air)/s]

    !parameters
    integer,  parameter :: top_lev = 1 ![FIXME: This may not be always true]
    real(wp), parameter :: seconds_in_a_day = 86400.0_wp !number of seconds in a day


    !---------------------------------------------------------------------------------------------------
    !"adj_tscale" is the amount of time, in seconds, in a day. This is used as a time scale for the
    !processes which are assumed to be occuring at a day's time scale. We compute its inverse as well
    !in advance.
    !---------------------------------------------------------------------------------------------------
    adj_tscale     = max( seconds_in_a_day, dt )       !A one day time scale (used in number adjustment process)
    adj_tscale_inv = 1.0_wp/(adj_tscale*close_to_one)  !inverse of the adjustment time scale

    !------------------------------
    !Extract relevant data
    !------------------------------

    !interstitial mass and number mixing ratios
    q_i => prognostics%interstitial_aerosols()
    n_i => prognostics%interstitial_num_mix_ratios()

    !cloud-borne mass and number mixing ratios
    q_c => prognostics%cloud_aerosols()
    n_c => prognostics%cloud_num_mix_ratios()

    !tendencies for interstitial mass mixing ratios
    didt => tendencies%interstitial_aerosols()

    !tendencies for interstitial number mixing ratios
    dnidt => tendencies%interstitial_num_mix_ratios()

    !tendencies for cloudborne mass mixing ratios
    dcdt => tendencies%cloud_aerosols()

    !tendencies for cloud-borne number mixing ratios
    dncdt => tendencies%cloud_num_mix_ratios()

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

          !a common factor
          cmn_factor = cmn_factor_nmodes(imode)

          !FIXME: size adjustment is done here based on volume to num ratios
          if (do_adjust) then


             !-----------------------------------------------------------------
             ! Do number adjustment for interstitial and activated particles
             !-----------------------------------------------------------------
             !Adjustments that are applied over time-scale dt (model time step in seconds):
             !(1) make numbers non-negative or
             !(2) make numbers zero when volume is zero
             !
             !
             !Adjustments that are applied over time-scale of a day (in seconds)
             !(3) bring numbers to within specified bounds
             ![Adjustment details are explained in the process]
             !-----------------------------------------------------------------


             !number tendencies to be updated by adjust_num_sizes subroutine
             inter_num_tend = dnidt(klev,imode)
             cld_num_tend   = dncdt(klev,imode)

             !NOTE: Only number tendencies (NOT mass mixing ratios) are updated in adjust_num_sizes
             call adjust_num_sizes (drv_a, drv_c, init_num_a, init_num_c, dt, adj_tscale_inv, & !input
                  v2nmin, v2nmax, v2nminrl, v2nmaxrl, &       !input
                  num_a, num_c, inter_num_tend, cld_num_tend )!output



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

          ! save number concentrations and dry volumes for explicit aitken <--> accum mode transfer
          ! which is the next step in the calcsize process
          if ( do_aitacc_transfer ) then
             if (imode == aitken_idx) then
                drv_a_aitsv(klev) = drv_a
                num_a_aitsv(klev) = num_a
                drv_c_aitsv(klev) = drv_c
                num_c_aitsv(klev) = num_c
             else if (imode == accum_idx) then
                drv_a_accsv(klev) = drv_a
                num_a_accsv(klev) = num_a
                drv_c_accsv(klev) = drv_c
                num_c_accsv(klev) = num_c
             end if
          end if
          drv_a_sv(klev,imode) = drv_a
          num_a_sv(klev,imode) = num_a
          drv_c_sv(klev,imode) = drv_c
          num_c_sv(klev,imode) = num_c

       enddo! klev
    enddo! imodes

    !------------------------------------------------------------------------------
    ! Overall logic for aitken<-->acculation transfer:
    ! ------------------------------------------------
    ! when the aitken mode mean size is too big, the largest
    !    aitken particles are transferred into the accum mode
    !    to reduce the aitken mode mean size
    ! when the accum mode mean size is too small, the smallest
    !    accum particles are transferred into the aitken mode
    !    to increase the accum mode mean size
    !------------------------------------------------------------------------------

    if ( do_aitacc_transfer ) then
       call aitken_accum_exchange( nlevs, top_lev, &
            aitken_idx,  accum_idx, adj_tscale_inv, &
            dt, q_i, q_c, &
            drv_a_aitsv, num_a_aitsv, drv_c_aitsv, num_c_aitsv,     &
            drv_a_accsv,num_a_accsv, drv_c_accsv, num_c_accsv,      &
            dgncur_a, v2ncur_a, dgncur_c, v2ncur_c, &!, dotend, dotendqqcw &
            didt, dcdt &
!update_mmr, &
       !     dt, pdel, state_q, state, &
       !     pbuf, &
       !     drv_a_aitsv, num_a_aitsv, drv_c_aitsv, num_c_aitsv, &
       !     drv_a_accsv,num_a_accsv, drv_c_accsv, num_c_accsv,  &
       !     dgncur_a, v2ncur_a, dgncur_c, v2ncur_c, dotend, dotendqqcw, &
       !     dqdt, dqqcwdt, qsrflx
)
    end if


!TODO: after aitken<->accum transfer, the rest of the code deals with history field output and updating the tendencies

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
  subroutine adjust_num_sizes (drv_a, drv_c, init_num_a, init_num_c, dt, adj_tscale_inv, & ! input
       v2nmin, v2nmax, v2nminrl, v2nmaxrl, & !input
       num_a, num_c, dqdt, dqqcwdt )!output

    implicit none

    !intent-ins
    real(wp), intent(in) :: drv_a, drv_c      !dry volumes [TODO:units]
    real(wp), intent(in) :: init_num_a, init_num_c    !initial number mixing ratios [TODO:units]
    real(wp), intent(in) :: dt                !time step [s]
    real(wp), intent(in) :: adj_tscale_inv    !inverse of "adjustment" time scale [1/s]
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
    real(wp) :: dtinv, frac_adj_in_dt, delnum_t3

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
    !mixing ratios will be within range in a day. "adj_tscale" represents that time scale.
    !"frac_adj_in_dt" represents the fraction of that "adj_tscale" covered in the current
    !time step "dt"
    !-----------------------------------------------------------------------------------------------
    frac_adj_in_dt = max( 0.0_wp, min( 1.0_wp, dt*adj_tscale_inv ) )

    !inverse of time step
    dtinv = 1.0_wp/(dt*close_to_one)

    ! If both interstitial (drv_a) and cloud borne (drv_c) dry volumes are zero (or less)
    ! adjust numbers(num_a and num_c respectively) for both of them to be zero for this mode and level
    if ((drv_a <= 0.0_wp) .and. (drv_c <= 0.0_wp)) then

       num_a   = 0.0_wp
       num_c   = 0.0_wp
       dqdt    = update_num_adj_tends(num_a, init_num_a, dtinv)
       dqqcwdt = update_num_adj_tends(num_c, init_num_c, dtinv)

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
    dqdt    = update_num_adj_tends(num_a, init_num_a, dtinv)
    dqqcwdt = update_num_adj_tends(num_c, init_num_c, dtinv)

  end subroutine adjust_num_sizes

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  pure function update_num_adj_tends(num, num0, dtinv) result (function_return)
    !update number mixing ratio tendencies

    real(wp) , intent(in) :: num, num0, dtinv

    real(wp) :: function_return

    function_return = (num - num0)*dtinv

  end function update_num_adj_tends

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  pure function min_max_bounded(drv, v2nmin, v2nmax, num) result (function_return)
    !Bounds the number mixing ratio in min/max bounds, so that it stays within bounds

    real(wp), intent(in) :: drv, v2nmin, v2nmax, num

    real(wp) :: function_return

    function_return = max( drv*v2nmin, min( drv*v2nmax, num ) )

  end function min_max_bounded


  subroutine  aitken_accum_exchange( nlevs, top_lev, &
       aitken_idx,  accum_idx, adj_tscale_inv, &
       dt, q_i, q_c, &
       drv_a_aitsv, num_a_aitsv, drv_c_aitsv, num_c_aitsv,     &
       drv_a_accsv,num_a_accsv, drv_c_accsv, num_c_accsv,      &
       dgncur_a, v2ncur_a, dgncur_c, v2ncur_c, & !, dotend, dotendqqcw, &
       dqdt, dqqcwdt &!, qsrflx
!ncol, lchnk, list_idx, update_mmr, &
!       dt, pdel, state_q, state, &
!       pbuf, &

!       dgncur_a, v2ncur_a, dgncur_c, v2ncur_c, dotend, dotendqqcw, &
!       dqdt, dqqcwdt, qsrflx
)

    !-----------------------------------------------------------------------------
    !Purpose: Exchange aerosols between aitken and accumulation modes based on new
    ! sizes
    !
    !Called by: modal_aero_calcsize_sub
    !Calls    : endrun
    !
    !Author: Richard Easter (Refactored by Balwinder Singh)
    !-----------------------------------------------------------------------------

    use haero_constants, only: gravity

    implicit none

    !inputs
    integer, intent(in) :: top_lev, nlevs         !for model level loop
    integer, intent(in) :: aitken_idx,  accum_idx !mode number for aitken and accumulation modes
    real(wp), intent(in) :: adj_tscale_inv
    real(wp), intent(in) :: q_i(:,:), q_c(:,:)     !interstitial and cldborne mix ratios [kg/kg(of air)]
    real(wp), intent(in) :: drv_a_aitsv(:), num_a_aitsv(:)
    real(wp), intent(in) :: drv_a_accsv(:), num_a_accsv(:)
    real(wp), intent(in) :: drv_c_aitsv(:), num_c_aitsv(:)
    real(wp), intent(in) :: drv_c_accsv(:), num_c_accsv(:)

    !logical,  intent(in) :: update_mmr
    real(wp), intent(in) :: dt
#if 0
    real(wp), intent(in) :: pdel(:,:)
    real(wp), intent(in) :: state_q(:,:,:) !state_q should only be used for list_idx==0
    type(physics_state), intent(in)    :: state       ! Physics state variables
    type(physics_buffer_desc), pointer :: pbuf(:)     ! physics buffer
#endif
    !outputs
    real(wp), intent(inout) :: dgncur_a(:,:) !TODO: Add comments here!
    real(wp), intent(inout) :: dgncur_c(:,:)
    real(wp), intent(inout) :: v2ncur_a(:,:)
    real(wp), intent(inout) :: v2ncur_c(:,:)
    !logical,  intent(inout) :: dotend(:), dotendqqcw(:)
    real(wp), intent(inout) :: dqdt(:,:), dqqcwdt(:,:)!, qsrflx(:,:,:,:)

    !local
    integer  :: iait, iacc
    integer  :: klev
    integer  :: s_spec_ind, e_spec_ind ! starting and ending index in population array for a mode
    integer  :: ixfer_ait2acc, ixfer_acc2ait
    integer  :: imode, jmode, aer_type, jsrflx, iq
    integer  :: lsfrm, lstoo, ispec_acc, idx

    real(wp) :: num_a, drv_a, num_c, drv_c
    real(wp) :: num_a_acc, num_c_acc
    real(wp) :: drv_a_acc, drv_c_acc
    real(wp) :: pdel_fac, num_t, drv_t
    real(wp) :: drv_a_noxf, drv_c_noxf, drv_t_noxf, dummwdens
    real(wp) :: num_t0, num_t_noxf
    real(wp) :: duma, cmn_factor, dgnum, sigmag
    real(wp) :: voltonumb_ait, voltonumb_acc, v2n_geomean
    real(wp) :: xfercoef_num_ait2acc, xfercoef_vol_ait2acc
    real(wp) :: xfercoef_num_acc2ait, xfercoef_vol_acc2ait
    real(wp) :: xfertend_num(2,2)
#if 0
    real(wp) :: xfercoef
    real(wp) :: xfertend
#endif
    logical  :: no_transfer_acc2ait(max_nspec)
#if 0
    character(len=32)  :: spec_name
    character(len=800) :: err_msg
#endif
    !FIXME: The following will become active once we create pairs in the init routine
    !if (npair_csizxf .le. 0)call endrun('npair_csizxf <= 0'//errmsg(__FILE__,__LINE__))

    !find out accumulation and aitken modes
    iacc = accum_idx
    iait = aitken_idx

    if(aitken_idx <= 0  .or. accum_idx <=0) then
       write(*,*)'Accumulation or the Aitken mode do not exist,', &
            ' Accu mode:',accum_idx,', Aitken mode:',aitken_idx,', a negative or zero mode is', &
            ' a non-existent mode'
       stop 1
    endif

    !FIXME: The following will become active once we create pairs in the init routine
    !if (modefrm_csizxf(list_idx) .ne. iait .or.modetoo_csizxf(list_idx) .ne. iacc) then
    !   write(err_msg,*)'modefrm/too_csizxf are wrong for radiation list:',list_idx,' ',errmsg(__FILE__,__LINE__)
    !   call endrun(trim(err_msg))
    !endif

    !FIXME: Do we need dotend capability?
    ! set dotend() for species that will be transferred
    ! for both mass and number
    !if(update_mmr) then
    !   do iq = 1, nspecfrm_csizxf(list_idx)
    !      lsfrm = lspecfrma_csizxf(iq,list_idx)
    !      lstoo = lspectooa_csizxf(iq,list_idx)
    !      if ((lsfrm > 0) .and. (lstoo > 0)) then
    !         dotend(lsfrm) = .true.
    !         dotend(lstoo) = .true.
    !      end if
          !for cloud borne aerosols
    !      lsfrm = lspecfrmc_csizxf(iq,list_idx)
    !      lstoo = lspectooc_csizxf(iq,list_idx)
    !      if ((lsfrm > 0) .and. (lstoo > 0)) then
    !         dotendqqcw(lsfrm) = .true.
    !         dotendqqcw(lstoo) = .true.
    !      end if
    !   end do
    !endif

    !------------------------------------------------------------------------
    ! Identify accum species cannot be transferred to aitken mode
    !
    ! Accumulation mode have more species than Aitken mode. Therefore, there
    ! will be some species which cannot be transferred from accumulation to
    ! Aitken mode as they don't exist in the Aitken mode
    !------------------------------------------------------------------------

    no_transfer_acc2ait(:) = .true. ! let us assume we are not going to move any species

    !balli-
    !Now compute species which can be transfered
    !Get # of species in accumulation mode
    call get_strt_end_spec_ind(population_offsets, iacc, s_spec_ind, e_spec_ind)

    !FIXME: do we need the following  call??
    !call rad_cnst_get_info(list_idx, iacc, nspec = nspec_acc) !output:nspec_acc

    !FIMXME: Following loop will work only when we map species in the init routine
    !do ispec_acc = 1, nspec_acc     !Go through all species within accumulation mode (only species, not number (e.g. num_a1))
    !   call rad_cnst_get_info(list_idx, iacc, ispec_acc,spec_name=spec_name) !output:spec_name
    !   call cnst_get_ind(spec_name, idx)
    !   do iq = 1, nspecfrm_csizxf(list_idx) !Go through all mapped species (and number) for the accumulation mode
    !      if (lspectooa_csizxf(iq,list_idx) == idx) then !compare idx with mapped species in the accumulation mode
    !         no_transfer_acc2ait(ispec_acc) = .false. ! species which can be tranferred
    !      end if
    !   end do
    !end do

    !------------------------------------------------------------------------
    ! Compute geometrically defined midpoint between aitken and accumulation
    ! modes. This value is used to decide whether or not particles are
    ! transfered between these modes
    !------------------------------------------------------------------------

    voltonumb_ait   = v2nnom_nmodes(iait) !volume to number for aitken mode
    voltonumb_acc   = v2nnom_nmodes(iacc) !volume to number for aitken mode

    ! v2n_geomean is voltonumb at the "geometrically-defined" mid-point
    ! between the aitken and accum modes
    v2n_geomean = sqrt(voltonumb_ait*voltonumb_acc) !


    ! loop over columns and levels
    do  klev = top_lev, nlevs

       !FIXME: we need state%pdel, which is currently not available
       !pdel_fac = pdel(icol,klev)/gravity   ! = rho*dz

       !Compute aitken->accumulation transfer
       call compute_coef_ait_acc_transfer(iacc, v2n_geomean, adj_tscale_inv, drv_a_aitsv(klev), & !input
            drv_c_aitsv (klev), num_a_aitsv(klev), num_c_aitsv(klev), voltonumb_acc,            & !input
            ixfer_ait2acc, xfercoef_num_ait2acc, xfercoef_vol_ait2acc, xfertend_num)              !output


       !----------------------------------------------------------------------------------------
       ! compute accum --> aitken transfer rates
       !
       ! accum may have some species (seasalt, dust, poa etc.) that are
       !    not in aitken mode
       ! so first divide the accum drv & num into not-transferred (noxf) species
       !    and transferred species, and use the transferred-species
       !    portion in what follows
       !----------------------------------------------------------------------------------------

       call compute_coef_acc_ait_transfer(iacc, klev, &
            v2n_geomean, adj_tscale_inv, q_i, q_c, drv_a_accsv(klev), drv_c_accsv(klev), num_a_accsv(klev),      &
            num_c_accsv(klev), no_transfer_acc2ait, voltonumb_ait,                                &
            drv_a_noxf, drv_c_noxf, ixfer_acc2ait, xfercoef_num_acc2ait, &
            xfercoef_vol_acc2ait, xfertend_num)

       ! jump to end-of-loop if no transfer is needed at current klev
       if (ixfer_ait2acc+ixfer_acc2ait > 0) then
          !
          ! compute new dgncur & v2ncur for aitken & accum modes
          !
          ! currently inactive (??? BSINGH: Not sure what this comment refers to...)

          !interstitial species
          duma = (xfertend_num(1,1) - xfertend_num(2,1))*dt    !diff in num from  ait->accum and accum->ait transfer
          num_a     = max( 0.0_wp, num_a_aitsv(klev) - duma ) !num removed/added from aitken mode
          num_a_acc = max( 0.0_wp, num_a_accsv(klev) + duma ) !num added/removed to accumulation mode

          duma = (drv_a_aitsv(klev)*xfercoef_vol_ait2acc -   &
               (drv_a_accsv(klev)-drv_a_noxf)*xfercoef_vol_acc2ait)*dt ! diff in volume transfer fomr ait->accum and accum->ait transfer
          drv_a     = max( 0.0_wp, drv_a_aitsv(klev) - duma ) !drv removed/added from aitken mode
          drv_a_acc = max( 0.0_wp, drv_a_accsv(klev) + duma ) !drv added/removed to accumulation mode

          !cloud borne species
          duma = (xfertend_num(1,2) - xfertend_num(2,2))*dt    !same as above for cloud borne aerosols

          num_c     = max( 0.0_wp, num_c_aitsv(klev) - duma )
          num_c_acc = max( 0.0_wp, num_c_accsv(klev) + duma )
          duma = (drv_c_aitsv(klev)*xfercoef_vol_ait2acc -   &
               (drv_c_accsv(klev)-drv_c_noxf)*xfercoef_vol_acc2ait)*dt
          drv_c     = max( 0.0_wp, drv_c_aitsv(klev) - duma )
          drv_c_acc = max( 0.0_wp, drv_c_accsv(klev) + duma )

          !interstitial species (aitken mode)
          call compute_new_sz_after_transfer(iait, drv_a, num_a, &
               dgncur_a(klev,iait), v2ncur_a(klev,iait))

          !cloud borne species (aitken mode)
          call compute_new_sz_after_transfer(iait, drv_c, num_c, &
               dgncur_c(klev,iait), v2ncur_c(klev,iait))

          num_a = num_a_acc
          drv_a = drv_a_acc
          num_c = num_c_acc
          drv_c = drv_c_acc

          !interstitial species (accumulation mode)
          call compute_new_sz_after_transfer(iacc, drv_a, num_a, &
               dgncur_a(klev,iacc), v2ncur_a(klev,iacc))

          !cloud borne species (accumulation mode)
          call compute_new_sz_after_transfer(iacc, drv_c, num_c, &
               dgncur_c(klev,iacc), v2ncur_c(klev,iacc))

          !------------------------------------------------------------------
          ! compute tendency amounts for aitken <--> accum transfer
          !------------------------------------------------------------------
#if 0
          !ASSUMPTION: "update_mmr" will only be true for the prognostic radiation list(i.e. list_idx=0, "radiation_climate")
          !If list_idx=0, it is okay to get specie mmr from state_q array. Therefore, state_q is used in update_tends_flx calls
          !below
          if(update_mmr) then
             ! jmode=1 does aitken-->accum
             if(ixfer_ait2acc > 0) then
                jmode = 1
                call update_tends_flx( klev, jmode, lspecfrma_csizxf, lspectooa_csizxf, &
                     lspecfrmc_csizxf, lspectooc_csizxf, xfertend_num, xfercoef_vol_ait2acc, q_i, &
                     pdel_fac, dqdt, dqqcwdt)!, qsrflx)
             endif

             !jmode=2 does accum-->aitken
             if(ixfer_acc2ait > 0) then
                jmode = 2
                !suboutine (update_tends_flx) is called but lspectooa and lspecfrma are switched and
                !xfercoef_vol_acc2ait is also used instead xfercoef_vol_ait2acc for accum->aitken transfer
                call update_tends_flx( klev, jmode, lspectooa_csizxf, lspecfrma_csizxf, &
                     lspectooc_csizxf, lspecfrmc_csizxf, xfertend_num, xfercoef_vol_acc2ait, q_i, &
                     pdel_fac, dqdt, dqqcwdt)!, qsrflx)
             endif
          endif !update_mmr
#endif
       end if !ixfer_ait2acc+ixfer_acc2ait > 0
!#endif
    end do !klev

    return
  end subroutine aitken_accum_exchange


  !---------------------------------------------------------------------------------------------

  subroutine compute_coef_ait_acc_transfer(iacc, v2n_geomean, adj_tscale_inv, drv_a_aitsv, &
       drv_c_aitsv, num_a_aitsv, num_c_aitsv,  voltonumb_acc, &
       ixfer_ait2acc, xfercoef_num_ait2acc, xfercoef_vol_ait2acc, xfertend_num)

    !------------------------------------------------------------
    ! Purpose: Computes coefficients for transfer from aitken to accumulation mode
    !
    ! Author: Richard Easter (Refactored by Balwinder Singh)
    !------------------------------------------------------------

    !intent ins
    integer,  intent(in) :: iacc
    real(wp), intent(in) :: v2n_geomean
    real(wp), intent(in) :: adj_tscale_inv
    real(wp), intent(in) :: drv_a_aitsv, drv_c_aitsv
    real(wp), intent(in) :: num_a_aitsv, num_c_aitsv, voltonumb_acc

    !intent outs
    integer,  intent(inout) :: ixfer_ait2acc
    real(wp), intent(inout) :: xfercoef_num_ait2acc, xfercoef_vol_ait2acc
    real(wp), intent(inout) :: xfertend_num(2,2)

    !local
    real(wp) :: drv_t, num_t
    real(wp) :: xferfrac_num_ait2acc, xferfrac_vol_ait2acc

    !initialize
    ixfer_ait2acc        = 0
    xfercoef_num_ait2acc = 0.0_wp
    xfercoef_vol_ait2acc = 0.0_wp
    xfertend_num(:,:)    = 0.0_wp

    ! compute aitken --> accum transfer rates

    drv_t = drv_a_aitsv + drv_c_aitsv
    num_t = num_a_aitsv + num_c_aitsv
    if (drv_t > 0.0_wp) then
       !if num is less than the mean value, we have large particles (keeping volume constant drv_t)
       !which needs to be moved to accumulation mode
       if (num_t < drv_t*v2n_geomean) then
          ixfer_ait2acc = 1
          if (num_t < drv_t*voltonumb_acc) then ! move all particles if number is smaller than the acc mean
             xferfrac_num_ait2acc = 1.0_wp
             xferfrac_vol_ait2acc = 1.0_wp
          else !otherwise scale the transfer
             xferfrac_vol_ait2acc = ((num_t/drv_t) - v2n_geomean)/   &
                  (voltonumb_acc - v2n_geomean)
             xferfrac_num_ait2acc = xferfrac_vol_ait2acc*   &
                  (drv_t*voltonumb_acc/num_t)
             !bound the transfer coefficients between 0 and 1
             if ((xferfrac_num_ait2acc <= 0.0_wp) .or.   &
                  (xferfrac_vol_ait2acc <= 0.0_wp)) then
                xferfrac_num_ait2acc = 0.0_wp
                xferfrac_vol_ait2acc = 0.0_wp
             else if ((xferfrac_num_ait2acc >= 1.0_wp) .or.   &
                  (xferfrac_vol_ait2acc >= 1.0_wp)) then
                xferfrac_num_ait2acc = 1.0_wp
                xferfrac_vol_ait2acc = 1.0_wp
             end if
          end if
          xfercoef_num_ait2acc = xferfrac_num_ait2acc*adj_tscale_inv
          xfercoef_vol_ait2acc = xferfrac_vol_ait2acc*adj_tscale_inv
          xfertend_num(1,1) = num_a_aitsv*xfercoef_num_ait2acc
          xfertend_num(1,2) = num_c_aitsv*xfercoef_num_ait2acc
       end if
    end if

  end subroutine compute_coef_ait_acc_transfer


  !---------------------------------------------------------------------------------------------

  subroutine compute_coef_acc_ait_transfer( iacc, klev, &
       v2n_geomean, adj_tscale_inv, q_i, q_c, drv_a_accsv, drv_c_accsv, num_a_accsv,      &
       num_c_accsv, no_transfer_acc2ait, voltonumb_ait,                                &
       drv_a_noxf, drv_c_noxf, ixfer_acc2ait, xfercoef_num_acc2ait, &
       xfercoef_vol_acc2ait, xfertend_num)

    !intent -ins
    integer,  intent(in) :: iacc, klev
    real(wp), intent(in) :: v2n_geomean
    real(wp), intent(in) :: adj_tscale_inv
    real(wp), intent(in) :: q_i(:,:), q_c(:,:)     !interstitial and cldborne mix ratios [kg/kg(of air)]
    real(wp), intent(in) :: drv_a_accsv, drv_c_accsv
    real(wp), intent(in) :: num_a_accsv, num_c_accsv, voltonumb_ait
    logical,  intent(in) :: no_transfer_acc2ait(:)

    !intent - outs
    integer,  intent(inout) :: ixfer_acc2ait
    real(wp), intent(inout) :: drv_a_noxf, drv_c_noxf
    real(wp), intent(inout) :: xfercoef_num_acc2ait, xfercoef_vol_acc2ait
    real(wp), intent(inout) :: xfertend_num(2,2)

    !local
    integer  :: ipop, ispec, s_spec_ind, e_spec_ind
    real(wp) :: drv_t, num_t, drv_t_noxf, num_t0
    real(wp) :: num_t_noxf
    real(wp) :: invdens
    real(wp) :: xferfrac_num_acc2ait, xferfrac_vol_acc2ait
    real(wp), parameter :: zero_div_fac = 1.0e-37_wp

    ixfer_acc2ait = 0
    xfercoef_num_acc2ait = 0.0_wp
    xfercoef_vol_acc2ait = 0.0_wp

    drv_t = drv_a_accsv + drv_c_accsv
    num_t = num_a_accsv + num_c_accsv
    drv_a_noxf = 0.0_wp
    drv_c_noxf = 0.0_wp
    if (drv_t > 0.0_wp) then
       !if number is larger than the mean, it means we have small particles (keeping volume constant drv_t),
       !we need to move particles to aitken mode
       if (num_t > drv_t*v2n_geomean) then
          !As there may be more species in the accumulation mode which are not present in the aitken mode,
          !we need to compute the num and volume only for the species which can be transferred
          call get_strt_end_spec_ind(population_offsets, iacc, s_spec_ind, e_spec_ind)
          ispec = 0 !running index for arrays which starts with specie index of 1 for this mode
          do ipop = s_spec_ind, e_spec_ind ! ipop is population index for q_i and q_c arrays
             ispec = ispec + 1

             if ( no_transfer_acc2ait(ispec) ) then !species which can't be transferred

                ! need qmass*invdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air
                invdens = 1.0_wp / spec_density(iacc,ispec)
                drv_a_noxf = drv_a_noxf + max(0.0_wp,q_i(klev, ipop))*invdens
                drv_c_noxf = drv_c_noxf + max(0.0_wp,q_c(klev, ipop))*invdens
             end if
          end do
          drv_t_noxf = drv_a_noxf + drv_c_noxf !total volume which can't be moved to the aitken mode
          num_t_noxf = drv_t_noxf*v2nmin_nmodes(iacc) !total number which can't be moved to the aitken mode
          num_t0 = num_t
          num_t = max( 0.0_wp, num_t - num_t_noxf )
          drv_t = max( 0.0_wp, drv_t - drv_t_noxf )
       end if
    end if

    if (drv_t > 0.0_wp) then
       !Find out if we need to transfer based on the new num_t
       if (num_t > drv_t*v2n_geomean) then
          ixfer_acc2ait = 1
          if (num_t > drv_t*voltonumb_ait) then! if number of larger than the aitken mean, move all particles
             xferfrac_num_acc2ait = 1.0_wp
             xferfrac_vol_acc2ait = 1.0_wp
          else ! scale the transfer
             xferfrac_vol_acc2ait = ((num_t/drv_t) - v2n_geomean)/   &
                  (voltonumb_ait - v2n_geomean)
             xferfrac_num_acc2ait = xferfrac_vol_acc2ait*   &
                  (drv_t*voltonumb_ait/num_t)
             !bound the transfer coefficients between 0 and 1
             if ((xferfrac_num_acc2ait <= 0.0_wp) .or.   &
                  (xferfrac_vol_acc2ait <= 0.0_wp)) then
                xferfrac_num_acc2ait = 0.0_wp
                xferfrac_vol_acc2ait = 0.0_wp
             else if ((xferfrac_num_acc2ait >= 1.0_wp) .or.   &
                  (xferfrac_vol_acc2ait >= 1.0_wp)) then
                xferfrac_num_acc2ait = 1.0_wp
                xferfrac_vol_acc2ait = 1.0_wp
             end if
          end if
          xferfrac_num_acc2ait = xferfrac_num_acc2ait*   &
               num_t/max( zero_div_fac, num_t0 )
          xfercoef_num_acc2ait = xferfrac_num_acc2ait*adj_tscale_inv
          xfercoef_vol_acc2ait = xferfrac_vol_acc2ait*adj_tscale_inv
          xfertend_num(2,1) = num_a_accsv*xfercoef_num_acc2ait
          xfertend_num(2,2) = num_c_accsv*xfercoef_num_acc2ait
       end if
    end if

  end subroutine compute_coef_acc_ait_transfer

  !----------------------------------------------------------------------

  subroutine compute_new_sz_after_transfer(imode, drv, num, &
       dgncur, v2ncur)

    implicit none

    !intent-ins
    integer,  intent(in) :: imode
    real(wp), intent(in) :: drv, num

    !intent-outs
    real(wp), intent(inout) :: dgncur, v2ncur

    !local
    real(wp) :: voltonumb, voltonumbhi, voltonumblo

    voltonumbhi = v2nmax_nmodes(imode)
    voltonumblo = v2nmin_nmodes(imode)
    voltonumb   = v2nnom_nmodes(imode)

    if (drv > 0.0_wp) then
       if (num <= drv*voltonumbhi) then
          dgncur = dgnmax_nmodes(imode)
          v2ncur = voltonumbhi
       else if (num >= drv*voltonumblo) then
          dgncur = dgnmin_nmodes(imode)
          v2ncur = voltonumblo
       else
          dgncur = drv/(cmn_factor_nmodes(imode)*num)**third
          v2ncur = num/drv
       end if
    else
       dgncur = dgnnom_nmodes(imode)
       v2ncur = voltonumb
    end if

  end subroutine compute_new_sz_after_transfer

  !------------------------------------------------------------------------------------------------

  subroutine update_tends_flx(klev, jmode, frm_spec_a, to_spec_a, &
       frm_spec_c, to_spec_c, xfertend_num, xfercoef, state_q, &
       pdel_fac, dqdt, dqqcwdt)!, qsrflx)

    implicit none

    !intent - ins
    integer,  intent(in) ::  klev, jmode
    integer,  intent(in) :: frm_spec_a(max_nspec, 0:maxpair_csizxf)
    integer,  intent(in) :: to_spec_a(max_nspec, 0:maxpair_csizxf)
    integer,  intent(in) :: frm_spec_c(max_nspec, 0:maxpair_csizxf)
    integer,  intent(in) :: to_spec_c(max_nspec, 0:maxpair_csizxf)

    real(wp), intent(in) :: xfertend_num(2,2)
    real(wp), intent(in) :: xfercoef
    real(wp), intent(in) :: state_q(:,:)
    real(wp), intent(in) :: pdel_fac

    !intent -inout
    real(wp), intent(inout) :: dqdt(:,:), dqqcwdt(:,:)!, qsrflx(:,:,:,:)

    !local
    integer  :: iq, lsfrm, lstoo, jsrflx
    real(wp) :: xfertend
    real(wp), pointer :: fldcw(:,:)    !specie mmr (cloud borne)

    character(len=800) :: err_msg


    jsrflx = jmode + 2

    !interstiatial species
    iq = 1 !iq = 1 is for num_* species
    lsfrm = frm_spec_a(iq,list_idx)
    lstoo = to_spec_a(iq,list_idx)
    if((lsfrm > 0) .and. (lstoo > 0)) then
       call update_num_tends( klev, jmode, jsrflx, lsfrm, lstoo, inter_aero, pdel_fac, xfertend_num, dqdt)!, qsrflx)
    endif
#if 0
    !FIXME: why loop starts from 2, add a comment?
    do iq = 2, nspecfrm_csizxf(list_idx)
       lsfrm = frm_spec_a(iq,list_idx)
       lstoo = to_spec_a(iq,list_idx)
       if((lsfrm > 0) .and. (lstoo > 0)) then
          xfertend = max(0.0_wp,state_q(klev,lsfrm))*xfercoef
          dqdt(klev,lsfrm) = dqdt(klev,lsfrm) - xfertend
          dqdt(klev,lstoo) = dqdt(klev,lstoo) + xfertend
          !qsrflx(lsfrm,jsrflx,inter_aero) = qsrflx(lsfrm,jsrflx,inter_aero) - xfertend*pdel_fac
          !qsrflx(lstoo,jsrflx,inter_aero) = qsrflx(lstoo,jsrflx,inter_aero) + xfertend*pdel_fac
       endif
    enddo
#endif
    !cloud borne apecies
    iq = 1 !number species
    lsfrm = frm_spec_c(iq,list_idx)
    lstoo = to_spec_c(iq,list_idx)

    if((lsfrm > 0) .and. (lstoo > 0)) then
       call update_num_tends( klev, jmode, jsrflx, lsfrm, lstoo, cld_brn_aero, pdel_fac, xfertend_num, dqqcwdt)!, qsrflx)
    endif
#if 0
    !mass species
    do iq = 2, nspecfrm_csizxf(list_idx)
       lsfrm = frm_spec_c(iq,list_idx)
       lstoo = to_spec_c(iq,list_idx)
       if((lsfrm > 0) .and. (lstoo > 0)) then
          !fldcw => qqcw_get_field(pbuf,lsfrm,lchnk)
          xfertend = max(0.0_wp,fldcw(klev))*xfercoef
          dqqcwdt(klev,lsfrm) = dqqcwdt(klev,lsfrm) - xfertend
          dqqcwdt(klev,lstoo) = dqqcwdt(klev,lstoo) + xfertend
          !qsrflx(lsfrm,jsrflx,cld_brn_aero) = qsrflx(lsfrm,jsrflx,cld_brn_aero) - xfertend*pdel_fac
          !qsrflx(lstoo,jsrflx,cld_brn_aero) = qsrflx(lstoo,jsrflx,cld_brn_aero) + xfertend*pdel_fac
       end if
    enddo
#endif
  end subroutine update_tends_flx

  !---------------------------------------------------------------------------------------------------------------------
  subroutine update_num_tends( klev, jmode, jsrflx, lsfrm, lstoo, aer_type, pdel_fac, xfertend_num, dqdt)!, qsrflx)
    !intent ins
    integer,  intent(in) ::  klev, jmode, jsrflx, lsfrm, lstoo, aer_type
    real(wp), intent(in) :: pdel_fac
    real(wp), intent(in) :: xfertend_num(:,:)

    !intent inouts
    real(wp), intent(inout) :: dqdt(:,:)
    !real(wp), intent(inout) :: qsrflx(:,:,:,:)

    !local
    real(wp) :: xfertend

    xfertend = xfertend_num(jmode,aer_type)
    dqdt(klev,lsfrm) = dqdt(klev,lsfrm) - xfertend
    dqdt(klev,lstoo) = dqdt(klev,lstoo) + xfertend
    !qsrflx(lsfrm,jsrflx,aer_type) = qsrflx(lsfrm,jsrflx,aer_type) - xfertend*pdel_fac
    !qsrflx(lstoo,jsrflx,aer_type) = qsrflx(lstoo,jsrflx,aer_type) + xfertend*pdel_fac

  end subroutine update_num_tends



  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  subroutine finalize()
    implicit none

    !local
    integer :: ierr

    !deallocate arrays

    deallocate(population_offsets, stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not deallocate population_offsets'
      stop 1
    endif

    deallocate(num_mode_species, stat=ierr)
    if (ierr .ne. 0) then
      print *, 'Could not deallocate num_mode_species'
      stop 1
    endif

    deallocate(spec_density, stat=ierr)
    if (ierr .ne. 0) then
       print *, 'Could not deallocate spec_density'
       stop 1
    endif

    deallocate(v2nmin_nmodes, v2nmax_nmodes, dgnmin_nmodes, dgnmax_nmodes, cmn_factor_nmodes, stat=ierr)
    if (ierr .ne. 0) then
       print *, 'Could not deallocate v2nmin_nmodes, v2nmax_nmodes, dgnmin_nmodes, dgnmax_nmodes, ', &
            'cmn_factor_nmodes with length ', nmodes
       stop 1
    endif

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
