module mam_calcsize

  use haero_precision, only: wp
  use haero, only: model_t, aerosol_species_t, gas_species_t, &
       prognostics_t, atmosphere_t, diagnostics_t, tendencies_t

  implicit none
  ! Module functions
  public :: init, &
            run, &
            finalize, &
            set_integer_param, &
            set_logical_param, &
            set_real_param

  !global parameters
  real(wp), parameter :: third   = 1.0_wp/3.0_wp

contains
  subroutine init(model)

    implicit none

    ! Arguments
    type(model_t), intent(in) :: model
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
    integer  :: imode, nmodes, nlevs, nspec, max_nspec

    integer  :: s_spec_ind, e_spec_ind ! species starting and ending index in the population (q_i and q_c) arrays for a mode
    integer  :: ierr, klev

    real(wp) :: v2nmin, v2nmax, dgnmin, dgnmax, cmn_factor !max and mins for diameter and volume to number ratios

    real(wp) :: drv_a, drv_c
    real(wp) :: num_a, num_c

    real(wp) :: dgncur_a(model%num_levels, model%num_modes) !interstitial particle diameter[m](FIXME: This should be diagnostic variable)
    real(wp) :: dgncur_c(model%num_levels, model%num_modes) !cldborne particle diameter [m](FIXME: This should be diagnostic variable)

    real(wp) :: v2ncur_a(model%num_levels, model%num_modes) !interstitial particle diameter[FIXME: units????]
    real(wp) :: v2ncur_c(model%num_levels, model%num_modes) !cldborne vol2num ratio [FIXME:units???]

    real(wp) :: dryvol_a(model%num_levels), dryvol_c(model%num_levels) !dry volume of a particle[m3/kg(of air)]

    real(wp), pointer, dimension(:,:) :: q_i    ! interstitial aerosol mix ratios [kg/kg(of air)]]
    real(wp), pointer, dimension(:,:) :: q_c    ! cldborne aerosol mix ratios [kg/kg(of air)]

    real(wp), pointer, dimension(:,:) :: n_i    ! interstitial aerosol number mixing ratios [#/kg(of air)]
    real(wp), pointer, dimension(:,:) :: n_c    ! cldborne aerosol number mixing ratios [#/kg(of air)]

    real(wp), allocatable, dimension(:) :: density !specie density array for each mode [kg/m3]

    !parameters
    integer, parameter :: top_lev = 1 ![FIXME: This may not be always true]

    !------------------------------
    !Extract relevant data
    !------------------------------

    !number of levels
    nlevs = model%num_levels

    !number of modes
    nmodes = model%num_modes

    !number of species in the mode with the max species
    max_nspec = maxval(model%num_mode_species(:))

    !interstitial mass and number mixing ratios
    q_i => prognostics%interstitial_aerosols()
    n_i => prognostics%interstitial_num_concs()

    !cloud-borne mass and number mixing ratios
    q_c => prognostics%cloudborne_aerosols()
    n_c => prognostics%cloudborne_num_concs()

    !allocate variable to store densities of all species in a mode
    !(use max_nspec to allocate, so as to avoid allocation in the nmodes loop below)
    allocate(density(max_nspec), STAT=ierr)
    if (ierr .ne. 0) then
       print *,'Could not allocate density array, error code=', ierr
       stop
    endif

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
       s_spec_ind = model%population_offsets(imode)       !start index
       e_spec_ind = model%population_offsets(imode+1) - 1 !end index of species for all (modes expect the last mode)

       if(imode.eq.nmodes) then ! for last mode
          e_spec_ind = model%num_populations !if imode==nmodes, end index is the total number of species
       endif

       nspec = model%num_mode_species(imode) !total number of species in mode "imode"

       !capture densities for each specie in this mode
       density(1:max_nspec) = huge(density) !initialize the whole array to a huge value [FIXME: NaN would be better than huge]
       density(1:nspec) = model%aero_species(imode, :)%density !assign density till nspec (as nspec can be different for each mode)

       !FIXME:Density needs to be fixed, the values are not right as they are all for the coarse mode currently, probably for simplification!

       call compute_dry_volume(imode, top_lev, nlevs, s_spec_ind, e_spec_ind, density, q_i, q_c, dryvol_a, dryvol_c)

       !compute upper and lower limits for volume to num (v2n) ratios and diameters (dgn)
       v2nmin = model%modes(imode)%min_vol_to_num_ratio()
       v2nmax = model%modes(imode)%max_vol_to_num_ratio()
       dgnmin = model%modes(imode)%min_diameter
       dgnmax = model%modes(imode)%max_diameter

       !FIXME: compute relaxed counterparts as well

       do klev = top_lev, nlevs

          drv_a = dryvol_a(klev)
          num_a = max( 0.0_wp, n_i(klev,imode))

          drv_c = dryvol_c(klev)
          num_c = max( 0.0_wp, n_c(klev,imode))

          !compute a common factor
          cmn_factor = exp(4.5_wp*log(model%modes(imode)%mean_std_dev)**2.0_wp)*pi_sixth

          !FIXME: size adjustment is done here based on volume to num ratios

          !FIXME: in (or better done after) the following update_diameter_and_vol2num calls, we need to update mmr as well
          !but we are currently skipping that update. That update wil require additional arguments

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
    deallocate(density)

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

  !----------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------

  subroutine finalize(model)
    implicit none

    ! Arguments
    type(model_t), intent(in) :: model
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
