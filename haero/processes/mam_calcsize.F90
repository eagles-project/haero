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

    !local variables
    integer  :: imode, nmodes, nlevs, nspec, max_nspec

    integer  :: s_spec_ind, e_spec_ind ! species starting and ending index in the population array for a mode
    integer  :: ierr

    real(wp) :: dgncur_a(model%num_levels, model%num_modes) !interstitial particle diameter[m](FIXME: This should be diagnostic variable)
    real(wp) :: dgncur_c(model%num_levels, model%num_modes) !cldborne particle diameter [m](FIXME: This should be diagnostic variable)

    real(wp) :: v2ncur_a(model%num_levels, model%num_modes) !interstitial particle diameter[FIXME: units????]
    real(wp) :: v2ncur_c(model%num_levels, model%num_modes) !cldborne vol2num ratio [FIXME:units???]

    real(wp) :: dryvol_a(model%num_levels), dryvol_c(model%num_levels) !dry volume of a particle[FIXME:units????]

    real(wp), pointer, dimension(:,:) :: q_i    ! interstitial aerosol mix ratios [FIXME:units???]
    real(wp), pointer, dimension(:,:) :: q_c    ! cldborne aerosol mix ratios [FIXME:units???]

    real(wp), allocatable, dimension(:) :: dens !specie density array for each mode

    !parameters
    integer, parameter :: top_lev = 1 ![FIXME: This may not be always true]

    !number of levels
    nlevs = model%num_levels

    !number of modes
    nmodes = model%num_modes

    !maximum number of species in the mode with max species
    max_nspec = maxval(model%num_mode_species(:))

    !interstitial mix ratios
    q_i => prognostics%interstitial_aerosols()

    !cloud-borne mix ratios
    q_c => prognostics%interstitial_aerosols()!prognostics%cloudborne_aerosols() !FIXME: THis should be cloud borne, i am getting an error currently

    !allocate variable to store densities of all species in a mode (use max_nspec to allocate, so as to avoid allocation in a nmodes loop)
    allocate(dens(max_nspec), STAT=ierr)
    if (ierr .ne. 0) then
       print *,'Could not allocate dens array, error code=', ierr
       stop
    endif

    print*,nmodes, nlevs, max_nspec

    !Loop through each mode and find particle diameter
    do imode = 1, nmodes

       !Initialize diameter(dgnum), volume to number ratios(v2ncur) and dry volume (dryvol) for both
       !interstitial and cloudborne aerosols

       !NOTE: In Haero we do not carry default dgnum, so we initialize dgncur_* and v2ncur_* to zero
       !That is why, we do not need to send "list_idx" as an argument
       call set_initial_sz_and_volumes (imode, top_lev, nlevs, dgncur_a, v2ncur_a, dryvol_a)

       !for cloud-borne aerosols
       call set_initial_sz_and_volumes (imode, top_lev, nlevs, dgncur_c, v2ncur_c, dryvol_c)

       !----------------------------------------------------------------------
       !Compute dry volume mixrats (aerosol diameter)
       !Current default: number mmr is prognosed
       !       Algorithm:calculate aerosol diameter from mass, number, and fixed sigmag
       !
       !sigmag ("sigma g") is "geometric standard deviation for aerosol mode"
       !
       !Volume = sum_over_components{ component_mass mixrat / density }
       !----------------------------------------------------------------------

       nspec = model%num_mode_species(imode) !total number of species in mode "imode"

       !find start and end index of species in this mode in the "population" array
       !The indices are same for interstitial and cloudborne species

       print*,'nspec:', nspec
       s_spec_ind = model%population_offsets(imode) ! start index
       if(imode.ne.nmodes) then ! for all modes expect the last mode
          e_spec_ind = model%population_offsets(imode+1) - 1! end index
       else !when imode == nmodes
          e_spec_ind = model%num_populations !if imode==nmodes, end index is the total number of species
       endif

       !capture densities for each specie in this mode
       dens(1:max_nspec) = 0.0_wp !initialize the whole array to zero
       dens(1:nspec) = model%aero_species(imode, :)%density !assign dens till nspec (as nspec can be different from max_nspec)

       call compute_dry_volume(imode, top_lev, nlevs, s_spec_ind, e_spec_ind, dens, q_i, q_c, dryvol_a, dryvol_c)

    enddo
    deallocate(dens)

  end subroutine run


  subroutine set_initial_sz_and_volumes(imode, top_lev, nlevs, dgncur, v2ncur, dryvol)

    !-----------------------------------------------------------------------------
    !Purpose: Set initial defaults for the dry diameter, volume to num
    ! and dry volume
    !
    !Called by: run
    !-----------------------------------------------------------------------------
    implicit none

    !inputs
    integer, intent(in) :: top_lev, nlevs !for model level loop
    integer, intent(in) :: imode   !mode index

    !outputs
    real(wp), intent(out) :: dgncur(:,:) !diameter
    real(wp), intent(out) :: v2ncur(:,:) !volume to number
    real(wp), intent(out) :: dryvol(:)   !dry volume

    !local variables
    integer  :: icol, klev
    real(wp) :: dgnum, sigmag, voltonumb


    do klev = top_lev, nlevs
       dgncur(klev,imode) = 0.0_wp !diameter
       v2ncur(klev,imode) = 0.0_wp !volume to number
       dryvol(klev)       = 0.0_wp !initialize dry vol
    end do

    return

  end subroutine set_initial_sz_and_volumes


  subroutine compute_dry_volume(imode, top_lev, nlevs, s_spec_ind, e_spec_ind, dens, q_i, q_c, dryvol_a, dryvol_c)

    !-----------------------------------------------------------------------------
    !Purpose: Compute initial dry volume based on mmr and specie density
    ! volume = mmr/density
    !
    !Called by: modal_aero_calcsize_sub
    !-----------------------------------------------------------------------------

    implicit none

    !inputs
    integer,  intent(in) :: top_lev, nlevs  !for model level loop
    integer,  intent(in) :: imode         !mode index
    integer,  intent(in) :: s_spec_ind, e_spec_ind !start and end indices of population array for this mode

    real(wp), intent(in) :: dens(:)
    real(wp), intent(in) :: q_i(:,:), q_c(:,:)       !interstitial and cldborne mix ratios [FIXME:units???]

    !in-outs
    real(wp), intent(inout) :: dryvol_a(:)                    ! interstital aerosol dry volume [FIXME:units??]
    real(wp), intent(inout) :: dryvol_c(:)                    ! cloud borne aerosol dry volume [FIXME:units??]

    !local vars
    integer  :: ispec, klev, spec_ind
    real(wp) :: dummwdens !density inverse


    character(len=32) :: spec_name

    spec_ind = 0 !species index for dens array goes from 1 to nspec
    do ispec = s_spec_ind, e_spec_ind

     spec_ind = spec_ind + 1 !increment index for dens array
     ! need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air
     dummwdens = 1.0_wp / dens(spec_ind) !inverse of density
     print*,'spec_ind:', spec_ind,s_spec_ind, e_spec_ind

     !compute dry volume as a function of space (i,k)
     do klev = top_lev, nlevs
        dryvol_a(klev) = dryvol_a(klev) + max(0.0_wp,q_i(ispec,klev))*dummwdens
        dryvol_c(klev) = dryvol_c(klev) + max(0.0_wp,q_c(ispec,klev))*dummwdens
        print*,q_i(ispec,klev), 1.0_wp - 10.0_wp*epsilon(1.0_wp)
     end do
  end do ! nspec loop


  end subroutine compute_dry_volume


  pure function compute_diameter(vol2num) result(diameter)

    use haero_precision, only: wp

    implicit none

    real(wp), intent(in) :: vol2num

    real(wp) :: diameter

    diameter = vol2num + 1

  end function compute_diameter


  subroutine finalize(model)
    implicit none

    ! Arguments
    type(model_t), intent(in) :: model
  end subroutine finalize


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


end module mam_calcsize
