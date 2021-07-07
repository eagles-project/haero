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
