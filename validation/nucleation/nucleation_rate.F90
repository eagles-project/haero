! This driver computes the binary or ternary nucleation rate for the given
! input.
module nucleation_rate_mod

  use skywalker
  implicit none

contains

  subroutine usage(exe_name)
    character(len=*), intent(in) :: exe_name

    print *, trim(exe_name), ": usage:"
    print *, trim(exe_name), " <input.yaml>"
    stop
  end subroutine

  subroutine run_vehkamaki2002(ensemble, pbl_method)

    use vehkamaki2002

    type(ensemble_t), intent(in) :: ensemble
    integer,          intent(in) :: pbl_method

    type(input_t)  :: input
    type(output_t) :: output
    real(swp)      :: c_h2so4, rel_hum, temp, x_crit, J

    do while (ensemble%next(input, output))
      ! Fetch input parameters
      c_h2so4 = input%get("c_h2so4")
      rel_hum = input%get("relative_humidity")
      temp = input%get("temperature")

      ! Compute the mole fraction of H2SO4 in a critical cluster, and from it
      ! the nucleation rate.
      x_crit = h2so4_critical_mole_fraction(c_h2so4, temp, rel_hum)
      J = nucleation_rate(c_h2so4, temp, rel_hum, x_crit)

      ! Write the computed nucleation rate.
      call output%set("nucleation_rate", J)
    end do
  end subroutine

  subroutine run_merikanto2007(ensemble, pbl_method)

    use merikanto2007

    type(ensemble_t), intent(in) :: ensemble
    integer,          intent(in) :: pbl_method

    type(input_t)  :: input
    type(output_t) :: output
    real(swp)      :: c_h2so4, xi_nh3, rel_hum, temp, log_J, J

    do while (ensemble%next(input, output))
      ! Fetch inputs
      c_h2so4 = input%get("c_h2so4")
      xi_nh3 = input%get("xi_nh3")
      rel_hum = input%get("relative_humidity")
      temp = input%get("temperature")

      ! Compute the nucleation rate.
      log_J = log_nucleation_rate(temp, rel_hum, c_h2so4, xi_nh3)
      J = exp(log_J)

      ! Write the computed nucleation rate.
      call output%set("nucleation_rate", J)
    end do
  end subroutine
end module

program nucleation_rate

  use nucleation_rate_mod
  use skywalker
  use validation

  character(len=255)      :: exe_name, input_file, output_file, val
  type(ensemble_result_t) :: load_result
  type(settings_t)        :: settings
  integer                 :: nuc_method, pbl_method
  type(ensemble_t)        :: ensemble

  ! Get the name of the executable.
  call get_command_argument(0, exe_name)

  if (command_argument_count() < 1) then
    call usage(exe_name)
  end if

  call get_command_argument(1, input_file)
  print *, trim(exe_name), ": reading ", trim(input_file)

  ! Load the ensemble. Any error encountered is fatal.
  load_result = load_ensemble(trim(input_file), "haero")

  ! Figure out settings for binary/ternary nucleation and planetary boundary
  ! layer treatment
  settings = load_result%settings
  if (settings%has("nucleation_method")) then
    val = settings%get("nucleation_method")
    read(val, '(i2)') nuc_method
  else
    nuc_method = 2
  end if

  if ((nuc_method /= 2) .and. (nuc_method /= 3)) then
    print *, "Invalid nucleation method: ", nuc_method
    stop
  end if

  if (settings%has("pbl_method")) then
    val = settings%get("pbl_method")
    read(val, '(i2)') pbl_method
  else
    pbl_method = 0
  end if
  if ((pbl_method < 0) .or. (pbl_method > 2)) then
    print *, "Invalid planetary boundary layer method: ", pbl_method
    stop
  end if

  ! Run the ensemble.
  ensemble = load_result%ensemble
  if (nuc_method == 2) then ! binary nucleation
    call run_vehkamaki2002(ensemble, pbl_method)
  else ! ternary nucleation
    call run_merikanto2007(ensemble, pbl_method)
  end if

  ! Write out a Python module.
  output_file = output_name(input_file)
  print *, trim(exe_name), ": writing ", trim(output_file)
  call ensemble%write(trim(output_file))

  ! Clean up.
  call ensemble%free()
end program
