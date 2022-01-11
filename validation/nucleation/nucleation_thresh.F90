! This driver computes the binary or ternary nucleation rate for the given
! input.
module nucleation_threshold_mod

  use skywalker
  implicit none

contains

  subroutine usage(exe_name)
    character(len=*), intent(in) :: exe_name

    print *, trim(exe_name), ": usage:"
    print *, trim(exe_name), " <input.yaml>"
    stop
  end subroutine

  subroutine compute_threshold(ensemble)

    use vehkamaki2002

    type(ensemble_t), intent(in) :: ensemble

    type(input_t)  :: input
    type(output_t) :: output
    real(swp)      :: temp, rel_hum, c_thresh

    do while (ensemble%next(input, output))
      ! Fetch input parameters
      rel_hum = input%get("relative_humidity")
      temp = input%get("temperature")

      ! Compute the threshold of H2SO4 above which nucleation occurs.
      c_thresh = h2so4_nucleation_threshold(temp, rel_hum)

      ! Write the computed nucleation rate.
      call output%set("nucleation_threshold", c_thresh)
    end do
  end subroutine

end module

program nucleation_threshold

  use nucleation_threshold_mod
  use skywalker
  use validation

  character(len=255)      :: exe_name, input_file, output_file
  type(ensemble_result_t) :: load_result
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

  ! Run the ensemble.
  ensemble = load_result%ensemble
  call compute_threshold(ensemble)

  ! Write out a Python module.
  output_file = output_name(input_file)
  print *, trim(exe_name), ": writing ", trim(output_file)
  call ensemble%write(trim(output_file))

  ! Clean up.
  call ensemble%free()
end program
