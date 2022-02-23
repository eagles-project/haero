! This program tests MAM's nucleation parameterizations using Skywalker


module uptkrates_mod
  use skywalker
  implicit none
contains
  subroutine fatal_error(message, line)
    character(len=*), intent(in) :: message
    integer :: line

    print *, message, line
    stop
  end subroutine

  function approx_equal(x, y) result(equal)
    real(swp), intent(in) :: x, y
    logical :: equal

    if (abs(x - y) < 1e-6) then
      equal = .true.
    else
      equal = .false.
    end if
  end function
end module uptkrates_mod
!--------------------------------------------------------------

! This macro halts the program if the predicate x isn't true.
#define assert(x) if (.not. (x)) call fatal_error("Assertion failed at line", __LINE__)
!--------------------------------------------------------------

!==============================================================
program skywalker_gasaerexch_gas_aer_uptkrates_1box1gas

  use uptkrates_mod

  ! skywalter

  use skywalker, only: ensemble_t, ensemble_result_t, settings_t, input_t, &
                       output_t, SW_SUCCESS, load_ensemble, input_array_result_t

  ! Subroutine to be tested
  use mam_gasaerexch,  only: gas_aer_uptkrates_1box1gas           

  implicit none

  ! Constants/parameters/flags that are input to the subroutine we will test
  integer, parameter :: ntot_amode = 4

  character(len=255)  :: exe_name, input_file, output_file
  integer             :: suffix, slash
  type(ensemble_result_t) :: load_result
  type(ensemble_t)    :: ensemble
  type(settings_t)    :: settings
  type(input_t)       :: input
  type(output_t)      :: output
  type(input_array_result_t):: values

  ! Tmp variables
  logical  :: l_condense_to_mode(ntot_amode)
  real(swp) :: temp             ! air temperature (K)
  real(swp) :: pmid             ! air pressure at model levels (Pa)
  real(swp) :: pstd             ! 101325 Pa
  real(swp) :: mw_gas           ! molecular weight of gas
  real(swp) :: mw_air           ! molecular weight of air
  real(swp) :: vol_molar_gas    !
  real(swp) :: vol_molar_air    !

  real(swp)  :: accom           ! accomodation coefficient (--)

  real(swp)  :: r_universal     ! universal gas constant
  real(swp)  :: pi              ! pi
  real(swp)  :: beta_inp        ! quadrature parameter (--)

  integer   :: n_mode          ! number of modes
  real(swp)  :: dgncur_awet(ntot_amode) ! mode-median wet diameter of number distribution (m)
  real(swp)  :: lnsg(ntot_amode)        ! ln( sigmag )  (--)

  real(swp)  :: aernum(ntot_amode)   ! aerosol number mixing ratio
  real(swp)  :: aircon             ! air molar concentration (kmol/m3)

  real(swp)  :: uptkaer(ntot_amode)  ! gas-to-aerosol mass transfer rates (1/s)
  real(swp)  :: test_uptkaer(ntot_amode)  ! used in the case the input has a solution 
  logical   :: has_solution
  integer   :: i

  !-------------------------------------------------------
  ! Read command line arguments--exit if they aren't given
  !-------------------------------------------------------
  call get_command_argument(0, exe_name)
  if (command_argument_count() >= 1) then
    call get_command_argument(1, input_file)
  else
    print *, trim(exe_name), ": No input file given!"
    stop
  end if

  ! Generate an output file name based on the name of the input file.
  suffix = index(trim(input_file), ".yaml")
  if (suffix > 0) then
    ! Look for the last slash in the filename.
    slash = index(trim(input_file), "/", back=.true.)
    if (slash == 0) slash = 1
    output_file = input_file(slash+1:suffix-1) // ".py"
  else
    print *, trim(exe_name), ": Invalid input filename (no .yaml suffix found): ", &
             trim(input_file)
    stop
  end if

  !-------------------------------------------------------
  ! Load the ensemble. Any error encountered is fatal.
  !-------------------------------------------------------
  print *, trim(exe_name), ": Loading ensemble from ", trim(input_file)
  load_result = load_ensemble(trim(input_file), "settings")
  ! Make sure everything is as it should be.
  if (load_result%error_code /= SW_SUCCESS) then
    print *, "skywalker_user_test_gasaerexch: ", trim(load_result%error_message)
    stop
  end if

  ! check settings
  settings = load_result%settings

  ! Parse parameters in the process section

  assert(settings%has ("name"))
  assert(trim(settings%get("name")) == "gas_aer_uptkrates_1box1gas")

  ! Check ensemble information

  ensemble = load_result%ensemble
  ! assert(ensemble%size >= 10)

  do while (ensemble%next(input, output))
    assert(input%has("temp"))
    assert(input%has_array("dgncur_awet"))
    assert(input%has_array("lnsg"))
    assert(input%has_array("aernum"))
    has_solution = input%has_array("uptkaer")
  end do

  !-------------------------------------------------------
  ! Process input, do calculations, and prepare output
  !-------------------------------------------------------
  n_mode          =     4
  accom           =     0.65000000000000002
  aircon          =     4.4055781358372036E-002
  beta_inp        =     0.0000000000000000
  pi              =     3.1415926535897931
  r_universal     =  8314.4675910000005
  mw_gas          =    98.078400000000002
  mw_air          =    28.966000000000001
  pmid            = 100000.00000000000
  pstd            = 101325.00000000000
  vol_molar_gas   =    42.880000000000003
  vol_molar_air   =    20.100000000000001
  do i = 1, ntot_amode
    l_condense_to_mode(i) = .true.
  end do
  do while (ensemble%next(input, output))
     ! Parse input 
     values      = input%get_array_param("dgncur_awet")
     dgncur_awet = values%values
     values      = input%get_array_param("lnsg")
     lnsg        = values%values
     values      = input%get_array_param("aernum")
     aernum      = values%values
     temp = input%get("temp")
     if (has_solution) then 
       values       = input%get_array_param("uptkaer")
       test_uptkaer = values%values
     end if
     call gas_aer_uptkrates_1box1gas(l_condense_to_mode,    &
                                     temp,                  &
                                     pmid,                  &
                                     pstd,                  &
                                     mw_gas,                &
                                     mw_air,                &
                                     vol_molar_gas,         &
                                     vol_molar_air,         &
                                     accom,                 &
                                     r_universal,           &
                                     pi,                    &
                                     beta_inp,              & 
                                     n_mode,                &
                                     dgncur_awet,           &
                                     lnsg,                  &
                                     aernum,                &
                                     aircon,                &
                                     uptkaer )
     ! Process output
     call output%set_array("uptkaer", uptkaer)
     if (has_solution) then 
       do i = 1, n_mode
          assert(approx_equal(test_uptkaer(i), uptkaer(i)))
       end do
     end if
  end do

  ! Now we write out a Python module containing the output data.
  print *, trim(exe_name), ": Writing output to ", trim(output_file)
  call ensemble%write(trim(output_file))

  ! Clean up.
  call ensemble%free()

  ! If we got here, the execution was successfull.
  print *, "PASS"

end program
