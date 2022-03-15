!> This module implements the stuff needed to provide namelist support to
!> Haero's box model.

module box_input

  use iso_c_binding
  use haero_precision, only: wp

  implicit none

  private

  public :: BoxInput, read_nml

  type, bind(c) :: BoxInput
    ! time_input
    integer(c_int)     mam_dt, mam_nstep, mam_output_intvl

    ! cntl_input
    integer(c_int)     mdo_gaschem, mdo_cloudchem, mdo_gasaerexch, &
                       mdo_test_siamgs19, mdo_rename, mdo_newnuc, mdo_coag, &
                       mdo_calcsize, newnuc_method_user_choice, &
                       pbl_nuc_wang2008_user_choice
    logical(c_bool)    cluster_growth_arg_write
    character(c_char)  cluster_growth_arg_write_fname(128)
    integer(c_int)     num_unit, frac_unit, gas_unit
    real(wp)           ic_perturb_factor

    ! met_input
    real(wp)           temp, press, RH_CLEA, hgt, cld_frac

    ! chem_input
    real(wp)           numc1, numc2, numc3, numc4, &
                       mfso41, mfpom1, mfsoa1, mfbc1, mfdst1, mfncl1, &
                       mfso42, mfsoa2, mfncl2, &
                       mfdst3, mfncl3, mfso43, mfbc3, mfpom3, mfsoa3, &
                       mfpom4, mfbc4, &
                       qso2, qh2so4, qsoag, num_factor, gas_factor, &
                       h2so4_chem_prod_rate
  end type BoxInput

contains

  subroutine read_nml(filename, input) bind(c)

    character(c_char), intent(in)  :: filename(255)
    type(BoxInput),    intent(out) :: input

    character(len=255) :: nml_file
    integer            :: nml_unit, i

    ! time_input
    integer(c_int)     mam_dt, mam_nstep, mam_output_intvl

    ! cntl_input
    integer(c_int)     mdo_gaschem, mdo_cloudchem, mdo_gasaerexch, &
                       mdo_test_siamgs19, mdo_rename, mdo_newnuc, mdo_coag, &
                       mdo_calcsize, newnuc_method_user_choice, &
                       pbl_nuc_wang2008_user_choice
    logical(c_bool)    cluster_growth_arg_write
    character(len=128) cluster_growth_arg_write_fname
    integer(c_int)     num_unit, frac_unit, gas_unit
    real(wp)           ic_perturb_factor

    ! met_input
    real(wp)           temp, press, RH_CLEA, hgt, cld_frac

    ! chem_input
    real(wp)           numc1, numc2, numc3, numc4, &
                       mfso41, mfpom1, mfsoa1, mfbc1, mfdst1, mfncl1, &
                       mfso42, mfsoa2, mfncl2, &
                       mfdst3, mfncl3, mfso43, mfbc3, mfpom3, mfsoa3, &
                       mfpom4, mfbc4, &
                       qso2, qh2so4, qsoag, num_factor, gas_factor, &
                       h2so4_chem_prod_rate

    namelist /time_input/ mam_dt, mam_nstep, mam_output_intvl
    namelist /cntl_input/ mdo_gaschem, mdo_cloudchem, &
                          mdo_gasaerexch, mdo_test_siamgs19, &
                          mdo_rename, mdo_newnuc, mdo_coag, &
                          mdo_calcsize, &
                          newnuc_method_user_choice,    &
                          pbl_nuc_wang2008_user_choice, &
                          cluster_growth_arg_write,   &
                          cluster_growth_arg_write_fname,   &
                          num_unit, frac_unit, gas_unit, &
                          ic_perturb_factor
    namelist /met_input/  temp, press, RH_CLEA, hgt, cld_frac
    namelist /chem_input/ numc1, numc2, numc3, numc4, &
                          mfso41, mfpom1, mfsoa1, &
                          mfbc1, mfdst1, mfncl1, &
                          mfso42, mfsoa2, mfncl2, &
                          mfdst3, mfncl3, mfso43, &
                          mfbc3, mfpom3, mfsoa3, &
                          mfpom4, mfbc4, &
                          qso2, qh2so4, qsoag, &
                          num_factor, gas_factor,   &
                          h2so4_chem_prod_rate

    do i=1,255
      nml_file(i:i) = filename(i)
    end do

    ! Read namelists
    nml_unit = 101
    open (UNIT = nml_unit, FILE = nml_file, STATUS = 'OLD')
    read (nml_unit, time_input)
    read (nml_unit, cntl_input)
    read (nml_unit, met_input)
    read (nml_unit, chem_input)
    close (nml_unit)

    ! Copy data into our input type.
#define set(x) input%x = x
    set(mam_dt)
    set(mam_nstep)
    set(mam_output_intvl)
    set(mdo_gaschem)
    set(mdo_cloudchem)
    set(mdo_gasaerexch)
    set(mdo_test_siamgs19)
    set(mdo_rename)
    set(mdo_newnuc)
    set(mdo_coag)
    set(mdo_calcsize)
    set(newnuc_method_user_choice)
    set(pbl_nuc_wang2008_user_choice)
    set(cluster_growth_arg_write)
    set(num_unit)
    set(frac_unit)
    set(gas_unit)
    set(ic_perturb_factor)
    set(temp)
    set(press)
    set(RH_CLEA)
    set(hgt)
    set(cld_frac)
    set(numc1)
    set(numc2)
    set(numc3)
    set(numc4)
    set(mfso41)
    set(mfpom1)
    set(mfsoa1)
    set(mfbc1)
    set(mfdst1)
    set(mfncl1)
    set(mfso42)
    set(mfsoa2)
    set(mfncl2)
    set(mfdst3)
    set(mfncl3)
    set(mfso43)
    set(mfbc3)
    set(mfpom3)
    set(mfsoa3)
    set(mfpom4)
    set(mfbc4)
    set(qso2)
    set(qh2so4)
    set(qsoag)
    set(num_factor)
    set(gas_factor)
    set(h2so4_chem_prod_rate)
#undef set

    do i=1,128
      input%cluster_growth_arg_write_fname(i) = cluster_growth_arg_write_fname(i:i)
    end do

  end subroutine

end module
