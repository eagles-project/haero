!> This module contains functions for converting between different aerosol
!> quantities in Fortran processes.
module haero_conversions

  use iso_c_binding

  implicit none

  private

  public :: relative_humidity_from_vapor_mixing_ratio

  interface

    real(c_real) function relative_humidity_from_vapor_mixing_ratio(qv, p, T) &
      bind(c, name="relative_humidity_from_vapor_mixing_ratio_f")
      use iso_c_binding, only: c_real
      real(c_real), value, intent(in) :: qv, p, T
    end function

  end interface

end module
