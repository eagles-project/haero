!> This module defines precision used by haero

module haero_precision

  use iso_c_binding

  implicit none

  private

  public :: wp

  !> Working precision real kind
  integer, parameter :: wp = c_real

end module haero_precision
