! chem_mods.F90
!    This F90 module file is a special version of the equivalent ACME (and CAM5) module.
!    It provides the functionality needed by the cambox offline code
!    that is used for development and testing of the modal aerosol module (MAM),
!    but (in most cases) not all the functionality of the equivalent ACME module.
!    Also, it may have been taken from a version of CAM5 that was older
!    than ACME-V0 (i.e., pre 2014).

      module chem_mods

      use shr_kind_mod, only: r8 => shr_kind_r8
      use modal_aero_config, only:  pcnst

      implicit none

      public

      integer, parameter :: imozart = 6
      integer, parameter :: gas_pcnst = pcnst - (imozart - 1)

      real(r8) :: adv_mass(gas_pcnst) = 0._r8

      end module chem_mods
