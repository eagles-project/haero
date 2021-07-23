#include "modal_aero_config.inc"

module modal_aero_microp_species
!------------------------------------------------------------------------------
! This module contains
!  - control variables and indices used for keeping track of 
!    gas- or condensed-phase chemical species inside 
!    MAM's aerosol microphysics parameterization suite;
!  - constants describing physical/chemical properties of those species;
!  - subroutines for assigning values to those variables at the beginning
!    of a simulation.
!
! Most of the variables and lines of calculation in this module were 
! placed in modal_aero_amicphys.F90 in EAMv1. This module was created 
! (together with modal_aero_microp_modes.F90 and 
! modal_aero_microp_host_mapping.F90) during the refactoring of MAM's 
! microphysics code as part of the EAGLES project.
!
! POC for the refactored code: Hui Wan, PNNL, Hui.Wan@pnnl.gov
!------------------------------------------------------------------------------

  use haero_precision, only: wp 
  use modal_aero_data,    only: nsoa, npoa, nbc

  implicit none

  public
  !---------------------
  ! gas species
  !---------------------
  ! dimension sizes

  integer, parameter :: max_gas = nsoa + 1   !used for determining dimension size of arrays
  integer :: ngas                            !actually registered # of gases

  ! species indices 

  integer :: igas_h2so4, igas_nh3, igas_hno3, igas_hcl

  integer :: igas_soag, igas_soagzz  ! first and last soag species

  ! physical/chemical properties

  real(wp) :: mw_gas(max_gas)         ! molecular weight
  real(wp) :: vol_molar_gas(max_gas)  ! molar volume
  real(wp) :: accom_coef_gas(max_gas) ! accommodation coefficient

  real(wp) :: fcvt_gas(max_gas)

  !---------------------
  ! aerosol species
  !---------------------
  ! dimension sizes

  integer, parameter :: max_aer = nsoa + npoa + nbc + 4 ! the +4 in max_aer are dst, ncl, so4, mom
  integer :: naer

  ! species indices
  !    when nsoa > 1, iaer_soa is index of the first soa species, iaer_soazz is the last soa species
  !    when nbc  > 1, iaer_bc  is index of the first bc  species
  !    when npom > 1, iaer_pom is index of the first pom species

  integer :: iaer_bc, iaer_dst, iaer_ncl, iaer_nh4, iaer_pom, iaer_soa, iaer_soazz, iaer_so4, &
             iaer_mpoly, iaer_mprot, iaer_mlip, iaer_mhum, iaer_mproc, iaer_mom, &
             iaer_no3, iaer_cl, iaer_ca, iaer_co3

  ! physical/chemical properties of species
  ! molecular weights of SOA and POA can be different from values in the host model

  real(wp) :: mwuse_soa(nsoa)
  real(wp) :: mwuse_poa(npoa)  

  real(wp) :: mw_aer(max_aer)     ! molecular weights of all aerosol species
  real(wp) :: dens_aer(max_aer)   ! density

  real(wp) :: hygro_aer(max_aer)  ! hygroscopicity

  real(wp) :: fac_m2v_aer(max_aer)        ! converts (mol-aero/mol-air) to (m3-aero/mol-air)
  real(wp) :: fac_eqvso4hyg_aer(max_aer)  ! converts a species volume to a volume of so4
!                                         !    (or nh4hso4) having same hygroscopicity
  real(wp) :: fac_m2v_eqvhyg_aer(max_aer) ! = fac_m2v_aer * fac_eqvso4hyg_aer

  real(wp) :: fcvt_aer(max_aer), fcvt_num, fcvt_wtr

  !------------------------------------
  ! Mapping between gases and aerosols
  !------------------------------------
  ! If a gas condenses, which aerosol species will it become?

  integer :: idx_gas_to_aer(max_gas)

contains

   subroutine set_gas_and_aer_names_and_indices( ntot_amode, name_gas, name_aerpfx )
!----------------------------------------------------------------------------------------------------------------
! Purpose:
! Register gas and aerosol species for the microphysics suite's internal bookkeeping.
! - Module variables that will be assigned values in this subroutine
!      ngas, naer
!      igas_soag, igas_soagzz, igas_h2so4, igas_nh3, igas_hno3, igas_hcl
!      iaer_bc, iaer_dst, iaer_ncl, iaer_nh4, iaer_pom, iaer_soa, iaer_so4, 
!      iaer_mpoly, iaer_mprot, iaer_mlip, iaer_mhum, iaer_mproc, iaer_mom, 
!      iaer_no3, iaer_cl, iaer_ca, iaer_co3
!
! - Module variablels with already assinged values that will be used in this subroutine
!      nsoa, npoa, nbc
!
! History
!   Original version by R. C. Easter
!   Packed into subroutine, slightly revised for generality, and moved to this module by Hui Wan, PNNL, 2021-02
!----------------------------------------------------------------------------------------------------------------
   use modal_aero_logging, only: iulog => iulog_main
   use abortutils,  only:  endrun

   integer, intent(in) :: ntot_amode
   character(len=16),intent(out) :: name_gas(max_gas), name_aerpfx(max_aer)

   ! local variables
   integer  :: jsoa, jsoag

   !-----------------------------

      name_gas    = "???"
      name_aerpfx = "???"

      idx_gas_to_aer(:) = -1

      igas_h2so4 = 0 ; igas_nh3 = 0
      igas_soag  = 0 ; igas_soagzz = -1

      iaer_bc  = 0 ; iaer_dst = 0
      iaer_ncl = 0 ; iaer_nh4 = 0
      iaer_pom = 0 

      iaer_soa = 0 ; iaer_soazz = -1

      iaer_so4 = 0
      iaer_no3 = 0 ; iaer_cl  = 0
      iaer_ca  = 0 ; iaer_co3 = 0
      iaer_mpoly = 0 ; iaer_mprot = 0
      iaer_mlip  = 0 ; iaer_mhum = 0
      iaer_mproc = 0 ; iaer_mom = 0

      ngas = 0
      naer = 0

      !----------------
      ! SOA and SOAG
      !----------------

      igas_soag = ngas + 1
      iaer_soa  = naer + 1

      if (nsoa == 1) then
         name_gas   (igas_soag) = 'SOAG'
         name_aerpfx(iaer_soa ) = 'soa'

      else if (nsoa == 2) then
         jsoag = igas_soag ; name_gas(jsoag) = 'SOAGa'  ; jsoa = iaer_soa ; name_aerpfx(jsoa) = 'soaa'  ! jsoa=1
         jsoag = jsoag+1   ; name_gas(jsoag) = 'SOAGb'  ; jsoa = jsoa+1   ; name_aerpfx(jsoa) = 'soab'  ! jsoa=2

      else if (nsoa == 6) then
         jsoag = igas_soag ; name_gas(jsoag) = 'SOAGa1' ; jsoa = iaer_soa ; name_aerpfx(jsoa) = 'soaa1' ! jsoa=1
         jsoag = jsoag+1   ; name_gas(jsoag) = 'SOAGa2' ; jsoa = jsoa+1   ; name_aerpfx(jsoa) = 'soaa2' ! jsoa=2
         jsoag = jsoag+1   ; name_gas(jsoag) = 'SOAGa3' ; jsoa = jsoa+1   ; name_aerpfx(jsoa) = 'soaa3' ! jsoa=3
         jsoag = jsoag+1   ; name_gas(jsoag) = 'SOAGb1' ; jsoa = jsoa+1   ; name_aerpfx(jsoa) = 'soab1' ! jsoa=4
         jsoag = jsoag+1   ; name_gas(jsoag) = 'SOAGb2' ; jsoa = jsoa+1   ; name_aerpfx(jsoa) = 'soab2' ! jsoa=5
         jsoag = jsoag+1   ; name_gas(jsoag) = 'SOAGb3' ; jsoa = jsoa+1   ; name_aerpfx(jsoa) = 'soab3' ! jsoa=6
      else
         call endrun( 'modal_aero_amicphys_init ERROR - bad nsoa' )
      end if
      ngas = ngas + nsoa ; igas_soagzz = ngas
      naer = naer + nsoa ; iaer_soazz  = naer

      ! if (nsoag == nsoa) then
           do jsoag = igas_soag,igas_soagzz
              idx_gas_to_aer(jsoag) = iaer_soa + (jsoag - igas_soag)
           end do
      ! else if (nsoa == 1) then
      !    idx_gas_to_aer(igas_soag:igas_soagzz) = iaer_soa 
      ! else
      !    call endrun( 'modal_aero_amicphys_init ERROR - check nsoa vs nsoag' )
      ! end if

      !------------------------------------
      ! Sulfuric acid gas and sulfate
      !------------------------------------
      igas_h2so4 = ngas + 1
      name_gas(igas_h2so4) = 'H2SO4'
      ngas = ngas + 1

      iaer_so4 = naer + 1
      name_aerpfx(iaer_so4) = 'so4'
      naer = naer + 1

      idx_gas_to_aer(igas_h2so4) = iaer_so4

      !------------------------------------
      ! NH3 and nitrate
      !------------------------------------
      if ( (ntot_amode==7) .or. &
           (ntot_amode==8) .or. &
           (ntot_amode==9) ) then

         igas_nh3 = ngas + 1
         name_gas(igas_nh3) = 'NH3'
         ngas = ngas + 1

         iaer_nh4 = naer + 1
         name_aerpfx(iaer_nh4) = 'nh4'
         naer = naer + 1

         idx_gas_to_aer(igas_nh3) = iaer_nh4
      else

         igas_nh3 = -999888777
         iaer_nh4 = -999888777
      end if

      !------------------------------------
      ! dummy indices
      !------------------------------------
      igas_hno3 = -999888777
      igas_hcl  = -999888777
      iaer_no3  = -999888777
      iaer_cl   = -999888777

      !------------------------------------
      ! POM
      !------------------------------------

      iaer_pom = naer + 1  ! (start-)index of POM
      if (npoa == 1) then
         name_aerpfx(iaer_pom) = 'pom'
         naer = naer + 1

      else if (npoa == 2) then
         name_aerpfx(iaer_pom)   = 'poma'
         name_aerpfx(iaer_pom+1) = 'pomb'

         naer = naer + 2  ! this is the total # of aerosol species registered so far
      else
         call endrun( 'modal_aero_amicphys_init ERROR - bad npoa' )
      end if

      !------------------------------------
      ! black carbon
      !------------------------------------
      iaer_bc = naer + 1
      if (nbc == 1) then
         name_aerpfx(iaer_bc) = 'bc'
         naer = naer + 1

      else if (nbc == 2) then
         name_aerpfx(iaer_bc)   = 'bca'
         name_aerpfx(iaer_bc+1) = 'bcb'
         naer = naer + 2  ! this is the total # of aerosol species registered so far
      else
         call endrun( 'modal_aero_amicphys_init ERROR - bad nbc' )
      end if

      !------------------------------------
      ! sea salt
      !------------------------------------
      iaer_ncl = naer + 1
      name_aerpfx(iaer_ncl) = 'ncl'
      naer = naer + 1

      !------------------------------------
      ! dust
      !------------------------------------
      iaer_dst = naer + 1
      name_aerpfx(iaer_dst) = 'dst'
      naer = naer + 1

      !------------------------------------
      ! dummy indices
      !------------------------------------
      iaer_ca   = -999888777
      iaer_co3  = -999888777

      !------------------------------------
      ! marine organic matter 
      !------------------------------------
      iaer_mom = naer + 1
      name_aerpfx(iaer_mom) = 'mom'
      naer = naer + 1

      !------------------------------------
      ! additional species 
      !------------------------------------
      if (ntot_amode==9) then
         iaer_mpoly = naer + 1 ; name_aerpfx(iaer_mpoly) = 'mpoly' ; naer = naer + 1
         iaer_mprot = naer + 1 ; name_aerpfx(iaer_mprot) = 'mprot' ; naer = naer + 1
         iaer_mlip  = naer + 1 ; name_aerpfx(iaer_mlip ) = 'mlip'  ; naer = naer + 1
         iaer_mhum  = naer + 1 ; name_aerpfx(iaer_mhum ) = 'mhum'  ; naer = naer + 1
         iaer_mproc = naer + 1 ; name_aerpfx(iaer_mproc) = 'mproc' ; naer = naer + 1
      else
         iaer_mpoly = -999888777
         iaer_mprot = -999888777
         iaer_mlip  = -999888777
         iaer_mhum  = -999888777
         iaer_mproc = -999888777
      end if

      !------------------------------------
      ! Sanity check
      !------------------------------------
      if ((ngas /= max_gas) .or. (naer /= max_aer)) then
         write(iulog,'(a,4i10)') 'ngas, max_gas, naer, max_aer', &
            ngas, max_gas, naer, max_aer
         call endrun( 'modal_aero_amicphys_init ERROR - bad ngas or naer' )
      end if

end subroutine set_gas_and_aer_names_and_indices

subroutine set_constants_for_gases( mwhost_gas )
!--------------------------------------------------------------------------------------------
! Purpose:
!  Set the following constants for gas species
!  - molecular weight
!  - molar volume
!  - fcvt_gas
! 
! Authors: 
!  Original version by R. Easter.
!  Slightly refactored and moved to here from modal_aero_amicphys_init by Hui Wan, Feb, 2021
!--------------------------------------------------------------------------------------------
! Arguments

  real(wp) :: mwhost_gas(max_gas)

! Local variables

  integer :: igas, i_th_soa
  logical :: is_soa

!---------------------------------

      mwuse_soa(:) = 150.0_wp
      mwuse_poa(:) = 150.0_wp

      ! Initialize arrays

      vol_molar_gas = 42.88_wp  ! value for h2so4
      mw_gas(:) = 1.0_wp
      fcvt_gas(:) = 1.0_wp

      ! Copy/map values from host model

      do igas = 1, ngas

         i_th_soa = igas - igas_soag + 1 
         is_soa = (igas>=igas_soag).and.(igas<=igas_soagzz) 

         ! molecular weight
         mw_gas(igas) = mwhost_gas(igas)                   ! copy over value from host model
         if (is_soa) mw_gas(igas) = mwuse_soa(i_th_soa)    ! reset value for soa
         fcvt_gas(igas) = mwhost_gas(igas)/mw_gas(igas)    ! ratio between host value and MAM value 

         ! molar volume
         if (is_soa) then
            vol_molar_gas(igas) = vol_molar_gas(igas_h2so4) * (mw_gas(igas)/98.0_wp)
         else if (igas == igas_nh3) then
            vol_molar_gas(igas) = 14.90_wp
         else if (igas == igas_hno3) then
            vol_molar_gas(igas) = 24.11_wp
         else if (igas == igas_hcl) then
            vol_molar_gas(igas) = 21.48_wp
         end if
      end do ! igas

      ! Accommodation coefficient (for condensation)

      accom_coef_gas = 0.65_wp  ! value for h2so4

  end subroutine set_constants_for_gases



  subroutine set_constants_for_aerosol_species( mwhost_aer, mwhost_dry_air, mwhost_h2o, &
                                                dens_aer_host, hygro_aer_host )

!--------------------------------------------------------------------------------------------
! Purpose:
!  Set various constants for aerosol species
! 
! Authors: 
!  Original version by R. Easter.
!  Slightly refactored and moved to here from modal_aero_amicphys_init by Hui Wan, Feb, 2021
!--------------------------------------------------------------------------------------------

      real(wp), intent(in) :: mwhost_aer(max_aer)
      real(wp), intent(in) :: mwhost_dry_air
      real(wp), intent(in) :: mwhost_h2o

      real(wp), intent(in) ::  dens_aer_host(max_aer)
      real(wp), intent(in) :: hygro_aer_host(max_aer)

      integer :: iaer

      logical :: is_soa, is_poa
      integer :: i_th_soa, i_th_poa


      !------------------------------------
      ! density

       dens_aer(1:naer) =  dens_aer_host(1:naer)  

      !--------------------------------------------
      ! molecular weights and related constants

      mw_aer(:) = 1.0_wp
      fcvt_aer(:) = 1.0_wp
      do iaer = 1,naer

         mw_aer(iaer) = mwhost_aer(iaer)  

         i_th_soa = iaer - iaer_soa + 1 ; is_soa = (i_th_soa>=1).and.(i_th_soa<=nsoa) 
         i_th_poa = iaer - iaer_pom + 1 ; is_poa = (i_th_poa>=1).and.(i_th_poa<=npoa) 
         
         if (is_soa) then
            mw_aer(iaer) = mwuse_soa(i_th_soa)
         else if (is_poa) then
            mw_aer(iaer) = mwuse_poa(i_th_poa)
         end if

         fcvt_aer(iaer) = mwhost_aer(iaer)/mw_aer(iaer)

      end do ! iaer

      ! aerosol number
      fcvt_num = 1.0_wp  ! leave number mix-ratios unchanged (#/kmol-air)

      ! factor for converting aerosol water mix-ratios from (kg/kg) to (mol/mol)
      fcvt_wtr = mwhost_dry_air/mwhost_h2o  

      ! factor for converting (kmol-AP/kmol-air) to (m3-AP/kmol-air)
      fac_m2v_aer(1:naer) = mw_aer(1:naer)/dens_aer(1:naer)

      !------------------------------------
      ! hygroscopicity

      hygro_aer(1:naer) = hygro_aer_host(1:naer)  

      ! equivalent hygroscopicity wrt sulfate

       fac_eqvso4hyg_aer(1:naer) =   hygro_aer(1:naer) / hygro_aer(iaer_so4)
      fac_m2v_eqvhyg_aer(1:naer) = fac_m2v_aer(1:naer) * fac_eqvso4hyg_aer(1:naer)

  end subroutine set_constants_for_aerosol_species


end module modal_aero_microp_species
