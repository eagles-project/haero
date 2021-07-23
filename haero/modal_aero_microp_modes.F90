#include "modal_aero_config.inc"

module modal_aero_microp_modes
!------------------------------------------------------------------------------
! Control variables, constants, and parameters that define/characterize
! the lognormal modes used by MAM's aerosol microphysics parameterization suite
! for representing the aerosol (particle) size distibution:
!
! Most of the variables in this module were declared in modal_aero_amicphys.F90
! in EAMv1. This module was created (together with modal_aero_microp_species.F90
! and modal_aero_microp_host_mapping.F90) during the refactoring of MAM's 
! microphysics code as part of the EAGLES project.
!
! POC for the refactored code: Hui Wan, PNNL, Hui.Wan@pnnl.gov
!------------------------------------------------------------------------------

  use haero_precision, only: wp
  use modal_aero_data,    only: ntot_amode

  implicit none

  public
  !---------------------------------------------------------------------------------------
  ! Number of lognormal modes used in the representation of the aerosol size distribution
  !---------------------------------------------------------------------------------------

  integer, parameter :: ntot_amode_extd = ntot_amode
  integer, parameter :: max_mode_fresh = 1
  integer, parameter :: max_mode = ntot_amode_extd + max_mode_fresh

  !---------------------------------------------------------------------------------------
  ! Fixed mode parameters 
  !---------------------------------------------------------------------------------------
  real(wp) :: sigmag_aer(max_mode)
  real(wp) :: alnsg_aer(max_mode)
  real(wp) :: dgnumhi_aer(max_mode), dgnumlo_aer(max_mode)

  ! factor for converting number geometric_mean diameter to volume-mean diameter
  real(wp) :: fcvt_dgnum_dvolmean(max_mode)

  ! The next variable, dgnum_aer, is the geometric dry mean diameter (m) of 
  ! the number distribution for aerosol mode m. According to Dick's comments in 
  ! the code, this variable is used only when aerosol number is not simulated.
  ! So this seems to be a legacy variable that can be removed.
  real(wp) :: dgnum_aer(max_mode)

  !---------------------------------------------------------------------------------------
  ! Pointers to specific modes
  !---------------------------------------------------------------------------------------
  ! Hui Wan 2021-02: the following are mode indices? should be named ixxx instead of nxxx?

  integer :: nacc, nait, npca, nufi, nmacc, nmait

  !---------------------------------------------------------------------------------------
  ! The next set of variables 
  !
  !needed by both aging and coag modules
 !integer, parameter :: max_agepair  = 1    
 !integer :: i_agepair_pca, i_agepair_macc, i_agepair_mait

contains


  subroutine set_additional_constants_for_modes
!------------------------------------------------------------------
! Purpose:
!  Calculate various constants for the lognormal modes
!
! History:
!  Original version by R. Easter
!  Moved here from modal_aero_amicphys.F90 by Hui Wan, PNNL, 2021-02
!------------------------------------------------------------------

     alnsg_aer(1:max_mode) = log(sigmag_aer(1:max_mode))

     ! factor for converting number geometric_mean diameter to volume-mean diameter
     fcvt_dgnum_dvolmean(1:max_mode) = exp( 1.5_wp*(alnsg_aer(1:max_mode)**2) )

  end subroutine set_additional_constants_for_modes

end module modal_aero_microp_modes
