#include "modal_aero_config.inc"

module modal_aero_microp_control

! !USES:
  use haero_precision, only: wp
  use chem_mods,       only:  gas_pcnst

  use modal_aero_microp_modes, only: max_mode

! use gauss_hermite
!==> JS END
  use modal_aero_data, only:  ntot_aspectype, ntot_amode

  implicit none

  public

! !PUBLIC DATA MEMBERS:
  type :: misc_vars_aa_type
! using this derived type reduces the number of changes needed to add more mosaic diagnostics to history
     real(wp) :: ncluster_tend_nnuc_1grid
  end type misc_vars_aa_type

  logical, public :: mosaic = .true. !BSINGH -  Added logical for mosaic model

  real(wp), public  :: n_so4_monolayers_pcage = huge(1.0_wp)
! number of so4(+nh4) monolayers needed to "age" a carbon particle

  real(wp), public :: dr_so4_monolayers_pcage = huge(1.0_wp)
! thickness of the so4 monolayers (m)
! for so4(+nh4), use bi-sulfate mw and 1.77 g/cm3,
!    --> 1 mol so4(+nh4)  = 65 cm^3 --> 1 molecule = (4.76e-10 m)^3
! aging criterion is approximate so do not try to distinguish
!    sulfuric acid, bisulfate, ammonium sulfate

#if ( defined( MAM_STANDALONE ) )
  integer, public :: cldy_rh_sameas_clear = 0
! this is only used for some specific box model tests

  real(wp), public :: alpha_astem_soa_boxtest = 0.05_wp
  integer, public :: niter_max_soa_boxtest = 1000
#endif

  integer, public :: mdo_gaexch_cldy_subarea = 0
! controls if gas condensation is done in cloudy subarea
!    1 = yes ; 0 = no

  integer, public :: gaexch_h2so4_uptake_optaa = 2
! controls treatment of h2so4 condensation in mam_gasaerexch_1subarea
!    1 = sequential   calc. of gas-chem prod then condensation loss
!    2 = simultaneous calc. of gas-chem prod and  condensation loss

  integer, public :: newnuc_h2so4_conc_optaa = 2
! controls treatment of h2so4 concentrationin mam_newnuc_1subarea
!    1 = use average value calculated in standard cam5.2.10 and earlier
!    2 = use average value calculated in mam_gasaerexch_1subarea
!   11 = use average of initial and final values from mam_gasaerexch_1subarea
!   12 = use final value from mam_gasaerexch_1subarea

  integer, public :: rename_method_optaa = 40
! controls renaming parameterization

  integer, public :: update_qaerwat = 0
  integer, public :: update_dgncur_a = 0
  integer, public :: update_dgncur_awet = 0
! controls updating of qaerwat
! controls updating of dgncur_a
! controls updating of dgncur_awet and wetdens_host

  real (wp) :: newnuc_adjust_factor_dnaitdt = 1.0_wp
  real (wp) :: newnuc_adjust_factor_pbl     = 1.0_wp


  integer, parameter :: maxsubarea  = 2
!==> JS changes to 5, add condensation alone
  integer, parameter :: nqtendaa    = 5     ! original value: 4
!==> JS END
  integer, parameter :: iqtend_cond = 1
  integer, parameter :: iqtend_rnam = 2
  integer, parameter :: iqtend_nnuc = 3
  integer, parameter :: iqtend_coag = 4
!==> JS ADD
  integer, parameter :: iqtend_cond_only = 5
!==> JS END
  integer, parameter :: nqqcwtendaa = 1
  integer, parameter :: iqqcwtend_rnam = 1
!==> JS changes to 5 for consistency
  integer, parameter :: iqqcwtend_match_iqtend(nqtendaa) = (/ 0, iqqcwtend_rnam, 0, 0, 0 /)
!==> JS END
  logical, parameter :: aging_include_seasalt = .false.
                      ! when .true., aging (by coagulation) includes contribution of seasalt
                      ! early versions of mam neglected the seasalt contribution

  integer :: mode_aging_optaa(max_mode)

  integer :: lun82,   lun97,   lun98,   lun13n,   lun15n
  logical :: ldiag82, ldiag97, ldiag98, ldiag13n, ldiag15n
  logical :: ldiagd1



! following were used in aging calcs but are no longer needed
!    fac_m2v_so4, fac_m2v_nh4, fac_m2v_soa(:)
!    fac_m2v_pcarbon(:)
!    soa_equivso4_factor(:)


!==> JS adds tendency for condensation alone
  character(len=8) :: suffix_q_coltendaa(nqtendaa) = &
     (/ '_sfgaex1', '_sfgaex2', '_sfnnuc1', '_sfcoag1', '_sfcond1' /)
!==> JS END
  character(len=8) :: suffix_qqcw_coltendaa(nqqcwtendaa) = &
                    '_sfgaex2'

  logical :: do_q_coltendaa(gas_pcnst,nqtendaa) = .false.
  logical :: do_qqcw_coltendaa(gas_pcnst,nqqcwtendaa) = .false.

end module modal_aero_microp_control
