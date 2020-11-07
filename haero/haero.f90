!> This module contains data structures that allow Fortran modules to access
!> data from the C++ model, prognostics, and diagnostics.
module haero

  use iso_c_binding

  implicit none

  !> Working precision real kind
  integer, parameter :: wp = c_double

  private
  public :: Model, Prognostics, Diagnostics

  !> This Fortran type is the equivalent of the C++ Mode struct.
  type, bind(c) :: Mode
    character(len=32) :: name
    real(wp) :: min_diameter
    real(wp) :: max_diameter
    real(wp) :: mean_std_dev
  end type

  !> This Fortran type is the equivalent of the C++ Species struct.
  type, bind(c) :: Species
    character(len=32) :: name
    character(len=8) :: symbol
  end type

  !> This Fortran type is the equivalent of the C++ Model class.
  type :: Model
  private
    type(c_ptr) :: ptr
  contains
    procedure :: from_c_ptr => m_from_c_ptr
    procedure :: modes => m_modes
    procedure :: aero_species => m_aero_species
    procedure :: aero_species_for_mode => m_aero_species_for_mode
    procedure :: gas_species => m_gas_species
    procedure :: num_columns => m_num_columns
    procedure :: num_levels => m_num_levels
  end type

  !> This type represents the set of prognostic variables for an aerosol
  !> model.
  type :: Prognostics
  private
    type(c_ptr) :: ptr
  contains
    procedure :: from_c_ptr => p_from_c_ptr
    procedure :: num_aero_modes => p_num_aero_modes
    procedure :: num_aero_species => p_num_aero_species
    procedure :: num_gas_species => p_num_gas_species
    procedure :: num_columns => p_num_columns
    procedure :: num_levels => p_num_levels
    procedure :: int_aero_mix_frac => p_int_aero_mix_frac
    procedure :: cld_aero_mix_frac => p_cld_aero_mix_frac
    procedure :: gas_mole_frac => p_gas_mole_frac
    procedure :: modal_num_densities => p_modal_num_densities
  end type

  !> This type represents the set of diagnostic variables for an aerosol
  !> model.
  type :: Diagnostics
  private
    type(c_ptr) :: ptr
  contains
    procedure :: from_c_ptr => d_from_c_ptr
    procedure :: num_aero_modes => d_num_aero_modes
    procedure :: num_aero_species => d_num_aero_species
    procedure :: num_gas_species => d_num_gas_species
    procedure :: num_columns => d_num_columns
    procedure :: num_levels => d_num_levels
    procedure :: has_var => d_has_var
    procedure :: var => d_var
    procedure :: has_modal_var => d_has_modal_var
    procedure :: modal_var => d_modal_var
  end type

contains

  !> Returns the number of aerosol modes in the given prognostics object.
  !> @param [in] p A Prognostics object.
  integer(c_int) function p_num_aero_modes(p)
    type(Prognostics), intent(in) :: p
    p_num_aero_modes = p_num_aero_modes_c(p%ptr)
  end function

  !> Returns the number of aerosol species in the given prognostics object.
  !> @param [in] p A Prognostics object.
  integer(c_int) function p_num_aero_species(p)
    type(Prognostics), intent(in) :: p
    p_num_aero_species = p_num_aero_species_c(p%ptr)
  end function

  !> Returns the number of gas species in the given prognostics object.
  !> @param [in] p A Prognostics object.
  integer(c_int) function p_num_gas_species(p)
    type(Prognostics), intent(in) :: p
    p_num_gas_species = p_num_gas_species_c(p%ptr)
  end function

  !> Returns the number of columns in the given prognostics object.
  !> @param [in] p A Prognostics object.
  integer(c_int) function p_num_columns(p)
    type(Prognostics), intent(in) :: p
    p_num_columns = p_num_columns_c(p%ptr)
  end function

  !> Returns the number of vertical levels in the given prognostics object.
  !> @param [in] p A Prognostics object.
  integer(c_int) function p_num_levels(p)
    type(Prognostics), intent(in) :: p
    p_num_levels = p_num_levels_c(p%ptr)
  end function

  !> Provides access to the interstitial aerosol mixing fractions array
  !> for the given mode in the given prognostics object.
  !> @param [in] p A pointer to a prognostics object.
  !> @param [in] mode An index identifying the desired mode.
  function p_int_aero_mix_frac(p, mode) result(data)
    type(Prognostics), intent(in)  :: p
    integer(c_int), intent(in) :: mode
    real(c_real), pointer, dimension(:,:,:), intent(out) :: data

    type(c_ptr) :: v_ptr
    integer(c_int), dimension(3) :: v_shape
    call p_get_int_aero_mix_frac_c(p%ptr, mode, v_ptr, data_ѕhape)
    call c_f_pointer(v_ptr, data, v_shape)
  end function

  !> Provides access to the cloud-borne aerosol mixing fractions array
  !> for the given mode in the given prognostics object.
  !> @param [in] p A pointer to a prognostics object.
  !> @param [in] mode An index identifying the desired mode.
  function p_cld_aero_mix_frac(p, mode) result(data)
    type(Prognostics), intent(in)  :: p
    integer(c_int), intent(in) :: mode
    real(c_real), pointer, dimension(:,:,:), intent(out) :: data

    type(c_ptr) :: v_ptr
    integer(c_int), dimension(3) :: v_shape
    call p_get_cld_aero_mix_frac_c(p%ptr, mode, v_ptr, data_ѕhape)
    call c_f_pointer(v_ptr, data, v_shape)
  end function

  !> Provides access to the gas mole fractions array for the given
  !> prognostics object.
  !> @param [in] p A Prognostics object.
  function p_gas_mole_frac(p) result(data)
    use iso_c_binding, only: c_ptr, c_int
    type(c_ptr), value, intent(in) :: p
    real(c_real), pointer, dimension(:,:,:), intent(out) :: data

    type(c_ptr) :: v_ptr
    integer(c_int), dimension(3) :: v_shape
    call p_get_gas_mole_frac_c(p%ptr, v_ptr, v_shape)
    call c_f_pointer(v_ptr, data, v_shape)
  end function

  !> Provides access to the modal number fractions array for the given
  !> prognostics object.
  !> @param [in] p A Prognostics object.
  function p_modal_num_densities(p) result(data)
    use iso_c_binding, only: c_ptr, c_int
    type(c_ptr), value, intent(in) :: p
    real(c_real), pointer, dimension(:,:,:), intent(out) :: data

    type(c_ptr) :: v_ptr
    integer(c_int), dimension(3) :: v_shape
    call p_get_modal_num_densities_c(p%ptr, v_ptr, data_ѕhape)
    call c_f_pointer(v_ptr, data, v_shape)
  end function

  interface

    integer(c_int) function p_num_aero_modes_c(p) bind(c)
      use iso_c_binding, only: c_int
      type(Prognostics), intent(in) :: p
    end function

    integer(c_int) function p_num_aero_species_c(p) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(Prognostics), value, intent(in) :: p
    end function

    integer(c_int) function p_num_gas_species_c(p) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
    end function

    integer(c_int) function p_num_columns_c(p) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
      integer(c_int), value, intent(out) :: num_columns
    end function

    integer(c_int) function p_num_levels_c(p) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
      type(c_int), value, intent(out) :: p_num_levels_c
    end function

    subroutine p_get_int_aero_mix_frac_c(p, mode, v_ptr, v_shape) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
      integer(int), value, intent(in) :: mode
      real(c_real), dimension(:,:,:), intent(out) :: v_ptr
      integer(c_int), dimension(3), intent(out) :: v_shape
    end function

    subroutine p_get_cld_aero_mix_frac_c(p, mode, v_ptr, v_shape) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
      integer(c_int), value, intent(in) :: mode
      real(c_real), dimension(:,:,:), intent(out) :: v_ptr
      integer(c_int), dimension(3), intent(out) :: v_shape
    end function

    subroutine p_get_gas_mole_frac_c(p, v_ptr, v_shape) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
      real(c_real), dimension(:,:,:), intent(out) :: v_ptr
      integer(c_int), dimension(3), intent(out) :: v_shape
    end function

    subroutine p_get_modal_num_densities_c(p, v_ptr, v_shape) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
      real(c_real), dimension(:,:,:), intent(out) :: v_ptr
      integer(c_int), dimension(3), intent(out) :: v_shape
    end function
  end interface

end module

