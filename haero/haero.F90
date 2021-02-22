!> This module contains data structures that allow Fortran modules to access
!> data from the C++ model, prognostics, and diagnostics.
module haero

  use iso_c_binding

  implicit none

  private

  public :: wp, mode_t, species_t, model_t, &
            prognostics_t, atmosphere_t, diagnostics_t, tendencies_t, &
            prognostics_from_c_ptr, atmosphere_from_c_ptr, &
            diagnostics_from_c_ptr, tendencies_from_c_ptr, model

  !> Working precision real kind
  integer, parameter :: wp = c_real

  !> This Fortran type is the equivalent of the C++ Mode struct.
  type :: mode_t
    !> Mode name
    character(len=:), allocatable :: name
    !> Minimum particle diameter
    real(wp) :: min_diameter
    !> Maximum particle diameter
    real(wp) :: max_diameter
    !> Geometric mean standard deviation
    real(wp) :: mean_std_dev
  end type

  !> This Fortran type is the equivalent of the C++ Species struct.
  type :: species_t
    !> Species name
    character(len=:), allocatable :: name
    !> Species symbol (abbreviation)
    character(len=:), allocatable :: symbol
    !> Molecular weight [g/mol]
    real(wp) :: molecular_wt
    !> Crystalization point [?]
    real(wp) :: crystal_pt
    !> Deliquenscence point [?]
    real(wp) :: deliques_pt
  end type

  !> This Fortran type is the equivalent of the C++ Model class. Exactly one
  !> read-only instance of a model is available to Fortran processes.
  type :: model_t
    !> The aerosol modes in the model, in indexed order.
    type(mode_t), dimension(:), allocatable :: modes
    !> The number of modes in the model. Equal to size(modes).
    integer :: num_modes
    !> The number of actual species that exist within each mode.
    integer, dimension(:), allocatable :: num_mode_species
    !> population index offsets for modes.
    integer, dimension(:), allocatable :: population_offsets
    !> The total number of distinct aerosol populations.
    integer :: num_populations
    !> The aerosol species within each mode. Indexed as (mode, species).
    type(species_t), dimension(:,:), allocatable :: aero_species
    !> The gas species in the model.
    type(species_t), dimension(:), allocatable :: gas_species
    !> The number of vertical levels in an atmospheric column.
    integer :: num_levels
  contains
    !> Given the index of an aerosol population, retrieve its mode and
    !> (modal) species indices.
    procedure :: get_mode_and_species => m_get_mode_and_species
    !> Given the name of a mode, retrieve its index.
    procedure :: mode_index => m_mode_index
    !> Given a mode index and the symbolic name of an aerosol species, retrieve
    !> its index within that mode
    procedure :: aerosol_index => m_aerosol_index
    !> Given mode and aerosol species indices, retrieve a population index
    !> that can be used to access aerosol data.
    procedure :: population_index => m_population_index
    !> Given the symbolic name of a gas, retrieve its index.
    procedure :: gas_index => m_gas_index
  end type

  !> The resident model instance, available to the single allowable C++ model
  !> instance that uses Fortran-backed processes.
  type(model_t) :: model

  !> This type represents the set of prognostic variables for an aerosol
  !> model.
  type :: prognostics_t
    type(c_ptr) :: ptr
  contains
    procedure :: interstitial_aerosols => p_int_aero_mix_frac
    procedure :: cloudborne_aerosols => p_cld_aero_mix_frac
    procedure :: gases => p_gases
    procedure :: modal_num_concs => p_modal_num_concs
  end type

  !> This type represents the set of atmospheric state variables for an
  !> aerosol model.
  type :: atmosphere_t
    type(c_ptr) :: ptr
  contains
    procedure :: temperature => a_temperature
    procedure :: pressure => a_pressure
    procedure :: relative_humidity => a_relative_humidity
    procedure :: height => a_height
  end type

  !> This type represents the set of diagnostic variables for an aerosol
  !> model.
  type :: diagnostics_t
    type(c_ptr) :: ptr
  contains
    procedure :: has_var => d_has_var
    procedure :: var => d_var
    procedure :: has_aerosol_var => d_has_aerosol_var
    procedure :: aerosol_var => d_aerosol_var
    procedure :: has_gas_var => d_has_gas_var
    procedure :: gas_var => d_gas_var
    procedure :: has_modal_var => d_has_modal_var
    procedure :: modal_var => d_modal_var
  end type

  !> This type represents a set of tendencies to be computed by a prognostic
  !> process.
  type :: tendencies_t
    type(c_ptr) :: ptr
  contains
    procedure :: interstitial_aerosols => t_int_aero_mix_frac
    procedure :: cloudborne_aerosols => t_cld_aero_mix_frac
    procedure :: gases => t_gases
    procedure :: modal_num_concs => t_modal_num_concs
  end type

  interface

    type(c_ptr) function p_int_aero_mix_frac_c(p) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
    end function

    type(c_ptr) function p_cld_aero_mix_frac_c(p) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
    end function

    type(c_ptr) function p_gases_c(p) bind(c)
      use iso_c_binding, only: c_ptr
      type(c_ptr), value, intent(in) :: p
    end function

    type(c_ptr) function p_modal_num_concs_c(p) bind(c)
      use iso_c_binding, only: c_ptr
      type(c_ptr), value, intent(in) :: p
    end function

    type(c_ptr) function a_temperature_c(a) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: a
    end function

    type(c_ptr) function a_pressure_c(a) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: a
    end function

    type(c_ptr) function a_relative_humidity_c(a) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: a
    end function

    type(c_ptr) function a_height_c(a) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: a
    end function

    integer(c_int) function d_has_var_c(d, name) bind(c)
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: d
      type(c_ptr), value, intent(in) :: name
    end function

    type(c_ptr) function d_var_c(p, token) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
      integer(c_int), value, intent(in) :: token
    end function

    integer(c_int) function d_has_aerosol_var_c(d, name) bind(c)
      use iso_c_binding, only: c_bool, c_ptr, c_int
      type(c_ptr), value, intent(in) :: d
      type(c_ptr), value, intent(in) :: name
    end function

    type(c_ptr) function d_aerosol_var_c(d, token) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: d
      integer(c_int), value, intent(in) :: token
    end function

    integer(c_int) function d_has_gas_var_c(d, name) bind(c)
      use iso_c_binding, only: c_bool, c_ptr, c_int
      type(c_ptr), value, intent(in) :: d
      type(c_ptr), value, intent(in) :: name
    end function

    type(c_ptr) function d_gas_var_c(d, token) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: d
      integer(c_int), value, intent(in) :: token
    end function

    integer(c_int) function d_has_modal_var_c(d, name) bind(c)
      use iso_c_binding, only: c_bool, c_ptr, c_int
      type(c_ptr), value, intent(in) :: d
      type(c_ptr), value, intent(in) :: name
    end function

    type(c_ptr) function d_modal_var_c(d, token) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: d
      integer(c_int), value, intent(in) :: token
    end function

    type(c_ptr) function t_int_aero_mix_frac_c(t) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: t
    end function

    type(c_ptr) function t_cld_aero_mix_frac_c(t) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: t
    end function

    type(c_ptr) function t_gases_c(t) bind(c)
      use iso_c_binding, only: c_ptr
      type(c_ptr), value, intent(in) :: t
    end function

    type(c_ptr) function t_modal_num_concs_c(t) bind(c)
      use iso_c_binding, only: c_ptr
      type(c_ptr), value, intent(in) :: t
    end function

  end interface

contains

  !> This helper function converts the given C string to a Fortran string.
  function c_to_f_string(c_string) result(f_string)
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), value, intent(in) :: c_string
    character(len=:), pointer      :: f_ptr
    character(len=:), allocatable  :: f_string
    integer(kind=c_size_t)         :: c_string_len

    interface
        function c_strlen(str_ptr) bind (c, name = "strlen" ) result(len)
        use, intrinsic :: iso_c_binding
            type(c_ptr), value     :: str_ptr
            integer(kind=c_size_t) :: len
        end function c_strlen
    end interface

    call c_f_pointer(c_string, f_ptr )
    c_string_len = c_strlen(c_string)

    f_string = f_ptr(1:c_string_len)
  end function c_to_f_string

  !> This helper function converts the given Fortran string to a C string.
  function f_to_c_string(f_string) result(c_string)
    use, intrinsic :: iso_c_binding
    implicit none
    character(len=*), target :: f_string
    character(len=:), pointer :: f_ptr
    type(c_ptr) :: c_string

    interface
        function new_c_string(f_str_ptr, f_str_len) bind (c) result(c_string)
        use, intrinsic :: iso_c_binding
            type(c_ptr), value :: f_str_ptr
            integer(c_int), value :: f_str_len
            type(c_ptr) :: c_string
        end function new_c_string
    end interface

    f_ptr => f_string
    c_string = new_c_string(c_loc(f_ptr), len(f_string))
  end function f_to_c_string

  ! Begin the process of initializing the Haero Fortran module.
  subroutine haerotran_begin_init() bind(c)
    ! Nothing here yet!
  end subroutine

  ! Set the number of modes in the global model.
  subroutine haerotran_set_num_modes(num_modes) bind(c)
    use iso_c_binding, only: c_int
    implicit none

    integer(c_int), value, intent(in) :: num_modes

    model%num_modes = num_modes
    allocate(model%modes(num_modes))
    allocate(model%num_mode_species(num_modes))
    allocate(model%population_offsets(num_modes+1))
    model%num_mode_species(:) = 0
    model%num_populations = 0
  end subroutine

  subroutine haerotran_set_max_mode_species(max_num_species) bind(c)
    use iso_c_binding, only: c_int
    implicit none

    integer(c_int), value, intent(in) :: max_num_species

    allocate(model%aero_species(size(model%modes), max_num_species))
  end subroutine

  subroutine haerotran_set_mode(mode, name, min_d, max_d, std_dev) bind(c)
    use iso_c_binding, only: c_int, c_ptr, c_real
    implicit none

    integer(c_int), value, intent(in) :: mode
    type(c_ptr), value, intent(in) :: name
    real(c_real), value, intent(in) :: min_d
    real(c_real), value, intent(in) :: max_d
    real(c_real), value, intent(in) :: std_dev

    model%modes(mode)%name = c_to_f_string(name)
    model%modes(mode)%min_diameter = min_d
    model%modes(mode)%max_diameter = max_d
    model%modes(mode)%mean_std_dev = std_dev

  end subroutine

  subroutine haerotran_set_aero_species(mode, species, name, symbol, &
    molecular_wt, crystal_pt, deliques_pt) bind(c)
    use iso_c_binding, only: c_int, c_ptr
    implicit none

    integer(c_int), value, intent(in) :: mode
    integer(c_int), value, intent(in) :: species
    type(c_ptr), value, intent(in) :: name
    type(c_ptr), value, intent(in) :: symbol
    real(c_real), value, intent(in) :: molecular_wt
    real(c_real), value, intent(in) :: crystal_pt
    real(c_real), value, intent(in) :: deliques_pt

    model%aero_species(mode, species)%name = c_to_f_string(name)
    model%aero_species(mode, species)%symbol = c_to_f_string(symbol)
    model%aero_species(mode, species)%molecular_wt = molecular_wt
    model%aero_species(mode, species)%crystal_pt = crystal_pt
    model%aero_species(mode, species)%deliques_pt = deliques_pt
    model%num_mode_species(mode) = max(species, model%num_mode_species(mode))
  end subroutine

  subroutine haerotran_set_num_gas_species(num_species) bind(c)
    use iso_c_binding, only: c_int
    implicit none

    integer(c_int), value, intent(in) :: num_species

    allocate(model%gas_species(num_species))
  end subroutine

  subroutine haerotran_set_gas_species(species, name, symbol, &
    molecular_wt, crystal_pt, deliques_pt) bind(c)
    use iso_c_binding, only: c_int, c_ptr
    implicit none

    integer(c_int), value, intent(in) :: species
    type(c_ptr), value, intent(in) :: name
    type(c_ptr), value, intent(in) :: symbol
    real(c_real), value, intent(in) :: molecular_wt
    real(c_real), value, intent(in) :: crystal_pt
    real(c_real), value, intent(in) :: deliques_pt

    model%gas_species(species)%name = c_to_f_string(name)
    model%gas_species(species)%symbol = c_to_f_string(symbol)
    model%gas_species(species)%molecular_wt = molecular_wt
    model%gas_species(species)%crystal_pt = crystal_pt
    model%gas_species(species)%deliques_pt = deliques_pt
  end subroutine

  subroutine haerotran_set_num_levels(num_levels) bind(c)
    use iso_c_binding, only: c_int
    implicit none

    integer(c_int), value, intent(in) :: num_levels

    model%num_levels = num_levels
  end subroutine

  ! Wrap up the process of initializing the Haero Fortran module.
  subroutine haerotran_end_init() bind(c)
    implicit none

    integer :: m

    model%population_offsets(1) = 0
    do m=1,model%num_modes
      model%population_offsets(m+1) = model%population_offsets(m) + model%num_mode_species(m)
    end do
    model%num_populations = model%population_offsets(model%num_modes+1)
  end subroutine

  ! This subroutine gets called when the C++ process exits.
  subroutine haerotran_finalize() bind(c)
    if (allocated(model%aero_species)) then
      deallocate(model%aero_species)
    end if
    if (allocated(model%gas_species)) then
      deallocate(model%gas_species)
    end if
    if (allocated(model%population_offsets)) then
      deallocate(model%population_offsets)
    end if
    if (allocated(model%num_mode_species)) then
      deallocate(model%num_mode_species)
    end if
    if (allocated(model%modes)) then
      deallocate(model%modes)
    end if
  end subroutine

  !> Extracts a prognostics_t variable from the given C pointer.
  function prognostics_from_c_ptr(ptr) result(retval)
    implicit none
    type(c_ptr), value, intent(in) :: ptr
    type(prognostics_t) :: retval

    retval%ptr = ptr
  end function

  !> Given an aerosol population index p, get the corresponding mode and
  !> aerosol species indices m and s.
  subroutine m_get_mode_and_species(model, p, m, s)
    class(model_t), intent(in)  :: model
    integer, intent(in)         :: p
    integer, intent(out)        :: m, s

    m = 1
    do while (model%population_offsets(m+1) < p)
      m = m + 1
    end do
    if (m == 1) then
      s = p
    else
      s = p - model%population_offsets(m)
    end if
  end subroutine

  !> Given the name of a mode, returns its index within the model.
  !> @param [in] m A pointer to a model object.
  !> @param [in] mode_name The name of the desired mode
  function m_mode_index(m, mode_name) result(mode_index)
    implicit none
    class(model_t),   intent(in) :: m
    character(len=*), intent(in) :: mode_name
    integer :: mode_index

    ! Find the mode index
    do mode_index = 1,m%num_modes
      if (m%modes(mode_index)%name == mode_name) then
        exit
      end if
    end do
  end function

  !> Returns the index of the aerosol species with the given (symbolic) name
  !> within the mode with the given index, or 0 if no such species is found.
  !> @param [in] m A pointer to a model object.
  !> @param [in] mode_index The index of the mode for the desired species
  !> @param [in] species_symbol The abbreviated symbolic name of the species
  function m_aerosol_index(m, mode_index, species_symbol) result(a_index)
    implicit none
    class(model_t),   intent(in) :: m
    integer,          intent(in) :: mode_index
    character(len=*), intent(in) :: species_symbol
    integer :: a_index

    ! If our mode index is invalid, return an invalid aerosol index.
    if (mode_index == 0) then
      a_index = 0
      return
    else
      do a_index = 1,size(m%aero_species, mode_index)
        if (m%aero_species(mode_index, a_index)%symbol == species_symbol) then
          return
        end if
      end do
    end if

    ! No such species
    a_index = 0
  end function

  !> Returns the index of the aerosol population corresponding to the given
  !> mode and aerosol species indices.
  !> @param [in] m A pointer to a model object.
  !> @param [in] mode_index The index of an aerosol mode
  !> @param [in] aero_index The index of an aerosol within the given mode
  function m_population_index(m, mode_index, aero_index) result(pop_index)
    implicit none
    class(model_t),   intent(in) :: m
    integer, intent(in) :: mode_index
    integer, intent(in) :: aero_index
    integer :: pop_index

    pop_index = m%population_offsets(mode_index) + aero_index - 1
  end function

  !> Returns the index of the gas species with the given (symbolic) name, or 0
  !> if no such species is found.
  !> @param [in] m A pointer to a model object.
  !> @param [in] species_symbol The abbreviated symbolic name of the gas species
  function m_gas_index(m, species_symbol) result(g_index)
    implicit none
    class(model_t),   intent(in) :: m
    character(len=*), intent(in) :: species_symbol
    integer :: g_index

    do g_index = 1,size(m%gas_species)
      if (m%gas_species(g_index)%symbol == species_symbol) then
        return
      end if
    end do

    ! No such species
    g_index = 0
  end function

  !> Provides access to the interstitial aerosol mixing fractions array
  !> for the given mode in the given prognostics object.
  !> @param [in] p A pointer to a prognostics object.
  !> @param [in] mode An index identifying the desired mode.
  function p_int_aero_mix_frac(p) result(retval)
    class(prognostics_t), intent(in)  :: p
    real(c_real), pointer, dimension(:,:) :: retval

    type(c_ptr) :: v_ptr
    v_ptr = p_int_aero_mix_frac_c(p%ptr)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels, model%num_populations])
  end function

  !> Provides access to the cloud-borne aerosol mixing fractions array
  !> for the given mode in the given prognostics object.
  !> @param [in] p A pointer to a prognostics object.
  !> @param [in] mode An index identifying the desired mode.
  function p_cld_aero_mix_frac(p) result(retval)
    class(prognostics_t), intent(in)  :: p
    real(c_real), pointer, dimension(:,:) :: retval

    type(c_ptr) :: v_ptr
    v_ptr = p_cld_aero_mix_frac_c(p%ptr)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels, model%num_populations])
  end function

  !> Provides access to the gas mole fractions array for the given
  !> prognostics object.
  !> @param [in] p A Prognostics object.
  function p_gases(p) result(retval)
    use iso_c_binding, only: c_ptr, c_int
    class(prognostics_t), intent(in) :: p
    real(c_real), pointer, dimension(:,:) :: retval

    type(c_ptr) :: v_ptr
    v_ptr = p_gases_c(p%ptr)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels, size(model%gas_species)])
  end function

  !> Provides access to the modal number fractions array for the given
  !> prognostics object.
  !> @param [in] p A Prognostics object.
  function p_modal_num_concs(p) result(retval)
    use iso_c_binding, only: c_ptr, c_int
    class(prognostics_t), intent(in) :: p
    real(c_real), pointer, dimension(:,:) :: retval

    type(c_ptr) :: v_ptr
    v_ptr = p_modal_num_concs_c(p%ptr)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels, size(model%modes)])
  end function

  !> Extracts an atmosphere_t variable from the given C pointer.
  function atmosphere_from_c_ptr(ptr) result(retval)
    implicit none
    type(c_ptr), value, intent(in) :: ptr
    type(atmosphere_t) :: retval

    retval%ptr = ptr
  end function

  !> Provides access to atmosphere temperature column data [K].
  !> @param [in] a A pointer to an atmosphere object.
  function a_temperature(a) result(retval)
    class(atmosphere_t), intent(in)  :: a
    real(c_real), pointer, dimension(:) :: retval

    type(c_ptr) :: v_ptr
    v_ptr = a_temperature_c(a%ptr)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels])
  end function

  !> Provides access to atmosphere pressure column data [Pa].
  !> @param [in] a A pointer to an atmosphere object.
  function a_pressure(a) result(retval)
    class(atmosphere_t), intent(in)  :: a
    real(c_real), pointer, dimension(:) :: retval

    type(c_ptr) :: v_ptr
    v_ptr = a_pressure_c(a%ptr)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels])
  end function

  !> Provides access to atmosphere relative humidity column data [-].
  !> @param [in] a A pointer to an atmosphere object.
  function a_relative_humidity(a) result(retval)
    class(atmosphere_t), intent(in)  :: a
    real(c_real), pointer, dimension(:) :: retval

    type(c_ptr) :: v_ptr
    v_ptr = a_relative_humidity_c(a%ptr)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels])
  end function

  !> Provides access to atmosphere height column data [Pa].
  !> @param [in] a A pointer to an atmosphere object.
  function a_height(a) result(retval)
    class(atmosphere_t), intent(in)  :: a
    real(c_real), pointer, dimension(:) :: retval

    type(c_ptr) :: v_ptr
    v_ptr = a_height_c(a%ptr)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels+1])
  end function

  !> Extracts a diagnostics_t variable from the given C pointer.
  function diagnostics_from_c_ptr(ptr) result(retval)
    implicit none
    type(c_ptr), value, intent(in) :: ptr
    type(diagnostics_t) :: retval

    retval%ptr = ptr
  end function

  !> Returns true if the given diagnostics object contains a (non-modal)
  !> variable with the given name, false otherwise.
  !> @param [in] d A pointer to a diagnostics object.
  !> @param [in] name The name of the desired variable.
  function d_has_var(d, name) result(retval)
    use iso_c_binding, only: c_ptr, c_bool
    class(diagnostics_t), intent(in) :: d
    character(len=*), intent(in) :: name
    integer :: retval

    type(c_ptr) :: c_name
    c_name = f_to_c_string(name)
    retval = d_has_var_c(d%ptr, c_name)
  end function

  !> Provides access to the given (non-modal) variable in the given
  !> diagnostics object.
  !> @param [in] d A pointer to a diagnostics object.
  !> @param [in] name The name of the desired variable.
  function d_var(d, name) result(retval)
    class(diagnostics_t), intent(in)  :: d
    character(len=*), intent(in) :: name
    real(wp), dimension(:), pointer :: retval
    integer(c_int) :: token

    type(c_ptr) :: c_name, v_ptr

    c_name = f_to_c_string(name)
    token = d_has_var_c(d%ptr, c_name)
    v_ptr = d_var_c(d%ptr, token)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels])
  end function

  !> Returns true if the given diagnostics object contains a modal aerosol
  !> variable with the given name, false otherwise.
  !> @param [in] d A pointer to a diagnostics object.
  !> @param [in] name The name of the desired aerosol variable.
  !> @param [in] mode The index of the desired mode.
  function d_has_aerosol_var(d, name) result(retval)
    use iso_c_binding, only: c_ptr, c_char, c_bool, c_int
    class(diagnostics_t), intent(in) :: d
    character(len=*), intent(in) :: name
    integer(c_int) :: retval

    type(c_ptr) :: c_name

    c_name = f_to_c_string(name)
    retval = d_has_aerosol_var_c(d%ptr, c_name)
  end function

  !> Provides access to the given (non-modal) variable in the given
  !> diagnostics object.
  !> @param [in] d A pointer to a diagnostics object.
  !> @param [in] name The name of the desired aerosol variable.
  !> @param [in] mode The index of the desired mode.
  function d_aerosol_var(d, name) result(retval)
    class(diagnostics_t), intent(in)  :: d
    character(len=*), intent(in) :: name
    real(wp), dimension(:,:), pointer :: retval
    integer(c_int) :: token

    type(c_ptr) :: c_name, v_ptr

    c_name = f_to_c_string(name)
    token = d_has_aerosol_var_c(d%ptr, c_name)
    v_ptr = d_aerosol_var_c(d%ptr, token)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels, model%num_populations])
  end function

  !> Returns true if the given diagnostics object contains a (non-modal)
  !> variable with the given name, false otherwise.
  !> @param [in] d A pointer to a diagnostics object.
  !> @param [in] name The name of the desired modal variable.
  function d_has_gas_var(d, name) result(retval)
    use iso_c_binding, only: c_ptr, c_bool
    class(diagnostics_t), intent(in) :: d
    character(len=*), intent(in) :: name
    integer(c_int) :: retval

    type(c_ptr) :: c_name
    c_name = f_to_c_string(name)
    retval = d_has_gas_var_c(d%ptr, c_name)
  end function

  !> Provides access to the given (non-modal) variable in the given
  !> diagnostics object.
  !> @param [in] d A pointer to a diagnostics object.
  !> @param [in] name The name of the desired modal variable.
  function d_gas_var(d, name) result(retval)
    class(diagnostics_t), intent(in)  :: d
    character(len=*), intent(in) :: name
    real(wp), dimension(:,:), pointer :: retval
    integer(c_int) :: token

    type(c_ptr) :: c_name, v_ptr

    c_name = f_to_c_string(name)
    token = d_has_gas_var_c(d%ptr, c_name)
    v_ptr = d_gas_var_c(d%ptr, token)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels, size(model%gas_species)])
  end function

  !> Returns true if the given diagnostics object contains a modal
  !> variable with the given name, false otherwise.
  !> @param [in] d A pointer to a diagnostics object.
  !> @param [in] name The name of the desired modal variable.
  function d_has_modal_var(d, name) result(retval)
    use iso_c_binding, only: c_ptr, c_bool
    class(diagnostics_t), intent(in) :: d
    character(len=*), intent(in) :: name
    integer(c_int) :: retval

    type(c_ptr) :: c_name

    c_name = f_to_c_string(name)
    retval = d_has_modal_var_c(d%ptr, c_name)
  end function

  !> Provides access to the given modal variable in the given
  !> diagnostics object.
  !> @param [in] d A pointer to a diagnostics object.
  !> @param [in] name The name of the desired modal variable.
  function d_modal_var(d, name) result(retval)
    class(diagnostics_t), intent(in)  :: d
    character(len=*), intent(in) :: name
    real(wp), dimension(:,:), pointer :: retval
    integer(c_int) :: token

    type(c_ptr) :: c_name, v_ptr

    c_name = f_to_c_string(name)
    token = d_has_gas_var_c(d%ptr, c_name)
    v_ptr = d_modal_var_c(d%ptr, token)
    call c_f_pointer(v_ptr, retval, [model%num_levels, size(model%modes)])
  end function

  !> Extracts a tendencies_t variable from the given C pointer.
  function tendencies_from_c_ptr(ptr) result(retval)
    implicit none
    type(c_ptr), value, intent(in) :: ptr
    type(tendencies_t) :: retval

    retval%ptr = ptr
  end function

  !> Provides access to the interstitial aerosol mixing fractions array
  !> for the given mode in the given tendencies object.
  !> @param [in] p A pointer to a prognostics object.
  !> @param [in] mode An index identifying the desired mode.
  function t_int_aero_mix_frac(t) result(retval)
    class(tendencies_t), intent(in) :: t
    real(c_real), pointer, dimension(:,:) :: retval

    type(c_ptr) :: v_ptr
    v_ptr = t_int_aero_mix_frac_c(t%ptr)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels, model%num_populations])
  end function

  !> Provides access to the cloud-borne aerosol mixing fractions array
  !> for the given mode in the given tendencies object.
  !> @param [in] t A pointer to a tendencies object.
  !> @param [in] mode An index identifying the desired mode.
  function t_cld_aero_mix_frac(t) result(retval)
    class(tendencies_t), intent(in)  :: t
    real(c_real), pointer, dimension(:,:) :: retval

    type(c_ptr) :: v_ptr
    v_ptr = t_cld_aero_mix_frac_c(t%ptr)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels, model%num_populations])
  end function

  !> Provides access to the gas mole fractions array for the given
  !> tendencies object.
  !> @param [in] p A pointer to a tendencies object.
  function t_gases(t) result(retval)
    use iso_c_binding, only: c_ptr, c_int
    class(tendencies_t), intent(in) :: t
    real(c_real), pointer, dimension(:,:) :: retval

    type(c_ptr) :: v_ptr
    v_ptr = t_gases_c(t%ptr)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels, size(model%gas_species)])
  end function

  !> Provides access to the modal number fractions array for the given
  !> tendencies object.
  !> @param [in] t A pointer to a tendencies object.
  function t_modal_num_concs(t) result(retval)
    use iso_c_binding, only: c_ptr, c_int
    class(tendencies_t), intent(in) :: t
    real(c_real), pointer, dimension(:,:) :: retval

    type(c_ptr) :: v_ptr
    v_ptr = t_modal_num_concs_c(t%ptr)
    call c_f_pointer(v_ptr, retval, shape=[model%num_levels, model%num_modes])
  end function

end module

