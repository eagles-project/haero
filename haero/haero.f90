!> This module contains data structures that allow Fortran modules to access
!> data from the C++ model, prognostics, and diagnostics.
module haero

  use iso_c_binding

  implicit none

  !> Working precision real kind
  integer, parameter :: wp = c_double

  private

  public :: mode_t, species_t, model_t, &
            prognostics_t, diagnostics_t, tendencies_t

  !> This Fortran type is the equivalent of the C++ Mode struct.
  type :: mode_t
    !> Mode name
    character(len=1), dimension(32) :: name
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
    character(len=1), dimension(32) :: name
    !> Species symbol (abbreviation)
    character(len=1), dimension(8)  :: symbol
  end type

  !> This Fortran type is the equivalent of the C++ Model class. Exactly one
  !> read-only instance of a model is available to Fortran processes.
  type :: model_t
    !> The aerosol modes in the model, in indexed order.
    type(mode_t), dimension(:), allocatable :: modes
    !> The aerosol species within each mode. Indexed as (mode, species).
    type(species_t), dimension(:,:), allocatable :: aero_species
    !> The gas species in the model.
    type(species_t), dimension(:), allocatable :: gas_species
    !> The number of columns in the model.
    integer :: num_columns
    !> The number of vertical levels in each column.
    integer :: num_levels
  end type

  !> The resident model instance, available to the single allowable C++ model
  !> instance that uses Fortran-backed processes.
  type(model_t) :: model

  !> This type represents the set of prognostic variables for an aerosol
  !> model.
  type :: prognostics_t
    type(c_ptr) :: ptr
  contains
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
  type :: diagnostics_t
    type(c_ptr) :: ptr
  contains
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

  !> This type represents a set of tendencies to be computed by a prognostic
  !> process.
  type :: tendencies_t
    type(c_ptr) :: ptr
  contains
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

  !> This type represents an interface for prognostic processes. One defines
  !> a prognostic process in Fortran by overridding each of the init, run, and
  !> destroy procedures below.
  type, abstract :: prognostic_process_t
  contains
    procedure(pp_init), deferred    :: init
    procedure(pp_run), deferred     :: run
    procedure(pp_destroy), deferred :: destroy
  end type

  !> This type represents an interface for diagnostic processes. One defines
  !> a diagnostic process in Fortran by overridding each of the init, update,
  !> and destroy procedures below.
  type, abstract :: diagnostic_process_t
  contains
    procedure(dp_init), deferred    :: init
    procedure(dp_update), deferred  :: update
    procedure(dp_destroy), deferred :: destroy
  end type

  interface
    integer(c_int) function p_num_aero_modes_c(p) bind(c)
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: p
    end function

    integer(c_int) function p_num_aero_species_c(p) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
    end function

    integer(c_int) function p_num_gas_species_c(p) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
    end function

    integer(c_int) function p_num_columns_c(p) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
    end function

    integer(c_int) function p_num_levels_c(p) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
    end function

    subroutine p_get_int_aero_mix_frac_c(p, mode, v_ptr, v_shape) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
      integer(c_int), value, intent(in) :: mode
      real(wp), dimension(:,:,:), intent(out) :: v_ptr
      integer(c_int), dimension(3), intent(out) :: v_shape
    end subroutine

    subroutine p_get_cld_aero_mix_frac_c(p, mode, v_ptr, v_shape) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
      integer(c_int), value, intent(in) :: mode
      real(wp), dimension(:,:,:), intent(out) :: v_ptr
      integer(c_int), dimension(3), intent(out) :: v_shape
    end subroutine

    subroutine p_get_gas_mole_frac_c(p, v_ptr, v_shape) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
      real(wp), dimension(:,:,:), intent(out) :: v_ptr
      integer(c_int), dimension(3), intent(out) :: v_shape
    end subroutine

    subroutine p_get_modal_num_densities_c(p, v_ptr, v_shape) bind(c)
      use iso_c_binding, only: c_ptr, c_int
      type(c_ptr), value, intent(in) :: p
      real(wp), dimension(:,:,:), intent(out) :: v_ptr
      integer(c_int), dimension(3), intent(out) :: v_shape
    end subroutine

    subroutine pp_init(process, model)
      import prognostic_process_t
      class(prognostic_process_t), intent(inout) :: process
      type(model_t), intent(in)                  :: model
    end subroutine

    subroutine pp_run(process, model, t, dt, prognostics, diagnostics, tendencies)
      import prognostic_process_t
      class(prognostic_process_t), intent(inout) :: process
      type(model_t), intent(in)                  :: model
      real(wp), intent(in)                       :: t
      real(wp), intent(in)                       :: dt
      type(prognostics_t), intent(in)            :: prognostics
      type(diagnostics_t), intent(in)            :: diagnostics
      type(tendencies_t), intent(inout)          :: tendencies
    end subroutine

    subroutine pp_destroy(process)
      import prognostic_process_t
      class(prognostic_process_t), intent(inout) :: process
    end subroutine

    subroutine dp_init(process, model)
      import diagnostic_process_t
      class(diagnostic_process_t), intent(inout) :: process
      type(model_t), intent(in)                  :: model
    end subroutine

    subroutine dp_run(process, model, t, prognostics, diagnostics)
      import diagnostic_process_t
      class(diagnostic_process_t), intent(inout) :: process
      type(model_t), intent(in)                  :: model
      real(wp), intent(in)                       :: t
      type(prognostics_t), intent(in)            :: prognostics
      type(diagnostics_t), intent(inout)         :: diagnostics
    end subroutine

    subroutine dp_destroy(process)
      import diagnostic_process_t
      class(diagnostic_process_t), intent(inout) :: process
    end subroutine

  end interface

contains

  ! Begin the process of initializing the Haero Fortran module.
  subroutine haerotran_begin_init() bind(c)
    ! Nothing here yet!
  end subroutine

  ! Set the number of modes in the global model.
  subroutine haerotran_set_num_modes(num_modes) bind(c)
    use iso_c_binding, only: c_int
    implicit none

    integer(c_int), value, intent(in) :: num_modes

    allocate(model%modes(num_modes))
  end subroutine

  subroutine haerotran_set_max_modal_aero_species(max_num_species) bind(c)
    use iso_c_binding, only: c_int
    implicit none

    integer(c_int), value, intent(in) :: max_num_species

    allocate(model%aero_species(size(model%modes), max_num_species))
  end subroutine

  subroutine haerotran_set_mode(mode, name, min_d, max_d, std_dev) bind(c)
    use iso_c_binding, only: c_char, c_double
    implicit none

    integer(c_int), intent(in) :: mode
    character(len=1), dimension(:), intent(in) :: name
    real(c_double), value, intent(in) :: min_d
    real(c_double), value, intent(in) :: max_d
    real(c_double), value, intent(in) :: std_dev

    integer :: i

    do i=1, min(32,size(name))
      model%modes(mode)%name(i:i) = name(i)
    end do
    model%modes(mode)%min_diameter = min_d
    model%modes(mode)%max_diameter = max_d
    model%modes(mode)%mean_std_dev = std_dev

  end subroutine

  subroutine haerotran_set_aero_species(mode, species, name, symbol) bind(c)
    use iso_c_binding, only: c_int, c_char
    implicit none

    integer(c_int), intent(in) :: mode
    integer(c_int), intent(in) :: species
    character(len=1), dimension(:), intent(in) :: name
    character(len=1), dimension(:), intent(in) :: symbol

    integer :: i

    do i=1, min(32,size(name))
      model%aero_species(mode, species)%name(i:i) = name(i)
    end do
    do i=1, min(8,size(symbol))
      model%aero_species(mode, species)%symbol(i:i) = symbol(i)
    end do
  end subroutine

  subroutine haerotran_set_num_gas_species(num_species) bind(c)
    use iso_c_binding, only: c_int
    implicit none

    integer(c_int), value, intent(in) :: num_species

    allocate(model%gas_species(num_species))
  end subroutine

  subroutine haerotran_set_gas_species(species, name, symbol) bind(c)
    use iso_c_binding, only: c_char, c_double
    implicit none

    integer(c_int), intent(in) :: species
    character(len=1), dimension(:), intent(in) :: name
    character(len=1), dimension(:), intent(in) :: symbol

    integer :: i

    do i=1, min(32,size(name))
      model%gas_species(species)%name(i:i) = name(i)
    end do
    do i=1, min(8,size(symbol))
      model%gas_species(species)%symbol(i:i) = symbol(i)
    end do
  end subroutine

  subroutine haerotran_set_num_columns(num_columns) bind(c)
    use iso_c_binding, only: c_int
    implicit none

    integer(c_int), value, intent(in) :: num_columns

    model%num_columns = num_columns
  end subroutine

  subroutine haerotran_set_num_levels(num_levels) bind(c)
    use iso_c_binding, only: c_int
    implicit none

    integer(c_int), value, intent(in) :: num_levels

    model%num_levels = num_levels
  end subroutine

  ! Wrap up the process of initializing the Haero Fortran module.
  subroutine haerotran_end_init() bind(c)
    ! Nothing here yet!
  end subroutine

  ! This subroutine gets called when the C++ process exits.
  subroutine haerotran_finalize() bind(c)
    if (allocated(model%modes)) then
      deallocate(model%aero_species)
      deallocate(model%gas_species)
      deallocate(model%modes)
    end if
  end subroutine

  ! Below are interoperable functions used by C++ processes to invoke Fortran
  ! process implementations.
  subroutine fortran_prognostic_process_init(c_process, c_model) bind(c)
    use iso_c_binding, only: c_ptr, c_f_pointer
    implicit none

    ! Arguments: C pointers to process and model data.
    type(c_ptr), intent(inout) :: c_process
    type(c_ptr), intent(in)    :: c_model

    ! Local variables: Fortran data types for arguments.
    type(prognostic_process_t), pointer :: process
    type(model_t)                       :: model

    ! Convert the C data to Fortran data.
    call c_f_pointer(c_process, process)
    model%ptr = c_model

    ! Initialize the process.
    call process%init(model)
  end subroutine

  subroutine fortran_prognostic_process_run(c_process, c_model, t, dt, &
                                            c_progs, c_diags, c_tends) bind(c)
    use iso_c_binding, only: c_ptr, c_f_pointer
    implicit none

    ! Arguments
    type(c_ptr), intent(inout) :: c_process
    type(c_ptr), intent(in)    :: c_model
    real(wp), intent(in)       :: t
    real(wp), intent(in)       :: dt
    type(c_ptr), intent(in)    :: c_progs
    type(c_ptr), intent(in)    :: c_diags
    type(c_ptr), intent(inout) :: c_tends

    ! Local variables: Fortran data types for arguments.
    type(prognostic_process_t), pointer :: process
    type(model_t)                       :: model
    type(prognostics_t)                 :: prognostics
    type(diagnostics_t)                 :: diagnostics
    type(tendencies_t)                  :: tendencies

    ! Convert C data to Fortran data.
    call c_f_pointer(c_process, process)
    model%ptr = c_model
    prognostics%ptr = c_progs
    diagnostics%ptr = c_diags
    tendencies%ptr = c_tends

    ! Run the process.
    call process%run(model, t, dt, prognostics, diagnostics, tendencies)
  end subroutine

  subroutine fortran_prognostic_process_destroy(c_process) bind(c)
    use iso_c_binding, only: c_ptr, c_f_pointer
    implicit none

    ! Argument: C pointers to process data.
    type(c_ptr), intent(inout) :: c_process

    ! Local variable: Fortran process data.
    type(prognostic_process_t), pointer :: process

    ! Convert the C data to Fortran data.
    call c_f_pointer(c_process, process)

    ! Destroy the process.
    call process%destroy()
  end subroutine

  subroutine fortran_diagnostic_process_init(c_process, c_model) bind(c)
    use iso_c_binding, only: c_ptr, c_f_pointer
    implicit none

    ! Arguments: C pointers to process and model data.
    type(c_ptr), intent(inout) :: c_process
    type(c_ptr), intent(in)    :: c_model

    ! Local variables: Fortran data types for arguments.
    type(diagnostic_process_t), pointer :: process
    type(model_t)                       :: model

    ! Convert the C data to Fortran data.
    call c_f_pointer(c_process, process)
    model%ptr = c_model

    ! Initialize the process.
    call process%init(model)
  end subroutine

  subroutine fortran_diagnostic_process_update(c_process, c_model, t, &
                                               c_progs, c_diags) bind(c)
    use iso_c_binding, only: c_ptr, c_f_pointer
    implicit none

    ! Arguments
    type(c_ptr), intent(inout) :: c_process
    type(c_ptr), intent(in)    :: c_model
    real(wp), intent(in)       :: t
    type(c_ptr), intent(in)    :: c_progs
    type(c_ptr), intent(inout) :: c_diags

    ! Local variables: Fortran data types for arguments.
    type(prognostic_process_t), pointer :: process
    type(model_t)                       :: model
    type(prognostics_t)                 :: prognostics
    type(diagnostics_t)                 :: diagnostics

    ! Convert C data to Fortran data.
    call c_f_pointer(c_process, process)
    model%ptr = c_model
    prognostics%ptr = c_progs
    diagnostics%ptr = c_diags

    ! Invoke the process.
    call process%update(model, t, prognostics, diagnostics)
  end subroutine

  subroutine fortran_diagnostic_process_destroy(c_process) bind(c)
    use iso_c_binding, only: c_ptr, c_f_pointer
    implicit none

    ! Argument: C pointers to process data.
    type(c_ptr), intent(inout) :: c_process

    ! Local variable: Fortran process data.
    type(diagnostic_process_t), pointer :: process

    ! Convert the C data to Fortran data.
    call c_f_pointer(c_process, process)

    ! Destroy the process.
    call process%destroy()
  end subroutine

  !> Returns the number of aerosol modes in the given prognostics object.
  !> @param [in] p A Prognostics object.
  integer(c_int) function p_num_aero_modes(p)
    type(prognostics_t), intent(in) :: p
    p_num_aero_modes = p_num_aero_modes_c(p%ptr)
  end function

  !> Returns the number of aerosol species in the given prognostics object.
  !> @param [in] p A Prognostics object.
  integer(c_int) function p_num_aero_species(p)
    type(prognostics_t), intent(in) :: p
    p_num_aero_species = p_num_aero_species_c(p%ptr)
  end function

  !> Returns the number of gas species in the given prognostics object.
  !> @param [in] p A Prognostics object.
  integer(c_int) function p_num_gas_species(p)
    type(prognostics_t), intent(in) :: p
    p_num_gas_species = p_num_gas_species_c(p%ptr)
  end function

  !> Returns the number of columns in the given prognostics object.
  !> @param [in] p A Prognostics object.
  integer(c_int) function p_num_columns(p)
    type(prognostics_t), intent(in) :: p
    p_num_columns = p_num_columns_c(p%ptr)
  end function

  !> Returns the number of vertical levels in the given prognostics object.
  !> @param [in] p A Prognostics object.
  integer(c_int) function p_num_levels(p)
    type(prognostics_t), intent(in) :: p
    p_num_levels = p_num_levels_c(p%ptr)
  end function

  !> Provides access to the interstitial aerosol mixing fractions array
  !> for the given mode in the given prognostics object.
  !> @param [in] p A pointer to a prognostics object.
  !> @param [in] mode An index identifying the desired mode.
  function p_int_aero_mix_frac(p, mode) result(data)
    type(prognostics_t), intent(in)  :: p
    integer(c_int), intent(in) :: mode
    real(wp), pointer, dimension(:,:,:), intent(out) :: data

    type(c_ptr) :: v_ptr
    integer(c_int), dimension(3) :: v_shape
    call p_get_int_aero_mix_frac_c(p%ptr, mode, v_ptr, data_shape)
    call c_f_pointer(v_ptr, data, v_shape)
  end function

  !> Provides access to the cloud-borne aerosol mixing fractions array
  !> for the given mode in the given prognostics object.
  !> @param [in] p A pointer to a prognostics object.
  !> @param [in] mode An index identifying the desired mode.
  function p_cld_aero_mix_frac(p, mode) result(data)
    type(prognostics_t), intent(in)  :: p
    integer(c_int), intent(in) :: mode
    real(wp), pointer, dimension(:,:,:), intent(out) :: data

    type(c_ptr) :: v_ptr
    integer(c_int), dimension(3) :: v_shape
    call p_get_cld_aero_mix_frac_c(p%ptr, mode, v_ptr, data_shape)
    call c_f_pointer(v_ptr, data, v_shape)
  end function

  !> Provides access to the gas mole fractions array for the given
  !> prognostics object.
  !> @param [in] p A Prognostics object.
  function p_gas_mole_frac(p) result(data)
    use iso_c_binding, only: c_ptr, c_int
    type(prognostics_t), value, intent(in) :: p
    real(wp), pointer, dimension(:,:,:), intent(out) :: data

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
    type(prognostics_t), value, intent(in) :: p
    real(wp), pointer, dimension(:,:,:), intent(out) :: data

    type(c_ptr) :: v_ptr
    integer(c_int), dimension(3) :: v_shape
    call p_get_modal_num_densities_c(p%ptr, v_ptr, data_shape)
    call c_f_pointer(v_ptr, data, v_shape)
  end function

end module

