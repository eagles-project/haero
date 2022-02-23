
# The Haero Library

## Overview

Haero is designed to provide a modal aerosol capability to an atmospheric model
written in C++ and/or Fortran. It makes no attempt to describe or evolve any
atmospheric phenomena outside of aerosols. Instead, Haero evolves the state of
aerosols within a modal aerosol model as part of a broader atmospheric
**host model**: a mathematically consistent description of the atmosphere.

To use Haero in your own host model, you write code to interact to construct
a modal aerosol system and invoke aerosol processes on that system. Haero gives
you all the flexibility and control you need to define how the aerosol processes
couple with and interoperate with the other processes in the host model. In this
sense, Haero is a set of building blocks you can use to construct the most
appropriate modal aerosol representation for your host model.

Haero provides an interface for running aerosol processes to evolve a set of
state data within a single atmospheric column. You can evolve several columns
in parallel, as long as state data is provided for each column.

All quantities in Haero are specified using the International System of Units
(SI) unless otherwise specified. In both code comments and documentation, we
place square brackets indicating units after the description of a quantity
(e.g. simulation time [s]).

### Aerosol Systems

Haero's representation of aerosols relies on a set of simple data structures
that define the assumptions underlying a specific modal aerosol system. These
elements are:

* **Modes**: statistical representations of aerosol particle populations
  organized by particle size
* **Species**: aerosol and gas molecules of interest. Each aerosol
  species belongs to a single aerosol mode and is tracked by mass and
  number. Gas species are assumed to be small and don't belong to modes.
* An **Aerosol Configuration**: the collection of modes and species of
  interest within a specific modal aerosol system to be simulated

These entities define the aerosol system of interest in Haero, and provide
any related **metadata** needed to make decisions about how an aerosol
processes does its work.

### Aerosol State Data

Haero deals with two distinct types of state variables:

* **Prognostic variables**: variables that are evolved in time according
  to a system of differential equations
* **Diagnostic variables**: variables that are algebraically related to
  other variables, whether those variables are prognostic or diagnostic

Prognostic variables are quantities that possess an initial state and are
evolved forward in time by their **tendencies** (time derivatives). It is
not possible to construct the value of a prognostic variable at time $t$ without
an initial condition at some time $t_0$ and a tendency defined over the period
$\left[t_0, t\right]$.

The concept of a "diagnostic" variable is more general than its name suggests.
The word *diagnostic* suggests that the variable is used only as an indicator
by a human attempting to "diagnose" some atmospheric condition. In fact, a
diagnostic variable can be any variable whose state can be constructed at any
instant in time, using only the relevant prognostic variables. In this sense,
a diagnostic variable serves as a "shared" variable that can be computed at
the appropriate time, and used by one or more aerosol processes.

Haero's aerosol state data lives in multi-dimensional arrays within "smart
containers":

* The `Prognostics` container contains prognostic state variables
  specific to aerosols
* The `Atmosphere` container contains a thermodynamic description
  of the atmosphere in which an aerosol system is embedded
* The `Diagnostics` container contains a registry of diagnostic
  variables shared amongst various aerosol processes, and made available
  for output

The arrays in these data structures are stored in [Kokkos Views](https://github.com/kokkos/kokkos/wiki/View).
Aerosol state data is allocated in C++, but available for use in Fortran for
implementing aerosol processes or for using Haero from within a Fortran host
model. Atmospheric state data is provided by the host model.

### Aerosol Processes

The aerosol "life cycle" consists of a set of complicated physical processes
involving many participants, with a wide range of length and time scales. These
different scales demand a degree of flexibility in how we evaluate changes to
the state of an aerosol system. For example, we expect to be able to resolve
processes whose time scale is similar to or larger than the time scale for
convection in the atmosphere, whereas processes with faster time scales must be
treated in some special way that accommodates a relaxation or equilibration
process.

In Haero, an "aerosol process" accepts a set of completely specified input
(prognostic and diagnostic variables, an atmospheric state, etc) and calculates
a set of tendencies for prognostic aerosol and gas variables.

An aerosol process consists of a set of parameterizations that encode simplifying
assumptions about a specific stage of the aerosol life cycle into an algorithm
that computes the relevant quantities of interest. The processes provided by
Haero, and their various parameterizations, are described in
[processes](processes.md).

These processes are the true assets of the Haero library. They can be
implemented in C++, in Fortran, or in both. This allows aerosol researchers to
make their latest parameterizations available in the Haero library, while
providing software engineers with a "future-proof" environment for optimizing
these and other parameterizations for DOE's Leadership Class Facilities.

That's an orbit-level view of the Haero library. Now let's take a closer look at
each of these aspects.

## Aerosol Systems in Haero

To model a specific aerosol system in Haero, we must answer some questions:

* How are aerosol particle sizes represented?
* What species are present in the system, and how are their sizes
  reflected in the particle size distribution?
* How is the state of an aerosol system represented mathematically?
* How many vertical levels are needed to resolve the profile of aerosols
  in the system?
* What are the relevent physical processes that describe the evolution
  of the aerosol system?

Each of these decisions greatly affects the nature of the system--systems with
different answers to these questions can have very different behavior. Haero
attempts to answer these questions with a few basic data structures.

## Particle Sizes: the Mode Type

We've seen how the dynamics of aerosols can be represented mathematically
by evolution equations for moments of modal distribution functions. Modes
simplify the description of aerosol particles in terms of their size: instead of
representing a population of particles with a distribution function
$n(V_p, \vec{x}, t)$ that varies continuously with the size of the particle, we
introduced $M$ discrete modes and declared that these modes partition the
population of aerosol particles in the sense of the modal assumption as given
by \refeq{modal_n}.

The essential information in a mode is the range of particle sizes it
encompasses, $[D_{\min}, D_{\max}]$, and its geometric standard deviation,
$\sigma_g$. In Haero's C++ interface, we represent an aerosol mode with the
`Mode` struct, whereas in Fortran we use the `mode_t` derived type:

=== "C++"
    ``` c++
    struct Mode {
      std::string name;  // a unique identifier for the mode
      Real min_diameter; // the mode's minimum particle diameter [m]
      Real max_diameter; // the mode's maximum particle diameter [m]
      Real mean_std_dev; // the geometric mean standard deviation for the mode [m]
    };
    ```
=== "Fortran"
    ``` fortran
    type :: mode_t
      ! Mode name
      character(len=:), allocatable :: name
      ! Minimum particle diameter [m]
      real(wp) :: min_diameter
      ! Maximum particle diameter [m]
      real(wp) :: max_diameter
      ! Geometric mean standard deviation [m]
      real(wp) :: mean_std_dev
    end type
    ```

The log-normal PDF for each mode is defined by two quantities, the (constant)
mean standard deviation contained within the `Mode` class, and the geometric
mean. The geometric mean is a variable function of the mass mixing ratios of all
species contained within the mode and the number mixing ratio of the mode
itself. Since these are time-dependent variables, they are not contained in the
`Mode` class, which is for (invariant) metadata only.

We also note that, mathematically, the log-normal size distributions of each
mode do not have a maximum or minimum; they are bounded by 0 and $\infty$.
The `min_diameter` and `max_diameter` member variables should therefore not be
interpreted in the context of the log-normal functions defined by e.g.,
\eqref{eq:log_normal_pdf_log}. Instead they represent the boundaries beyond
which particles are better described by a different mode. These minimum and
maximum sizes are used to trigger redistribution of particle mass and number
mixing ratios between modes.

In principle, a Haero calculation can support any number of modes, but care
must be taken to ensure that the modal assumptions remain valid, and that
the parametrizations selected can accommodate the given modes.

## Aerosol Species: the Species Type

A particle species is a specifically-identified molecular assembly with
a number of relevant physical properties. The fundamental description of a
species includes

* a descriptive name (e.g. `sulfate`)
* a symbolic name (e.g. `SO4`, for sulfate)
* information about the chemical properties of the species

Each aerosol mode consists of one or more particle species. Additionally,
gas particles also come in different species. Aerosol particles and gas
particles have physical properties that are described respectively by the
`AerosolSpecies` and `GasSpecies` types in C++, and the `aerosol_species_t` and
`gas_species_t` derived types in Fortran.

In Haero, we represent this information in the following way:

=== "C++"
    ``` c++
    struct AerosolSpecies {
      std::string name;          // full species name
      std::string symbol;        // abbreviated symbolic name
      Real molecular weight;     // molecular weight [g/mol]
      Real crystalization_point; // crystalization point [?]
      Real deliquescence_point;  // deliquenscence point [?]
    };

    struct GasSpecies {
      std::string name;          // full species name
      std::string symbol;        // abbreviated symbolic name
      Real molecular weight;     // molecular weight [g/mol]
    };
    ```

=== "Fortran"
    ``` fortran
    type :: aerosol_species_t
      ! Species name
      character(len=:), allocatable :: name
      ! Species symbol (abbreviation)
      character(len=:), allocatable :: symbol
      ! Molecular weight [g/mol]
      real(wp) :: molecular_wt
      ! Crystalization point [?]
      real(wp) :: crystal_pt
      ! Deliquenscence point [?]
      real(wp) :: deliques_pt
    end type

    type :: gas_species_t
      ! Species name
      character(len=:), allocatable :: name
      ! Species symbol (abbreviation)
      character(len=:), allocatable :: symbol
      ! Molecular weight [g/mol]
      real(wp) :: molecular_wt
    end type
    ```

## Species and Their Sizes: the Modal Aerosol Configuration Type

We have data types that express particle sizes and particle species. Now we
need something that relates these two pieces of information. In other words,
we need a way to express how particles of a specific aerosol species are allowed
to grow and shrink, and how that activity is reflected in our representation
of particle sizes.

In the past, aerosol models have often elected to fix the modes and aerosol/gas
species that they treat, in order to simplify code development. Haero instead
allows a researcher to select these modes and species at runtime, allowing a
far greater family of aerosol models to be represented.

=== "C++"
    ``` c++
    class ModalAerosolConfig final {
      public:
      // Constructor -- creates a new modal aerosol configuration given all relevant
      // data.
      ModalAerosolConfig(const std::vector<Mode>& aerosol_modes,
                         const std::vector<AerosolSpecies>& aerosol_species,
                         const std::map<std::string, std::vector<std::string> >& mode_species,
                         const std::vector<AerosolSpecies>& gas_species);

      // The list of aerosol modes.
      std::vector<Mode> aerosol_modes;

      // The list of aerosol species.
      std::vector<Species> aerosol_species;

      // The list of gas species.
      std::vector<Species> gas_species;

      // The total number of distinct aerosol species populations in the
      // system, counting appearances of one species in different modes separately.
      int num_aerosol_populations;

      // Returns the list of aerosol species associated with the system with the
      // given mode index.
      std::vector<Species> aerosol_species_for_mode(int mode_index) const;
    };
    ```
=== "Fortran"
    ``` fortran
    type :: modal_aerosol_config_t
      ! The aerosol modes in the model, in indexed order.
      type(mode_t), dimension(:), allocatable :: aerosol_modes
      ! The number of modes in the model. Equal to size(aerosol_modes).
      integer :: num_aerosol_modes
      ! The number of actual species that exist within each mode.
      integer, dimension(:), allocatable :: num_mode_species
      ! population index offsets for modes.
      integer, dimension(:), allocatable :: population_offsets
      ! The total number of distinct aerosol populations.
      integer :: num_aerosol_populations
      ! The aerosol species within each mode. Indexed as (mode, species).
      type(aerosol_species_t), dimension(:,:), allocatable :: aerosol_species
      ! The gas species in the model.
      type(gas_species_t), dimension(:), allocatable :: gas_species
      ! The number of gases in the model. Equal to size(gas_species).
      integer :: num_gases
    contains
      ! Returns the maximum number of aerosol species found in any aerosol mode.
      procedure :: max_species_per_mode => m_max_species_per_mode
      ! Given the index of an aerosol population, retrieve its mode and
      ! (modal) species indices.
      procedure :: get_mode_and_species => m_get_mode_and_species
      ! Given the name of a mode, retrieve its index.
      procedure :: aerosol_mode_index => m_aerosol_mode_index
      ! Given a mode index and the symbolic name of an aerosol species, retrieve
      ! its index within that mode
      procedure :: aerosol_species_index => m_aerosol_species_index
      ! Given mode and aerosol species indices, retrieve a population index
      ! that can be used to access aerosol data.
      procedure :: population_index => m_population_index
      ! Given the symbolic name of a gas, retrieve its index.
      procedure :: gas_index => m_gas_index
    end type
    ```

Once you have a modal aerosol configuration, you can answer the first two
questions at the beginning of this section. Next, we look at how the state of
an aerosol system is represented mathematically.

## Aerosol and Atmospheric State: Container Types

Once you've described the constituents of your aerosol system with a modal
aerosol configuration, you can create state variables for that system. The state
of an aerosol system is defined by the following prognostic variables within the
`Prognostics` data structure:

* **aerosol modal mass mix fraction** $q_{m,s}$:
  the ratio of aerosol mass to dry air mass for aerosol species $s$
  occupying mode $m$ $q_{m,s}$ [kg aerosol species $s$ /kg dry air]
* **gas mass mix fraction** $q_g$: the ratio of the mass of gas species $g$ to
  dry air mass [kg gas species $g$/kg dry air]
* **modal number concentrations** $n_m$: the total number of particles per unit
  mass of dry air in the mode $m$ [\# /kg dry air]

The state of the atmosphere (expressed in averaged thermodynamic quantities
like pressure and temperature) greatly affects the behavior of aerosols, so
this atmosphere state information is made available in the `Atmosphere` data
structure.

Finally, the system can use a set of diagnostic variables, stored in the
`Diagnostics` data structure, that depend on a set of aerosol processes
(which are discussed in a later section).

The `Prognostics`, `Atmosphere`, and `Diagnostics` containers store state data
in multidimensional arrays allocated in C++ but made available to both C++ and
Fortran. The data for each array is stored within a Kokkos `View`.

### Digression: Kokkos Views as Multidimensional Arrays

The C++ programming language has lots of features, but remarkably it includes no
mechanism for allocating multidimensional arrays at runtime. The Kokkos C++
library fills this gap by providing a data structure called a `View`. A
`View` is essentially an interface that allows a C++ programmer to treat a
chunk of memory like a multidimensional array.

A `View` has a rank and a set of dimensions, just like an allocatable
Fortran array. You access a `View` in the same way that you'd access a
Fortran array, except that

1. A Kokkos `View` uses row-major indexing instead of Fortran's column-major
   indexing
2. A Kokkos `View` uses 0-based indexing instead of Fortran's 1-based indexing

So for a rank-3 view `f` that you access in C++ as `f(i,j,k)`, you would access
the corresponding array element in Fortran as `f(k-1,j-1,i-1)`. Clear as mud?
Welcome to mixed language development!

#### Packs and Vectorization}

Haero uses Views that consist of `Pack` objects instead of floating point
numbers. A `Pack` (or just "pack" is a contiguous array of numbers that allows a
modern CPU or GPU to vectorize calculations using special instructions.

When executing a vector instruction, a processor performs arithmetic on more
than one number at a time within a mathematical expression. In many cases,
vectorizing expressions can produce significantly faster code. The cost of this
optimization is that a pack represents several numbers, not one. This can make
it tricky to reason about the physical quantities stored in a pack.

Haero uses packs with a size (number of contiguously stored numbers) set at
compile time by the CMake variable `HAERO_PACK_SIZE`. To simplify the
process of reasoning about packs, Haero uses these objects one way only:
in Haero a pack stores data for `HAERO_PACK_SIZE` vertical levels in
a column. Thus, a pack contains data for exactly one variable (with the same
units and physical interpretation) whose values span one or more vertical
levels.

This is the easiest way for Haero to support vectorization. It does mean,
however, that the number of vertical levels in a column differs in general from
the number of packs spanning a vertical level. For example, a column of data
with 72 vertical levels running in a Haero build with a `HAERO_PACK_SIZE`
of 2 contains $36 = 72 / 2 $ packs in its vertical extent.

#### Haero-Specific Views

Because Haero is concerned with arrays having very specific dimensions, we
define some named types that correspond to views/arrays that span specific
spaces:

| View Name           | Rank | Description | C++ | Fortran |
| ------------------- | ---- | ----------- | --- | ------- |
| `ColumnView`        |    1 | Maps a vertical level index $k$ to a pack | `v(k)` | `v(k)` |
| `SpeciesColumnView` |    2 | Maps a population index $p$ and a vertical level index $k$ to a pack | `v(p,k)` | `v(k,p)` |
| `ModeColumnView`    |    2 | Maps a mode index $m$ and a vertical level index $k$ to a pack | `v(m,k)` | `v(k,m)` |

The `Prognostics`, `Atmosphere`, and `Diagnostics` containers described below
make use of these named types.

### Prognostics Type

The `Prognostics` type provides access to prognostic variables that
describe aerosols in a modal description. Here's the essential information for
the C++ and Fortran interfaces (abbreviated for brevity---see the full
interfaces in `haero/prognostics.hpp` and `haero/haero.F90`):

=== "C++"
    ``` c++
    class Prognostics final {
      public:
      // Returns the number of aerosol modes in the system.
      int num_aerosol_modes() const;

      // Returns the number of aerosol species in the mode with the given index.
      int num_aerosol_species(int mode_index) const;

      // Returns the number of gas species in the system.
      int num_gas_species() const;

      // Returns the number of vertical levels in the system.
      int num_levels() const;

      // Returns the view storing interstitial aerosol species mass mixing fraction
      // data.
      const SpeciesColumnView& interstitial_aerosols() const;

      // Returns the view storing cloud-borne aerosol species mass mixing fraction
      // data.
      const SpeciesColumnView& cloudborne_aerosols() const;

      // Returns the view storing mass mixing fraction data for gas species.
      const SpeciesColumnView& gases() const;

      // Returns the view storing modal number concentrations.
      const ModeColumnView& modal_num_concs() const;

      // Scales the given set of tendencies and adds it into this state, summing
      // the values of the prognostic variables in place.
      void scale_and_add(Real scale_factor, const Tendencies& tendencies);
    };
    ```
=== "Fortran"
    ``` fortran
    type :: prognostics_t
    contains
      ! Access to interstitial aerosol mix fractions array (no dummy arguments)
      procedure :: interstitial_aerosols => p_int_aero_mix_frac
      ! Access to cloudborne aerosol mix fractions array (no dummy arguments)
      procedure :: cloudborne_aerosols => p_cld_aero_mix_frac
      ! Access to gas mix fractions array (no dummy arguments)
      procedure :: gases => p_gas_mix_frac
      ! Access to modal number concentrations array (no dummy arguments)
      procedure :: modes => p_modal_num_concs
    end type
    ```

Typically, you never modify a `Prognostics` variable directly. Instead, you
compute a set of tendencies in a `Tendencies` variable and accumulate them into
your `Prognostics` variable by calling `scale_and_add`.

### Atmosphere Type

The `Atmosphere` type stores a fixed set of state variables that describe the
atmosphere, such as

* temperature [K]
* pressure [Pa]
* relative humidity [-]
* heights at level interfaces [m]

Each of these variables are stored in `ColumnView` objects whose memory
is managed by the host model. Here's how the interfaces look:

=== "C++"
    ``` c++
    ```
=== "Fortran"
    ``` fortran
    ```

### Diagnostics Type

The `Diagnostics` type stores a dynamically-determined set of diagnostic
variables that correspond to the specific parameterizations available
to a specific aerosol system. The variables are identified by unique tokens that
can be retrieved by name.

=== "C++"
    ``` c++
    class Diagnostics final {
      public:
      // Returns the number of aerosol modes in the system.
      int num_aerosol_modes() const;

      // Returns the number of aerosol species in the mode with the given index.
      int num_aerosol_species(int mode_index) const;

      // Returns the number of gas species in the system.
      int num_gas_species() const;

      // Returns the number of vertical levels in the system.
      int num_levels() const;

      // Returns a unique token that identifies the given (non-modal) variable
      // within this object. Returns VAR_NOT_FOUND if this variable does not exist.
      Token find_var(const std::string& name) const;

      // Returns the view storing the diagnostic variable with a name corresponding
      // to the given token. If such a variable does not exist, this throws an
      // exception.
      ColumnView& var(const Token token);

      // Returns a unique token that identifies the given modal aerosol variable
      // within this object. Returns VAR_NOT_FOUND if this variable does not exist.
      Token find_aerosol_var(const std::string& name) const;

      // Returns the view storing the modal aerosol diagnostic variable with a name
      // corresponding to the given token. If such a variable does not exist, this
      // throws an exception.
      SpeciesColumnView& aerosol_var(const Token token);

      // Returns a unique token that identifies the given gas variable within this
      // object. Returns VAR_NOT_FOUND if this variable does not exist.
      Token find_gas_var(const std::string& name) const;

      // Returns the view storing the gas diagnostic variable with a name
      // corresponding to the given token. If such a variable does not exist, this
      // throws an exception.
      SpeciesColumnView& gas_var(const Token token);

      // Returns a unique token that identifies the given modal variable within
      // this object. Returns VAR_NOT_FOUND if this variable does not exist.
      Token find_modal_var(const std::string& name) const;

      // Returns the view storing the mode-specific diagnostic variable with a name
      // corresponding to the given token. If such a variable does not exist, this
      // throws an exception.
      ModeColumnView& modal_var(const Token token);
    };
    ```
=== "Fortran"
    ``` fortran
    type :: diagnostics_t
    contains
      ! Returns a token that can be used to retrieve a variable with the given
      ! name from a diagnostics object, or var_not_found (-1) if no such variable
      ! exists.
      procedure :: find_var(name) -> token
      ! Provides access to the given (non-modal) variable in the given
      ! diagnostics object, given its token
      procedure :: var(token) -> array pointer
      ! Returns a token that can be used to retrieve an aerosol variable with the
      ! given name and mode from a diagnostics object, or var_not_found (-1) if no such
      ! variable exists.
      procedure :: find_aerosol_var(name) -> token
      ! Provides access to the given (non-modal) variable in the given
      ! diagnostics object, given its token.
      procedure :: aerosol_var(token) -> array pointer
      ! Returns a token that can be used to retrieve a gas variable with the
      ! given name from a diagnostics object, or var_not_found (-1) if no such
      ! variable exists.
      procedure :: find_gas_var(name) -> token
      ! Provides access to the given gas variable in the given diagnostics object,
      ! given its token.
      procedure :: gas_var(token) -> array pointer
      ! Returns a token that can be used to retrieve a modal variable with the
      ! given name from a diagnostics object, or var_not_found (-1) if no such
      ! variable exists.
      procedure :: has_modal_var(name) -> token
      ! Provides access to the given modal variable in the given
      ! diagnostics object, given its token.
      procedure :: modal_var(token) -> array pointer
    end type
    ```

At this point, you might wonder how a `Diagnostics` variable knows which
variables it needs. In fact, the `Diagnostics` type provides functions
for creating variables that it needs when it needs them.

For examples of how the `Prognostics`, `Atmosphere`, and `Diagnostics` types
are used in practice, take a look at one of the existing aerosol process
implementations.

## Aerosol Processes in Haero

The aerosol life cycle consists of several important and distinct physical
processes. Haero offers a data structures that makes it very easy to implement
such a process. Because the structure of a given process doesn't depend on the
details of its implementation, we can define an abstract interface to simplify
its implementation. Instead of designing a new process from the ground up every
time you want to add new functionality to Haero, you can simply implement a
small number of functions (or subroutines) that define the behavior of a
process, and let the Haero library handle the details of how these processes are
created and used (and where they run).

For detailed descriptions of the specific processes provided by Haero, take a
look at the [Aerosol Processes](processes.md) section. You can find examples of
source code for Haero's processes in the `haero/processes` subdirectory.

### The Aerosol Process Interface

An aerosol process has three behaviors which must be defined by any
implementation. Each of these behaviors is implemented in a C++ function or
a Fortran subroutine.

* **initialization**: the process must be able to allocate any resources
  it needs to do its work. These resources include temporary work arrays,
  look-up tables, and quantities that need to be precomputed. State data
  is not managed by processes, so it's not included in process
  initialization. If nothing needs to be done for initialization, its
  function or subroutine body can be empty.
* **running**: the process must know how to "run". In other words, it
  must define a procedure for computing tendencies for a relevant set of
  prognostic variables given their current values at a specific simulation time,
  along with the current values of any diagnostic variables. The function or
  subroutine that implements this behavior does not apply these tendencies to
  any prognostic variables---it simply computes the tendencies and returns.
* **finalization**: at the end of a simulation program, when the aerosol
  system is destroyed, the process must free all resources it allocated
  in its initialization. If no resources are allocated, the function or
  subroutine body implementing finalization can be empty.

In addition, the process may support named parameters that can be set to
specific values. Some examples of these kinds of parameters are

* Integer-valued parameters that select one of several supported algorithms
* Boolean flags for enabling or disabling features
* Real-valued scale factors for quantities based on tuning or assumptions
* String-valued parameters (just in case they're helpful)

Haero provides an object-oriented approach for implementing a process in terms
of this simple interface. In an object-oriented approach, an abstract interface
is encoded in a "base class"--a data type that declares the necessary
functions and subroutines. Then any implementation of this interface is defined
in a *derived class*: a type derived from that base class.

Haero uses the object-oriented features of C++ for process development. All
aerosol process implementations are derived from a C++ base class called
`AerosolProcess`. This is true regardless of whether you implement the
process using C++ or Fortran.

The `AerosolProcess` provides the following interface (see
`haero/aerosol_process.hpp` for more details):

=== "AerosolProcess"
    ``` c++
    class AerosolProcess {
     public:

      // Constructor, called by all AerosolProcess subclasses.
      explicit AerosolProcess(const std::string& name):
       name_(name) {}

      // Destructor.
      virtual ~AerosolProcess() {}

      // Initializes the process with the aerosol configuration.
      void init(const ModalAerosolConfig& config);

      // Runs the process at the given time with the given aerosol data and the
      // given Kokkos thread team.
      void run(const TeamType& team, Real t, Real dt,
               const Prognostics& prognostics, const Atmosphere& atmosphere,
               const Diagnostics& diagnostics, Tendencies& tendencies) const;

      // Set named integer, boolean, and real-valued parameters.
      void set_param(const std::string& name, int value);
      void set_param(const std::string& name, bool value);
      void set_param(const std::string& name, Real value);

      // On host: copies this aerosol process to the device, returning a
      // pointer to the copy.
      AerosolProcess* copy_to_device() const;

      // On host: call this static method to delete a copy of the process
      // that has been created on a device.
      static void delete_on_device(AerosolProcess* device_process);

     protected:

      // Override this method if your aerosol process needs to be initialized
      // with information about the system. The default implementation does nothing.
      virtual void init_(const ModalAerosolConfig& config) {}

      // Override this method to implement the aerosol process using the specific
      // parameterization for the subclass.
      virtual void run_(Real t, Real dt,
                        const Prognostics& prognostics,
                        const Atmosphere& atmosphere,
                        const Diagnostics& diagnostics,
                        Tendencies& tendencies) const = 0;

      // Override these methods to set a parameter to a given value based on its
      // name.
      virtual void set_param_(const std::string& name, int value) {}
      virtual void set_param_(const std::string& name, bool value) {}
      virtual void set_param_(const std::string& name, Real value) {}
    };
    ```

In addition to the "constructor" function used to create an instance of a
`AerosolProcess` and the interface functions for initializing, running, and
setting parameters, the interface declares `protected` methods with
underscores after their names. These are the methods you must override in
order to define the behaviors of the aerosol process. Make sure you document
any supported parameters (recognized by your `set_param` methods).

The constructor accepts a single argument: a string containing the name of the
aerosol process. This can be helpful for debugging.

#### Digression: running aerosol processes on a GPU

Haero is designed to allow aerosol physics to be computed on CPUs or GPUs, with
different levels of parallelism. Running code on a GPU is tricky, because the
data it uses must be copied to memory allocated on the GPU itself. In fact,
a process object *itself* must be allocated on the GPU in order for the code
to run there.

The process of allocating this memory on the GPU is esoteric and
confusing. Haero solves this problem by inserting an intermediary class between
your derived class and the `AerosolProcess` class. This intermediary class is
named `DeviceAerosolProcess`. It uses C++'s [curiously recurring template pattern](https://en.cppreference.com/w/cpp/language/crtp)
to add all the necessary logic for your class to run on a GPU.

The way it works is this:

* An object for your process class is allocated and initialized (via `init`) on
  the CPU by an atmospheric host model.
* The host model calls the `copy_to_device` method to obtain a copy of the
  object that lives on the GPU. This method uses a copy constructor defined by
  your process class to copy itself from the CPU to the GPU.
* The host model invokes your GPU-resident object's `run` method within a
  Kokkos parallel dispatch as needed, passing it a [Kokkos thread team](https://github.com/kokkos/kokkos/wiki/HierarchicalParallelism#82-thread-teams)
  that determines the number of threads available witin the method.
* When the calculation is finished, the host model calls the `delete_on_device`
  static method, passing it the GPU-resident object to deallocate it from the
  GPU.

The most important thing to remember here is that `init` is called on the CPU,
whether or not you intend to run your process on the GPU. You must use the
`init` method to record any information from the `ModalAerosolConfig` object
that defines your simulation, because `ModalAerosolConfig` variables cannot
reside on the GPU. Parameters are also set (using `set_param`) on the CPU, not
the GPU.

Let's explore how we might implement an aerosol process in C++ and in Fortran.
Here we describe only the steps needed to implement the process itself. You must
also test your process to make sure it behaves the way you think it does! The
procedure for testing an aerosol process is described in the [Testing](testing.md)
section.

### C++ aerosol processes

In C++, all you have to do in order to implement an aerosol process is to define
a class with a specific name. For concreteness, let's examine a process named
`SimpleNucleationProcess`, which lives in `haero/processes/simple_nucleation_process.hpp`
and `haero/processes/simple_nucleation_process.cpp`.

To allow your process to reside on a CPU or GPU, your derive your class from
the `DeviceAerosolProcess`. This class accepts a single template parameter:
your class. So in our example, you would derive `SimpleNucleationProcess` from
`DeviceAerosolProcess<SimpleNucleationProcess>`. This curiously recursive trick,
in which the type of the intermediary class depends on the type of its
descendent, gives the Curiously Recurring Template Pattern its name. Let's not
worry about how it works for now.

Before we go any further, some terminology: a C++ class derived from a base
class is called a **subclass** of that base class. So `MyProcess` is a
subclass of `DeviceAerosolProcess<SimpleNucleationProcess>`, which is itself a
subclass of the base class `AerosolProcess`.

To create a C++ implemention for an aerosol process:

1. Create a header file that declares your subclass. This header file must
   declare a class constructor, a copy constructor, a destructor, and the
   overridable `init_` and `run_` functions. If you want to support
   configurable parameters, declare whatever versions of `set_param_` you need. Implement your
2. Implement your copy constructor, your destructor, and your `run_` method
   in the header file, declaring each with the `KOKKOS_INLINE_FUNCTION` macro.
   This allows them to be called on a GPU.
3. Create a source file containing implementations for the remaining methods
   (the constructor, the `init_` method, and any `set_param_` methods you need).
   See `haero/processes/simple_nucleation_process.cpp`, for example.
4. Add your source file to the set of source files in the `PROCESS_SOURCES`
   variable in `haero/processes/CMakeLists.txt`.
5. Write one or more tests for your new aerosol process. The [Testing](testing.md)
   section provides details about how to do this.

### Fortran aerosol processes

A Fortran aerosol process implementation consists of a Fortran module that
contains `init`, `run`, and `finalize` subroutines that
implement the same functionality as their C++ counterparts:

=== "Fortran aerosol module"
    ``` fortran
    module MODULE_NAME
      use haero, only: wp, modal_aerosol_config_t, prognostics_t, atmosphere_t, &
                       diagnostics_t, tendencies_t
      ...
      implicit none

      ! Aerosol process interface subroutines
      public :: init, run, finalize

      ! Parameter setting subroutines
      public :: set_integer_param, set_logical_param, set_real_param

      ! Module variables, including settable parameters
      integer :: my_option
      logical :: my_flag
      real(wp) :: my_scale_factor

      ! SAVE keyword for retaining module variables
      save

    contains

    subroutine init(config)
      implicit none

      ! Arguments
      type(modal_aerosol_config_t), intent(in) :: config

      ...
    end subroutine

    subroutine run(t, dt, prognostics, atmosphere, diagnostics, tendencies)
      implicit none

      ! Arguments
      real(wp), value, intent(in)       :: t
      real(wp), value, intent(in)       :: dt
      type(prognostics_t), intent(in)   :: prognostics
      type(atmosphere_t), intent(in)    :: atmosphere
      type(diagnostics_t), intent(in)   :: diagnostics
      type(tendencies_t), intent(inout) :: tendencies

      ...
    end subroutine

    subroutine finalize()
      implicit none

      ...
    end subroutine

    subroutine set_integer_param(name, val)
      implicit none

      ! Arguments
      character(len=*), intent(in) :: name
      integer, intent(in)          :: val

      if (trim(name) == "my_option") then
        my_option = val
      end if
    end subroutine

    subroutine set_logical_param(name, val)
      implicit none

      ! Arguments
      character(len=*), intent(in) :: name
      logical, intent(in)          :: val

      if (trim(name) == "my_flag") then
        my_flag = val
      end if
    end subroutine

    subroutine set_real_param(name, val)
      implicit none

      ! Arguments
      character(len=*), intent(in) :: name
      real(wp), intent(in)         :: val

      if (trim(name) == "my_scale_factor") then
        my_scale_factor = val
      end if
    end subroutine

    ...

    end module
    ```

In Fortran, you **must** implement subroutines for `set_integer_param`,
`set_logical_param`, and `set_real_param`, even if your process
doesn't support settable parameters.

To implement an aerosol process in Fortran, you create such a module in a
Fortran source file and then declare it as a `faerosol_process` in the
`CMakeLists.txt` file within the `haero/processes` subdirectory.
There are instructions on how to declare your process in that file.

For example, the Fortran MAM nucleation process is implemented in a Fortran
module named `mam_nucleation`, implemented in the source file
`haero/processes/mam_nucleation.F90`. In
`haero/processes/CMakeLists.txt`, it's declared as a Fortran aerosol
process the following way:

=== "haero/processes/CMakeLists.txt"
    ```
    if (HAERO_FORTRAN)
      ...
      faerosol_process(MAMNucleationFProcess NucleationProcess mam_nucleation)
      ...
    endif()
    ```

The declaration takes three arguments:

* The name of a C++ class that will be created that exposes the Fortran
  implementation of the process
* The C++ enumerated type that identifies what kind of aerosol process
  is being implemented
* The name of the Fortran module that implements the process

The build system automatically generates a C++ class for the Fortran module that
allows the process to be used in a Haero simulation. The name of this C++ class
is important for exposing it to the Haero library.

Once you've done these things and rebuilt Haero, your new aerosol process
implementation is available for use.

### Diagnostic Functions

Aerosol processes compute tendencies for prognostic variables. But how are
diagnostic variables updated? These variables are quite different in nature from
their prognostic counterparts:

* they depend algebraically on prognostic variables and other diagnostic
  variables
* they are updated in place instead of being evolved in time by
  differential equations
* they are often shared/needed by several distinct aerosol processes, and
  at various points in time, depending on a given process ordering

Because of these considerations, it's not clear that we can update more than a
single diagnostic variable at once. To do so implies a knowledge of how the
aerosol processes are invoked during a time step, and Haero does not make any
such decision on behalf of a host model. This means we must provide a diagnostic
variable update mechanism that is flexible but easy to understand. This
mechanism, which updates a single diagnostic variable, is called a
**diagnostic function**.

Unlike aerosol processes, which have multiple behaviors associated with
initialization, finalization, and the execution of the process itself, a
diagnostic function only ever does one thing: it takes a set of input and
uses it to update its diagnostic variable. In other words, it has no internal
state of its own---it's just a function that you can call whenever you need
to update a particular diagnostic variable.



