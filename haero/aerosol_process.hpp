#ifndef HAERO_AEROSOL_PROCESS_HPP
#define HAERO_AEROSOL_PROCESS_HPP
#include <memory>

#include "haero/atmosphere.hpp"
#include "haero/diagnostics.hpp"
#include "haero/modal_aerosol_config.hpp"
#include "haero/prognostics.hpp"
#include "haero/region_of_validity.hpp"
#include "haero/tendencies.hpp"

namespace haero {

/// @enum AerosolProcessType
/// This enumerated type lists all relevant physical aerosol processes.
enum AerosolProcessType {
  ActivationProcess,
  CloudBorneWetRemovalProcess,
  CoagulationProcess,
  CondensationProcess,
  DryDepositionProcess,
  EmissionsProcess,
  InterstitialWetRemovalProcess,
  NucleationProcess,
  GasAerosolExchangeProcess,
  CalcsizeProcess,
  RenameProcess,
  ResuspensionProcess,
  Terminator  // This isn't a real process type--it only terminates the enum!
};

/// @class AerosolProcess
/// This type is an interface (base class) to an aerosol process quantified by
/// a parametrization that computes a set of tendencies for prognostic variables
/// within an aerosol system. Each subclass of this type implements a particular
/// **implementation** of a specific **parametrization** for a particular
/// **process**.
///
/// To make these ideas more complete, consider the following examples of
/// important **physical processes** in the ProcessType above.
///
/// Each of these processes has one or more **parametrizations**--mathematical
/// models that quantify the outcomes of these processes in specific
/// circumstances. For example, **surface emissions** may be parametrized by:
/// * time series input from a file
/// * a low-order polynomial in the time variable
/// * estimated from, e.g., an agent-based representation of human activity.
///
/// Finally, each of these parametrizations can have one or more
/// **implementations**. For example, every parametrization can have a Fortran
/// implementation that runs only on CPUs, as well as a C++ implementation that
/// can run on CPUs and GPUs.
///
/// The AerosolProcess class provides an interface for all
/// implementations of all parametrizations for all physical processes that
/// compute tendencies for aerosol systems.
class AerosolProcess {
 public:
  /// Constructor, called by all AerosolProcess subclasses.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  AerosolProcess(AerosolProcessType type, const std::string& name)
      : type_(type), name_(name), validity_region_() {}

  /// Destructor.
  KOKKOS_INLINE_FUNCTION
  virtual ~AerosolProcess() {}

  /// Default constructor is disabled.
  AerosolProcess() = delete;

  /// Copy constructor, for use in moving host instance to device.
  KOKKOS_INLINE_FUNCTION
  AerosolProcess(const AerosolProcess& pp)
      : type_(pp.type_),
        name_(pp.name_),
        validity_region_(pp.validity_region_) {}

  /// AerosolProcess objects are not assignable.
  AerosolProcess& operator=(const AerosolProcess&) = delete;

  //------------------------------------------------------------------------
  //                                Accessors
  //------------------------------------------------------------------------

  /// Returns the type of physical process modeled by this process.
  AerosolProcessType type() const { return type_; }

  /// Returns the name of this process/parametrization/realization.
  std::string name() const { return name_.label(); }

  /// Returns the region of validity for this aerosol process.
  const RegionOfValidity& region_of_validity() const {
    return validity_region_;
  }

  /// Returns the region of validity for this aerosol process (non-const).
  RegionOfValidity& region_of_validity() { return validity_region_; }

  //------------------------------------------------------------------------
  //                            Public Interface
  //------------------------------------------------------------------------

  /// On host: performs any system-specific process initialization.
  /// @param [in] config The aerosol configuration describing the aerosol
  ///                    system to which this process belongs.
  void init(const ModalAerosolConfig& config) {
    // This method must be called on the host.
    init_(config);
  }

  /// On device: Validates input aerosol and atmosphere data, returning true if
  /// all data falls within this process's region of validity, and false if not.
  /// @param [in] prognostics The prognostic variables used by and affected by
  ///                         this process.
  /// @param [in] atmosphere The atmosphere state variables used by this
  ///                        process.
  KOKKOS_INLINE_FUNCTION
  bool validate(const Prognostics& prognostics,
                const Atmosphere& atmosphere) const {
    return validity_region_.contains(atmosphere, prognostics);
  }

  /// On device: runs the aerosol process at a given time with the given data.
  /// @param [in] team The Kokkos team used to run this process in a parallel
  ///                  dispatch.
  /// @param [in] t The simulation time at which this process is being invoked
  ///               (in seconds).
  /// @param [in] dt The simulation time interval ("timestep size") over which
  ///                this process occurs.
  /// @param [in] prognostics The prognostic variables used by and affected by
  ///                         this process.
  /// @param [in] atmosphere The atmosphere state variables used by this
  ///                        process.
  /// @param [in] diagnostics The prognostic variables used by and affected by
  ///                         this process.
  /// @param [out] tendencies A container that stores time derivatives for
  ///                         prognostic variables evolved by this process.
  KOKKOS_FUNCTION
  void run(const TeamType& team, Real t, Real dt,
           const Prognostics& prognostics, const Atmosphere& atmosphere,
           const Diagnostics& diagnostics, const Tendencies& tendencies) const {
    // This method must be called on the device.
    run_(team, t, dt, prognostics, atmosphere, diagnostics, tendencies);
  }

  /// On host: Sets a named integer value for this aerosol process.
  /// @param [in] name The name of the parameter to set
  /// @param [in] value The parameter's value
  void set_param(const std::string& name, int value) {
    // This method must be called on the host.
    set_param_(name, value);
  }

  /// On host: Sets a named boolean value for this aerosol process.
  /// @param [in] name The name of the parameter to set
  /// @param [in] value The parameter's value
  void set_param(const std::string& name, bool value) {
    // This method must be called on the host.
    set_param_(name, value);
  }

  /// On host: Sets a named real value for this aerosol process.
  /// @param [in] name The name of the parameter to set
  /// @param [in] value The parameter's value
  void set_param(const std::string& name, Real value) {
    // This method must be called on the host.
    set_param_(name, value);
  }

  /// On host: Sets a named string value for this aerosol process.
  /// @param [in] name The name of the parameter to set
  /// @param [in] value The parameter's value
  void set_param(const std::string& name, const std::string& value) {
    // This method must be called on the host.
    set_param_(name, value);
  }

  /// On host: returns a vector of strings containing the names of diagnostic
  /// variables required by this aerosol process in order to compute its
  /// tendencies. The Model that runs this aerosol process checks that these
  /// diagnostics variables are present before executing the process.
  std::vector<std::string> required_diagnostics() const {
    // This method must be called on the host.
    return required_diagnostics_();
  }

  /// On host: copies this aerosol process to the device, returning a
  /// pointer to the copy.
  AerosolProcess* copy_to_device() const { return copy_to_device_(); }

  /// On host: call this static method to delete a copy of the process
  /// that has been created on a device.
  /// @param [inout] device_process A pointer to a process created with
  ///                copy_to_device().
  KOKKOS_INLINE_FUNCTION
  static void delete_on_device(AerosolProcess* device_process) {
    Kokkos::parallel_for(
        "delete", 1,
        KOKKOS_LAMBDA(const int) { device_process->~AerosolProcess(); });
    Kokkos::kokkos_free<MemorySpace>((void*)device_process);
  }

 protected:
  /// On host: override this method to perform system-specific initialization
  /// for the aerosol process. By default, this does nothing.
  virtual void init_(const ModalAerosolConfig& config) {}

  /// On device: override this method to run the aerosol process using its
  /// specific parameterizations for a subclass.
  /// @param [in] team The Kokkos team used to run this process in a parallel
  ///                  dispatch. Not used by Fortran process implementations.
  /// @param [in] t The simulation time at which this process is being invoked
  ///               (in seconds).
  /// @param [in] dt The simulation time interval ("timestep size") over which
  ///                this process occurs.
  /// @param [in] prognostics The prognostic variables used by and affected by
  ///                         this process.
  /// @param [in] atmosphere The atmosphere state variables used by this
  ///                        process.
  /// @param [in] diagnostics The prognostic variables used by and affected by
  ///                         this process.
  /// @param [out] tendencies A container that stores time derivatives for
  ///                         prognostic variables evolved by this process.
  KOKKOS_FUNCTION
  virtual void run_(const TeamType& team, Real t, Real dt,
                    const Prognostics& prognostics,
                    const Atmosphere& atmosphere,
                    const Diagnostics& diagnostics,
                    const Tendencies& tendencies) const = 0;

  /// Override this method to return a vector of strings containing the names
  /// of diagnostic variables required by this aerosol process in order to
  /// compute its tendencies. By default, an aerosol process does not require
  /// any diagnostic variables.
  virtual std::vector<std::string> required_diagnostics_() const {
    return std::vector<std::string>();
  }

  /// Override these methods to allow a host model to set named parameters
  /// for this aerosol process.
  virtual void set_param_(const std::string& name, int value) {}
  virtual void set_param_(const std::string& name, bool value) {}
  virtual void set_param_(const std::string& name, Real value) {}
  virtual void set_param_(const std::string& name, const std::string& value) {}

  /// This gets overridden by the AerosolProcessOnDevice middleware class.
  virtual AerosolProcess* copy_to_device_() const = 0;

 private:
  const AerosolProcessType type_;
  // Use View as a struct to store a string and allows copy to device.
  // Since std::string can not be used, it was either this or a char *
  const Kokkos::View<int> name_;
  RegionOfValidity validity_region_;
};

/// @class DeviceAerosolProcess
/// This "middleware" class provides all its subclasses with the ability to copy
/// itself from the host to the device. All actual aerosol processes derive from
/// this type with their own type as the template parameter, in accordance with
/// the C++ "curiously recurring template pattern" (see
/// https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern).
template <typename Subclass>
class DeviceAerosolProcess : public AerosolProcess {
 public:
  /// Constructor, called by all DeviceAerosolProcess subclasses.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  DeviceAerosolProcess(AerosolProcessType type, const std::string& name)
      : AerosolProcess(type, name) {}

 protected:
  AerosolProcess* copy_to_device_() const override {
    const std::string debug_name = name();
    Subclass* process =
        static_cast<Subclass*>(Kokkos::kokkos_malloc<MemorySpace>(
            debug_name + "_malloc", sizeof(Subclass)));

    // Copy this object (including our virtual table) into the storage using
    // a lambda capture.
    const auto* this_process = dynamic_cast<const Subclass*>(this);
    Kokkos::parallel_for(
        debug_name + "_copy", 1,
        KOKKOS_LAMBDA(const int) { new (process) Subclass(*this_process); });
    return process;
  }
};

/// @class NullAerosolProcess
/// This AerosolProcess represents a null process in which no tendencies
/// are computed. It can be used for all prognostic processes that have been
/// disabled.
class NullAerosolProcess : public DeviceAerosolProcess<NullAerosolProcess> {
 public:
  /// Constructor: constructs a null aerosol process of the given type.
  /// @param [in] type The type of aerosol process.
  explicit NullAerosolProcess(AerosolProcessType type)
      : DeviceAerosolProcess<NullAerosolProcess>(
            type, "Null prognostic aerosol process") {}

  // Overrides
  KOKKOS_FUNCTION
  void run_(const TeamType& team, Real t, Real dt,
            const Prognostics& prognostics, const Atmosphere& atmosphere,
            const Diagnostics& diagnostics,
            const Tendencies& tendencies) const override {}
};

#if HAERO_FORTRAN

/// @class FAerosolProcess
/// This AerosolProcess makes it easier to implement an aerosol process in
/// Fortran by wrapping the run method in a Fortran call that creates the
/// proper Fortran proxies for the model, prognostics, diagnostics, and
/// tendencies.
class FAerosolProcess : public DeviceAerosolProcess<FAerosolProcess> {
 public:
  /// A pointer to a Fortran subroutine that initializes a process.
  /// Since the model is already available to the Fortran process via Haero's
  /// Fortran helper module, this type of function takes no arguments.
  typedef void (*InitProcessSubroutine)(void);

  /// A pointer to a Fortran subroutine that runs a process.
  typedef void (*RunProcessSubroutine)(Real t, Real dt, void* progs, void* atm,
                                       void* diags, void* tends);

  /// A pointer to a Fortran subroutine that finalizes a process,
  /// deallocating any resources allocated in the process's init function.
  typedef void (*FinalizeProcessSubroutine)(void);

  /// Pointers to Fortran subroutines that set parameters for a process.
  typedef void (*SetIntegerParamSubroutine)(const char* name, int value);
  typedef void (*SetLogicalParamSubroutine)(const char* name, bool value);
  typedef void (*SetRealParamSubroutine)(const char* name, Real value);

  /// Constructor.
  /// @param [in] type The type of aerosol process modeled by the subclass.
  /// @param [in] name A descriptive name that captures the aerosol process,
  ///                  its underlying parametrization, and its implementation.
  /// @param [in] init_process A pointer to an interoperable Fortran subroutine
  ///                          that initializes the Fortran-backed process.
  /// @param [in] run_process A pointer to an interoperable Fortran subroutine
  ///                         that runs the Fortran-backed process.
  /// @param [in] finalize_process A pointer to an interoperable Fortran
  /// subroutine
  ///                              that frees resources for the Fortran-backed
  ///                              process.
  /// @param [in] set_integer_param A pointer to an interoperable Fortran
  ///                               subroutine that sets an integer parameter
  ///                               for the Fortran-backed process.
  /// @param [in] set_logical_param A pointer to an interoperable Fortran
  ///                               subroutine that sets a logical parameter
  ///                               for the Fortran-backed process.
  /// @param [in] set_real_param A pointer to an interoperable Fortran
  ///                            subroutine that sets a real-valued parameter
  ///                            for the Fortran-backed process.
  FAerosolProcess(AerosolProcessType type, const std::string& name,
                  InitProcessSubroutine init_process,
                  RunProcessSubroutine run_process,
                  FinalizeProcessSubroutine finalize_process,
                  SetIntegerParamSubroutine set_integer_param,
                  SetLogicalParamSubroutine set_logical_param,
                  SetRealParamSubroutine set_real_param)
      : DeviceAerosolProcess<FAerosolProcess>(type, name),
        init_process_(init_process),
        run_process_(run_process),
        finalize_process_(finalize_process),
        set_integer_param_(set_integer_param),
        set_logical_param_(set_logical_param),
        set_real_param_(set_real_param),
        initialized_(false) {}

  /// Copy constructor.
  FAerosolProcess(const FAerosolProcess& pp)
      : DeviceAerosolProcess<FAerosolProcess>(pp),
        init_process_(pp.init_process_),
        run_process_(pp.run_process_),
        finalize_process_(pp.finalize_process_),
        set_integer_param_(pp.set_integer_param_),
        set_logical_param_(pp.set_logical_param_),
        set_real_param_(pp.set_real_param_),
        initialized_(pp.initialized_) {}

  /// Destructor.
  ~FAerosolProcess() {
    if (initialized_) {
      finalize_process_();
      initialized_ = false;
    }
  }

 protected:
  // Overrides.
  void init_(const ModalAerosolConfig& modal_aerosol_config) override {
    if (not initialized_) {
      init_process_();
      initialized_ = true;
    }
  }

  void run_(const TeamType& team, Real t, Real dt,
            const Prognostics& prognostics, const Atmosphere& atmosphere,
            const Diagnostics& diagnostics,
            const Tendencies& tendencies) const override {
    // Set tendencies to zero.
    auto nc_tendencies = const_cast<Tendencies&>(tendencies);
    nc_tendencies.scale(0.0);

    // Call the Fortran subroutine for the process.
    run_process_(t, dt, (void*)&prognostics, (void*)&atmosphere,
                 (void*)&diagnostics, (void*)&tendencies);
  }

  void set_param_(const std::string& name, int value) override {
    set_integer_param_(name.c_str(), value);
  }

  void set_param_(const std::string& name, bool value) override {
    set_logical_param_(name.c_str(), value);
  }

  void set_param_(const std::string& name, Real value) override {
    set_real_param_(name.c_str(), value);
  }

 private:
  // Pointers to Fortran subroutines.
  InitProcessSubroutine init_process_;
  RunProcessSubroutine run_process_;
  FinalizeProcessSubroutine finalize_process_;

  SetIntegerParamSubroutine set_integer_param_;
  SetLogicalParamSubroutine set_logical_param_;
  SetRealParamSubroutine set_real_param_;

  // Has the process been initialized?
  bool initialized_;
};

#endif  // HAERO_FORTRAN

}  // namespace haero

#endif
