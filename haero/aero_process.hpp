#ifndef HAERO_AERO_PROCESS_HPP
#define HAERO_AERO_PROCESS_HPP
#include <memory>

#include "haero/atmosphere.hpp"

namespace haero {

/// @class AeroProcess
/// This type defines the interface for a specific process in the aerosol
/// lifecycle, backed by a specific implementation, the structure of which is
/// defined by a specific "aerosol configuration".
template <typename AeroConfig, typename AeroImpl>
class AeroProcess final {
 public:

  // Types derived from template parameters.
  using Config      = AeroConfig;
  using Prognostics = typename AeroConfig::Prognostics;
  using Diagnostics = typename AeroConfig::Diagnostics;
  using Tendencies  = typename AeroConfig::Prognostics; // look!!
  using Impl        = AeroImpl;

  /// Constructs an instance of an aerosol process with the given name,
  /// associated with the given aerosol configuration.
  /// @param [in] name A descriptive name that uniquely identifies a specific
  ///                  implementation of a process in the aerosol life cycle.
  /// @param [in] config An instance of a type that defines any metadata needed
  ///                    by this process's implementation.
  AeroProcess(const std::string& name, const Config& config)
      : name_(name), config_(config), impl_() {}

  /// Destructor.
  KOKKOS_INLINE_FUNCTION ~AeroProcess() {}

  /// Default constructor is disabled.
  AeroProcess() = delete;

  // Deep copying is forbidden.
  AeroProcess(const AeroProcess&) = delete;
  AeroProcess& operator=(const AeroProcess&) = delete;

  //------------------------------------------------------------------------
  //                          Accessors (host only)
  //------------------------------------------------------------------------

  /// Returns the name of this process.
  std::string name() const { return name_.label(); }

  /// Returns the metadata associated with this process.
  const Config& config() const {
    return config_;
  }

  //------------------------------------------------------------------------
  //                            Public Interface
  //------------------------------------------------------------------------

  /// On host: performs any system-specific process initialization.
  /// Initialization is typically performed after process parameters are set
  /// with set_param().
  /// @param [in] config The aerosol configuration describing the aerosol
  ///                    system to which this process belongs.
  void init() {
    // Pass the configuration data to the implementation to initialize it.
    impl_.init(config_);
  }

  /// On host or device: Validates input aerosol and atmosphere data, returning
  /// true if all data is physically consistent (whatever that means), and false
  /// if not.
  /// @param [in] atmosphere Atmosphere state variables with which to validate.
  /// @param [in] prognostics A collection of aerosol prognostic variables to be
  ///                         validated.
  KOKKOS_INLINE_FUNCTION
  bool validate(const Atmosphere& atmosphere,
                const Prognostics& prognostics) const {
    return impl_.validate(config_, atmosphere, prognostics);
  }

  /// On host or device: runs the aerosol process at a given time with the given
  /// data.
  /// @param [in]    team The Kokkos team used to run this process in a parallel
  ///                     dispatch.
  /// @param [in]    t The simulation time at which this process is being
  ///                  invoked (in seconds).
  /// @param [in]    dt The simulation time interval ("timestep size") over
  ///                   which this process occurs.
  /// @param [in]    atmosphere The atmosphere state variables used by this
  ///                           process.
  /// @param [in]    prognostics An array containing aerosol tracer data to be
  ///                            evolved.
  /// @param [inout] diagnostics An array that can store aerosol diagnostic
  ///                            data computed or updated by this process.
  /// @param [out]   tendencies An array analogous to prognostics that
  ///                           stores computed tendencies.
  KOKKOS_INLINE_FUNCTION
  void compute_tendencies(const TeamType& team, Real t, Real dt,
                          const Atmosphere& atmosphere,
                          const Prognostics& prognostics,
                          Diagnostics& diagnostics,
                          Tendencies& tendencies) const {
    impl_.compute_tendencies(config_, team, t, dt, atmosphere, prognostics,
                             diagnostics, tendencies);
  }

 private:
  const std::string name_;
  Config            config_;
  Impl              impl_;
};

}  // namespace haero

#endif
