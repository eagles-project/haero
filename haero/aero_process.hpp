#ifndef HAERO_AEROSOL_PROCESS_HPP
#define HAERO_AEROSOL_PROCESS_HPP
#include <memory>

#include "haero/atmosphere.hpp"

namespace haero {

/// @class AerosolProcess
/// This type defines the interface for a specific process in the aerosol
/// lifecycle, backed by a specific implementation, the structure of which is
/// defined by a specific "aerosol configuration".
template <typename AerosolConfig, typename AerosolImpl>
class AerosolProcess final {
 public:

  using ConfigType = AerosolConfig;
  using ImplType = AerosolImpl;

  /// Constructs an instance of an aerosol process with the given name,
  /// associated with the given aerosol configuration.
  /// @param [in] name A descriptive name that uniquely identifies a specific
  ///                  implementation of a process in the aerosol life cycle.
  /// @param [in] config An instance of a type that defines any metadata needed
  ///                    by this process's implementation.
  AerosolProcess(const std::string& name, const ConfigType &config)
      : name_(name), config_(config), impl_() {}

  /// Destructor.
  KOKKOS_INLINE_FUNCTION ~AerosolProcess() {}

  /// Default constructor is disabled.
  AerosolProcess() = delete;

  // Deep copying is forbidden.
  AerosolProcess(const AerosolProcess& pp) = delete;
  AerosolProcess& operator=(const AerosolProcess&) = delete;

  //------------------------------------------------------------------------
  //                                Accessors
  //------------------------------------------------------------------------

  /// Returns the name of this process.
  std::string name() const { return name_.label(); }

  /// Returns the metadata associated with this process.
  const ConfigType& config() const {
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
  void init(const ConfigType& config) {
    // Pass the configuration data to the implementation.
    impl_.configure(config_);
  }

  /// On host or device: Validates input aerosol and atmosphere data, returning
  /// true if all data is physically consistent (whatever that means), and false
  /// if not.
  /// @param [in] atmosphere Atmosphere state variables with which to validate.
  /// @param [in] aero_tracers An array containing aerosol tracer data to be
  ///                          validated.
  KOKKOS_INLINE_FUNCTION
  bool validate(const Atmosphere& atmosphere,
                const TracersView& aero_tracers) const {
    return impl_.validate(atmosphere, aerosol_tracers);
  }

  /// On host or device: runs the aerosol process at a given time with the given
  /// data.
  /// @param [in] team The Kokkos team used to run this process in a parallel
  ///                  dispatch.
  /// @param [in] t The simulation time at which this process is being invoked
  ///               (in seconds).
  /// @param [in] dt The simulation time interval ("timestep size") over which
  ///                this process occurs.
  /// @param [in] atmosphere The atmosphere state variables used by this
  ///                        process.
  /// @param [in] aero_tracers An array containing aerosol tracer data to be
  ///                          evolved.
  /// @param [in] aero_diagnostics An array that can store aerosol diagnostic
  ///                               data computed or updated by this process.
  /// @param [out] aero_tendencies An array analogous to aero_tracers that
  ///                              stores computed tendencies for those tracers.
  KOKKOS_FUNCTION
  void compute_tendencies(const TeamType& team, Real t, Real dt,
                          const Atmosphere& atmosphere,
                          const TracersView& aero_tracers,
                          const DiagnosticsView& aero_diagnostics,
                          const TracersView& aero_tendencies) const {
    // This method must be called on the device.
    impl_.compute_tendencies(team, t, dt, atmosphere, aero_tracers,
                             aero_diagnostics, aero_tendencies);
  }

 private:
  const std::string name_;
  ImplType impl_;
};

}  // namespace haero

#endif
