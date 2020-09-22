#ifndef HAERO_MODE_HPP
#define HAERO_MODE_HPP

#include <string>
#include "haero/haero.hpp"

namespace haero {

/// This struct represents an aerosol particle mode and contains all associated
/// metadata. It is not polymorphic, so don't derive any subclass from it.
struct Mode final {
  public:

  /// Creates a new aerosol particle mode.
  /// @param [in] name A unique name for this mode.
  /// @param [in] min_diameter The minimum diameter for particles that belong
  ///                          to this mode.
  /// @param [in] max_diameter The maximum diameter for particles that belong
  ///                          to this mode.
  /// @param [in] mean_std_dev The geometric standard deviation for this mode.
  Mode(const std::string& name,
       Real min_diameter,
       Real max_diameter,
       Real mean_std_dev):
    name(name), min_diameter(min_diameter), max_diameter(max_diameter),
    mean_std_dev(mean_std_dev) {}

  /// A unique name for this mode.
  std::string name;

  /// The minimum diameter for particles that belong to this mode.
  Real min_diameter;

  /// The maximum diameter for particles that belong to this mode.
  Real max_diameter;

  /// The geometric mean standard deviation for this mode.
  Real mean_std_dev;
};

}

#endif
