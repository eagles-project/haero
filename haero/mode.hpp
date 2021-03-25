#ifndef HAERO_MODE_HPP
#define HAERO_MODE_HPP

#include "haero/haero.hpp"

#include <string>
#include <vector>

namespace haero {

/// @struct Mode
/// This struct represents an aerosol particle mode and contains all associated
/// metadata. It is not polymorphic, so don't derive any subclass from it.
struct Mode final {
  public:

  // Default constructor needed to resize Kokkos Views on device before deep copy.
  KOKKOS_INLINE_FUNCTION
  Mode() :
    min_diameter(0),
    max_diameter(0),
    mean_std_dev(0)
  { 
    name_view[0]='\0';
  }
  /// Creates a new aerosol particle mode.
  /// @param [in] name A unique name for this mode.
  /// @param [in] min_diameter The minimum diameter for particles that belong
  ///                          to this mode.
  /// @param [in] max_diameter The maximum diameter for particles that belong
  ///                          to this mode.
  /// @param [in] mean_std_dev The geometric standard deviation for this mode.
  /// @param [in] rh_crystal TODO
  /// @param [in] rh_deliques TODO
  Mode(const std::string& name,
       Real min_diameter,
       Real max_diameter,
       Real mean_std_dev):
    min_diameter(min_diameter),
    max_diameter(max_diameter),
    mean_std_dev(mean_std_dev)
  {
    EKAT_ASSERT(name.size() < 100);
    strncpy(name_view, name.c_str(), 100);
  }

  KOKKOS_INLINE_FUNCTION
  Mode(const Mode &m):
    min_diameter(m.min_diameter),
    max_diameter(m.max_diameter),
    mean_std_dev(m.mean_std_dev)
  {
    for (int i=0; i<100; ++i) 
      name_view[i] = m.name_view[i];
  }

  /// Constructor can be called on device.
  KOKKOS_INLINE_FUNCTION
  ~Mode() {}

  /// A unique name for this mode.
  std::string name() const { return std::string(name_view); }

  /// The minimum diameter for particles that belong to this mode.
  Real min_diameter;

  /// The maximum diameter for particles that belong to this mode.
  Real max_diameter;

  /// The geometric mean standard deviation for this mode.
  Real mean_std_dev;

private:
  char name_view[100];
};

/// This factory function constructs a set of modes corresponding to the
/// legacy MAM4 model. Four modes are included:
/// 1. Aitken,         with D_min = 0.0087, D_max = 0.052, sigma = 1.6,
/// 2. Accumulation,   with D_min = 0.0535, D_max =  0.44, sigma = 1.8,
/// 3. Coarse,         with D_min = 1.0,    D_max =  4.0,  sigma = 1.8,
/// 4. Primary carbon, with D_min = 0.01,   D_max =  0.1,  sigma = 1.6.
inline std::vector<Mode> create_mam4_modes() {
  return std::vector<Mode>({
    Mode("aitken", 0.0087, 0.052, 1.6),
    Mode("accumulation", 0.0535, 0.44, 1.8),
    Mode("coarse", 1.0, 0.01, 1.8),
    Mode("primary_carbon", 0.01, 0.1, 1.6)
  });
}

}

#endif
