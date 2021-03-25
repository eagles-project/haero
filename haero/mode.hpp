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
    mean_std_dev(0),
    name_view(),
    crystallization_pt(0),
    deliquesence_pt(0) {}
  /// Creates a new aerosol particle mode.
  /// @param [in] name A unique name for this mode.
  /// @param [in] min_diameter The minimum diameter for particles that belong
  ///                          to this mode [m].
  /// @param [in] max_diameter The maximum diameter for particles that belong
  ///                          to this mode [m].
  /// @param [in] mean_std_dev The geometric standard deviation for this mode.
  /// @param [in] crystal_pt The crystallization point of the mode
  /// @param [in] deliq_pt The deliquescence point of the mode
  Mode(const std::string& name,
       Real min_diameter,
       Real max_diameter,
       Real mean_std_dev,
       Real crystal_pt,
       Real deliq_pt):
    name_view(name)
    min_diameter(min_diameter),
    max_diameter(max_diameter),
    mean_std_dev(mean_std_dev),
    deliquesence_pt(deliq_pt),
    crystallization_pt(crystal_pt)
     {}

  KOKKOS_INLINE_FUNCTION
  Mode(const Mode &m):
    name_view(m.name_view),
    min_diameter(m.min_diameter),
    max_diameter(m.max_diameter),
    mean_std_dev(m.mean_std_dev),
    deliquesence_pt(m.deliquesence_pt),
    crystallization_pt(m.crystallization_pt)
    {}

  /// Constructor can be called on device.
  KOKKOS_INLINE_FUNCTION
  ~Mode() {}

  /// A unique name for this mode.
  std::string name() const { return name_view.label(); }

  /// The minimum diameter for particles that belong to this mode.
  Real min_diameter;

  /// The maximum diameter for particles that belong to this mode.
  Real max_diameter;

  /// The geometric mean standard deviation for this mode.
  Real mean_std_dev;

  /// The deliquescence point (rel. humidity) for this mode.
  Real deliquesence_pt;

  /// The crystallization point (rel. humidity) for this mode.
  Real crystallization_pt;

private:
  Kokkos::View<int>  name_view;
};


inline std::vector<Mode> create_es3mv1_modes() {
  /// Legacy MAM4 used the same constant crystallization and deliquescence values for all
  /// modes & species.  See links for additional discussion:
  /// https://eagles-project.atlassian.net/wiki/spaces/Computation/pages/1125515265/Aerosol+species+and+mode+data
  /// https://eagles-project.atlassian.net/wiki/spaces/Computation/pages/354877515/Module+verifications
  /// These data are found on Anvil in  /lcrc/group/acme/ccsm-data/inputdata/atm/cam/physprops/
  static constexpr Real rh_crystal = 0.35;
  static constexpr Real rh_deliq = 0.8;
  const std::vector<std::string> mode_names = {"aitken", "accumulation", "coarse", "primary_carbon"};
  const std::vector<Real> mode_mean_diam =    {1.1e-7,    2.6e-8,         2e-6,     5e-8};
  const std::vector<Real> mode_min_diam =     {5.35e-8,   8.7e-9,         1e-6,     1e-8};
  const std::vector<Real> mode_max_diam =     {4.4e-7,    5.2e-8,         4e-6,     1e-7};

  std::vector<Mode> result(4);
  for (int i=0; i<4; ++i) {
    result[i] = Mode(mode_names[i], mode_min_diam[i], mode_max_diam[i], mode_mean_diam[i], rh_crystal, rh_deliq);
  }
  return result;
}

}

#endif
