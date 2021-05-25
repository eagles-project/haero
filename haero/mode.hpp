#ifndef HAERO_MODE_HPP
#define HAERO_MODE_HPP

#include <cmath>
#include <string>
#include <vector>

#include "ekat/ekat_pack.hpp"
#include "haero/haero.hpp"
#include "haero/math_helpers.hpp"
#include "haero/physical_constants.hpp"

namespace haero {

/// @struct Mode
/// This struct represents an aerosol particle mode and contains all associated
/// metadata. It is not polymorphic, so don't derive any subclass from it.
struct Mode final {
  static const int NAME_LEN = 128;

 public:
  // Default constructor needed to resize Kokkos Views on device before deep
  // copy.
  KOKKOS_INLINE_FUNCTION
  Mode()
      : min_diameter(0),
        max_diameter(0),
        mean_std_dev(1),
        deliquesence_pt(0),
        crystallization_pt(0) {
    name_view[0] = '\0';
  }
  /// Creates a new aerosol particle mode.
  /// @param [in] name A unique name for this mode.
  /// @param [in] min_diameter The minimum diameter for particles that belong
  ///                          to this mode [m].
  /// @param [in] max_diameter The maximum diameter for particles that belong
  ///                          to this mode [m].
  /// @param [in] mean_std_dev The geometric standard deviation for this mode.
  /// @param [in] crystal_pt The crystallization point of the mode
  /// @param [in] deliq_pt The deliquescence point of the mode
  Mode(const std::string &name, Real min_diam, Real max_diam, Real sigma,
       Real deliq_pt, Real crystal_pt)
      : min_diameter(min_diam),
        max_diameter(max_diam),
        mean_std_dev(sigma),
        deliquesence_pt(deliq_pt),
        crystallization_pt(crystal_pt) {
    EKAT_ASSERT(name.size() < NAME_LEN);
    strncpy(name_view, name.c_str(), NAME_LEN);
  }

  KOKKOS_INLINE_FUNCTION
  Mode(const Mode &m)
      : min_diameter(m.min_diameter),
        max_diameter(m.max_diameter),
        mean_std_dev(m.mean_std_dev),
        deliquesence_pt(m.deliquesence_pt),
        crystallization_pt(m.crystallization_pt) {
    for (int i = 0; i < NAME_LEN; ++i) name_view[i] = m.name_view[i];
  }

  KOKKOS_INLINE_FUNCTION
  Mode &operator=(const Mode &m) {
    min_diameter = m.min_diameter;
    max_diameter = m.max_diameter;
    mean_std_dev = m.mean_std_dev;
    deliquesence_pt = m.deliquesence_pt;
    crystallization_pt = m.crystallization_pt;
    for (int i = 0; i < NAME_LEN; ++i) name_view[i] = m.name_view[i];
    return *this;
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

  /// The deliquescence point (rel. humidity) for this mode.
  Real deliquesence_pt;

  /// The crystallization point (rel. humidity) for this mode.
  Real crystallization_pt;

  /** @brief This function returns the minimum volume to number ratio,
      which is computed using the maximum diameter(units:meters) and
      modal standard deviation

      @return modal minimum volume to number ratio [m^-3] FIXME: Check the units
     again
  */
  template <typename T>
  KOKKOS_INLINE_FUNCTION T min_vol_to_num_ratio() {
    return 1 / (constants::pi_sixth * (cube(max_diameter)) *
                exp(4.5 * square(log(mean_std_dev))));
  }

  /** @brief This function returns the maximum volume to number ratio,
      which is computed using the minimum diameter(units:meters) and
      modal standard deviation

      @return modal maximum volume to number ratio [m^-3] FIXME: Check the units
     again
  */
  template <typename T>
  KOKKOS_INLINE_FUNCTION T max_vol_to_num_ratio() {
    return 1 / (constants::pi_sixth * (cube(min_diameter)) *
                exp(4.5 * square(log(mean_std_dev))));
  }

 private:
  char name_view[NAME_LEN];
};

/** @brief This function returns the modal geometric mean particle diameter,
  given the mode's mean volume (3rd log-normal moment) and the modal standard
  deviation.

  @param mode_mean_particle_volume mean particle volume for mode [m^3]
  @param log_sigma natural log of the mode's geometric mean std. dev.
  @return modal mean particle diameter (~ 1st log-normal moment) [m]
*/
template <typename T>
KOKKOS_INLINE_FUNCTION T modal_mean_particle_diameter(
    const T mode_mean_particle_volume, const Real log_sigma) {
  return cbrt(constants::pi_sixth * mode_mean_particle_volume) *
         exp(-1.5 * square(log_sigma));
}

inline std::vector<Mode> create_mam4_modes() {
  /// Legacy MAM4 used the same constant crystallization and deliquescence
  /// values for all modes & species.  See links for additional discussion:
  /// https://eagles-project.atlassian.net/wiki/spaces/Computation/pages/1125515265/Aerosol+species+and+mode+data
  /// https://eagles-project.atlassian.net/wiki/spaces/Computation/pages/354877515/Module+verifications
  /// These data are found on Anvil in
  /// /lcrc/group/acme/ccsm-data/inputdata/atm/cam/physprops/
  static constexpr Real rh_crystal = 0.35;
  static constexpr Real rh_deliq = 0.8;
  const std::vector<std::string> mode_names = {"accumulation", "aitken",
                                               "coarse", "primary_carbon"};
  const std::vector<Real> mode_min_diam = {5.35e-8, 8.7e-9, 1e-6, 1e-8};
  const std::vector<Real> mode_max_diam = {4.4e-7, 5.2e-8, 4e-6, 1e-7};
  const std::vector<Real> mode_std_dev = {1.8, 1.6, 1.8, 1.6};
  std::vector<Mode> result(4);
  for (int i = 0; i < 4; ++i) {
    result[i] = Mode(mode_names[i], mode_min_diam[i], mode_max_diam[i],
                     mode_std_dev[i], rh_crystal, rh_deliq);
  }
  return result;
}

}  // namespace haero

#endif
