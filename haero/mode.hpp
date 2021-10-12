#ifndef HAERO_MODE_HPP
#define HAERO_MODE_HPP

#include <cmath>
#include <string>
#include <vector>

#include "ekat/ekat_pack.hpp"
#include "haero/constants.hpp"
#include "haero/haero.hpp"
#include "haero/math.hpp"

namespace haero {

/// @struct Mode
/// This struct represents an aerosol particle mode and contains all associated
/// metadata. By definition, these metadata are immutable (constant in time).
/// The struct is not polymorphic, so don't derive any subclass from it.
///
/// This class represents the log-normal distribution that defines the mode via
/// the mean_std_dev member variable.  The other parameter necessary to define
/// the log-normal function is a variable (a function of mass- and number-
/// mixing ratios) and is not included in this class.
///
/// The member variables min_diameter and max_diameter do not define the bounds
/// of the log-normal distribution (which, matematically, are 0 and positive
/// infinity).  Rather, these min/max values are used to trigger a mass and
/// number redistribution elsewhere in the code; they signify the bounds beyond
/// which particles are considered to better belong in a different mode.
///
/// Variable nom_diameter is the nominal geometeric mean diameter [m]
/// of particles in a mode
///
/// Crystalization and deliquesence refer to the non-cloud water uptake process,
/// by which liquid water condenses into aerosol droplets.  They are relative
/// humidity values.  When the environmental relative humidity lies below the
/// cyrstalization point, water uptake does not occur.  When it lies between the
/// crystallization and deliquesence point, water uptake does occur, but not at
/// its maximum rate.   When the environmental relative humidty exceeds the
/// deliquescence_pt, particles achieve their maximum amount of liquid water.
///
struct Mode final {
  static const int NAME_LEN = 128;

 public:
  // Default constructor needed to resize Kokkos Views on device before deep
  // copy.
  KOKKOS_INLINE_FUNCTION
  Mode()
      : min_diameter(0),
        nom_diameter(0),
        max_diameter(0),
        mean_std_dev(1),
        crystallization_pt(0),
        deliquescence_pt(0) {
    name_view[0] = '\0';
  }
  /// Creates a new aerosol particle mode.
  /// @param [in] name A unique name for this mode.
  /// @param [in] min_diam The minimum diameter for particles that belong
  ///                      to this mode [m].
  /// @param [in] nom_diam The nominal diameter for particles that belong
  ///                      to this mode [m].
  /// @param [in] max_diam The maximum diameter for particles that belong
  ///                      to this mode [m].
  /// @param [in] sigma    The geometric standard deviation for this mode.
  /// @param [in] crystal_pt The crystallization point of the mode
  /// @param [in] deliq_pt The deliquescence point of the mode
  Mode(const std::string &name, Real min_diam, Real nom_diam, Real max_diam,
       Real sigma, Real crystal_pt, Real deliq_pt)
      : min_diameter(min_diam),
        nom_diameter(nom_diam),
        max_diameter(max_diam),
        mean_std_dev(sigma),
        crystallization_pt(crystal_pt),
        deliquescence_pt(deliq_pt) {
    EKAT_ASSERT(max_diam > min_diam);
    EKAT_ASSERT(nom_diam > min_diam);
    EKAT_ASSERT(max_diam > nom_diam);
    EKAT_ASSERT(deliq_pt > crystal_pt);
    EKAT_ASSERT(sigma >= 1);
    EKAT_ASSERT(name.size() < NAME_LEN);
    strncpy(name_view, name.c_str(), NAME_LEN);
  }

  KOKKOS_INLINE_FUNCTION
  Mode(const Mode &m)
      : min_diameter(m.min_diameter),
        nom_diameter(m.nom_diameter),
        max_diameter(m.max_diameter),
        mean_std_dev(m.mean_std_dev),
        crystallization_pt(m.crystallization_pt),
        deliquescence_pt(m.deliquescence_pt) {
    for (int i = 0; i < NAME_LEN; ++i) name_view[i] = m.name_view[i];
  }

  KOKKOS_INLINE_FUNCTION
  Mode &operator=(const Mode &m) {
    min_diameter = m.min_diameter;
    nom_diameter = m.nom_diameter;
    max_diameter = m.max_diameter;
    mean_std_dev = m.mean_std_dev;
    crystallization_pt = m.crystallization_pt;
    deliquescence_pt = m.deliquescence_pt;
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

  /// The nominal diameter for particles that belong to this mode.
  Real nom_diameter;

  /// The maximum diameter for particles that belong to this mode.
  Real max_diameter;

  /// The geometric mean standard deviation for this mode.
  Real mean_std_dev;

  /// The crystallization point (rel. humidity) for this mode.
  Real crystallization_pt;

  /// The deliquescence point (rel. humidity) for this mode.
  Real deliquescence_pt;

  /** @brief This function returns the modal geometric mean particle diameter,
  given the mode's mean volume (~ to 3rd log-normal moment) and the modal
  standard deviation.

  @param mode_mean_particle_volume mean particle volume for mode [m^3 per
  particle]
  @return modal mean particle diameter [m per particle]
*/
  template <typename T>
  KOKKOS_INLINE_FUNCTION T
  mean_particle_diameter_from_volume(const T mode_mean_particle_volume) const {
    const Real pio6 = Constants::pi_sixth;
    return cbrt(mode_mean_particle_volume / pio6) *
           exp(-1.5 * square(log(mean_std_dev)));
  }

  /** @brief This function is the inverse of
    modal_mean_particle_diameter_from_volume; given the modal mean geometric
    diamaeter, it returns the corresponding volume.

    @param [in] geom_diam geometric mean diameter [m per particle]
    @return mean volume [m^3 per particle]
  */
  template <typename T>
  KOKKOS_INLINE_FUNCTION T
  mean_particle_volume_from_diameter(const T geom_diam) const {
    const Real pio6 = Constants::pi_sixth;
    return cube(geom_diam) * exp(4.5 * square(log(mean_std_dev))) * pio6;
  }

  /** @brief This function returns the minimum volume to number ratio,
      which is computed using the maximum diameter(units:meters) and
      modal standard deviation

      @return modal minimum volume to number ratio [m^-3] FIXME: Check the units
     again
  */
  template <typename T>
  KOKKOS_INLINE_FUNCTION T min_vol_to_num_ratio() {
    return 1 / mean_particle_volume_from_diameter(max_diameter);
  }

  /** @brief This function returns the maximum volume to number ratio,
      which is computed using the minimum diameter(units:meters) and
      modal standard deviation

      @return modal maximum volume to number ratio [m^-3] FIXME: Check the units
     again
  */
  template <typename T>
  KOKKOS_INLINE_FUNCTION T max_vol_to_num_ratio() {
    return 1 / mean_particle_volume_from_diameter(min_diameter);
  }

 private:
  char name_view[NAME_LEN];
};

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
  const std::vector<Real> mode_nom_diam = {1.1e-07, 2.6e-08, 2e-06, 5e-08};
  const std::vector<Real> mode_max_diam = {4.4e-7, 5.2e-8, 4e-6, 1e-7};
  const std::vector<Real> mode_std_dev = {1.8, 1.6, 1.8, 1.6};
  std::vector<Mode> result(4);
  for (int i = 0; i < 4; ++i) {
    result[i] = Mode(mode_names[i], mode_min_diam[i], mode_nom_diam[i],
                     mode_max_diam[i], mode_std_dev[i], rh_crystal, rh_deliq);
  }
  return result;
}

}  // namespace haero

#endif
