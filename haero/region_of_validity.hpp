#ifndef HAERO_REGION_OF_VALIDITY_HPP
#define HAERO_REGION_OF_VALIDITY_HPP

#include <algorithm>

#include "haero/atmosphere.hpp"
#include "haero/haero.hpp"
#include "haero/prognostics.hpp"

namespace haero {

/// @class RegionOfValidity
/// This type expresses the "region of validity" for a given aerosol process
/// in terms of bounds on atmospheric and aerosol prognostic states.
class RegionOfValidity final {
 public:
  using Bounds = Kokkos::pair<Real, Real>;

  /// Constructor
  RegionOfValidity()
      : temp_bounds({0, 500}),
        rel_hum_bounds({0, 1}),
        int_aero_indices_("Interstitial aerosol indices", 0),
        cld_aero_indices_("Cloudborne aerosol indices", 0),
        int_aero_bounds_("Interstitial aerosol bounds", 0),
        cld_aero_bounds_("Cloudborne aerosol bounds", 0),
        int_n_indices_("Interstitial number conc indices", 0),
        cld_n_indices_("Cloudborne number conc indices", 0),
        int_n_bounds_("Interstitial number conc bounds", 0),
        cld_n_bounds_("Cloudborne number conc bounds", 0),
        gas_indices_("Gas indices", 0),
        gas_bounds_("Gas bounds", 0) {}

  /// Destructor.
  KOKKOS_FUNCTION
  ~RegionOfValidity() {}

  /// Minimum and maximum bounds on atmospheric temperature [K]
  Bounds temp_bounds;
  /// Minimum and maximum bounds on relative humidity [-]
  Bounds rel_hum_bounds;

  /// Adds a set of bounds for the interstitial aerosol species with the given
  /// population index. Callable from host only.
  /// @param [in] pop_index The index of the aerosol species within a mode as
  ///                       indexed by a population index in a prognostics
  ///                       container
  /// @param [in] min_mmr The minimum valid mass mixing ratio [kg aero/kg air]
  /// @param [in] max_mmr The maximum valid mass mixing ratio [kg aero/kg air]
  void set_interstitial_aerosol_mmr_bounds(int pop_index, Real min_mmr,
                                           Real max_mmr) {
    insert_bounds_at_index_(pop_index, min_mmr, max_mmr, int_aero_indices_,
                            int_aero_bounds_);
  }

  /// Adds a set of bounds for the cloudborne aerosol species with the given
  /// population index. Callable from host only.
  /// @param [in] pop_index The index of the aerosol species within a mode as
  ///                       indexed by a population index in a prognostics
  ///                       container
  /// @param [in] min_mmr The minimum valid mass mixing ratio [kg aero/kg air]
  /// @param [in] max_mmr The maximum valid mass mixing ratio [kg aero/kg air]
  void set_cloud_aerosol_mmr_bounds(int pop_index, Real min_mmr, Real max_mmr) {
    insert_bounds_at_index_(pop_index, min_mmr, max_mmr, int_aero_indices_,
                            int_aero_bounds_);
  }

  /// Adds a set of bounds for the interstitial aerosol number densities for the
  /// mode with the given index. Callable from host only.
  /// @param [in] mode_index The index of the mode of interest
  /// @param [in] min_n The minimum valid number concentration [m-3]
  /// @param [in] max_n The maximum valid number concentration [m-3]
  void set_interstitial_aerosol_num_conc_bounds(int mode_index, Real min_n,
                                                Real max_n) {
    insert_bounds_at_index_(mode_index, min_n, max_n, int_aero_indices_,
                            int_aero_bounds_);
  }

  /// Adds a set of bounds for the cloudborne aerosol number densities for the
  /// mode with the given index. Callable from host only.
  /// @param [in] mode_index The index of the mode of interest
  /// @param [in] min_n The minimum valid number concentration [m-3]
  /// @param [in] max_n The maximum valid number concentration [m-3]
  void set_cloud_aerosol_num_conc_bounds(int mode_index, Real min_n,
                                         Real max_n) {
    insert_bounds_at_index_(mode_index, min_n, max_n, cld_aero_indices_,
                            cld_aero_bounds_);
  }

  /// Adds a set of bounds for the gas species with the given index. Callable
  /// from host only.
  /// @param [in] gas_index The index of the gas as it appears in the set of
  ///                       gases in a prognostics container
  /// @param [in] min_mmr The minimum valid mass mixing ratio [kg gas/kg air]
  /// @param [in] max_mmr The maximum valid mass mixing ratio [kg gas/kg air]
  void set_gas_mmr_bounds(int gas_index, Real min_mmr, Real max_mmr) {
    insert_bounds_at_index_(gas_index, min_mmr, max_mmr, gas_indices_,
                            gas_bounds_);
  }

  /// Minimum and maximum bounds on the mass mixing ratio for the interstitial
  /// aerosol corresponding to the given population index. Callable from host
  /// only.
  const Bounds interstitial_aerosol_mmr_bounds(int pop_index) const {
    return get_bounds_(int_aero_indices_, int_aero_bounds_, pop_index, 0.0,
                       1.0);
  }

  /// Minimum and maximum bounds on the mass mixing ratio for the cloudborne
  /// aerosol corresponding to the given population index. Callable from host
  /// only.
  const Bounds cloud_aerosol_mmr_bounds(int pop_index) const {
    return get_bounds_(cld_aero_indices_, cld_aero_bounds_, pop_index, 0.0,
                       1.0);
  }

  /// Minimum and maximum bounds on the number concentration for the
  /// interstitial aerosols in the mode with the given index. Callable from host
  /// only.
  const Bounds interstitial_aerosol_num_conc_bounds(int mode_index) const {
    return get_bounds_(int_n_indices_, int_n_bounds_, mode_index, 0.0, 1e20);
  }

  /// Minimum and maximum bounds on the number concentration for the cloudborne
  /// aerosols in the mode with the given index. Callable from host only.
  KOKKOS_INLINE_FUNCTION
  const Bounds cloud_aerosol_num_conc_bounds(int mode_index) const {
    return get_bounds_(cld_n_indices_, cld_n_bounds_, mode_index, 0.0, 1e20);
  }

  /// Minimum and maximum bounds on the mass mixing ratio for the gas species
  /// corresponding to the given index. Callable from host only.
  KOKKOS_INLINE_FUNCTION
  const Bounds gas_mmr_bounds(int gas_index) const {
    return get_bounds_(gas_indices_, gas_bounds_, gas_index, 0.0, 1.0);
  }

  /// Returns true if all state data in the given atmosphere container falls
  /// within this region of validity, false otherwise.
  KOKKOS_INLINE_FUNCTION
  bool contains(const Atmosphere& atmosphere) const {
    int violations = 0;
    int num_levels = atmosphere.temperature.extent(0);
    Kokkos::parallel_reduce(
        "RegionOfValidity::contains(atm)", num_levels,
        KOKKOS_LAMBDA(const int k, int& violation) {
          const auto& T = atmosphere.temperature(k);
          auto invalid_T = haero::MaskType((T < temp_bounds.first) or
                                           (T > temp_bounds.second));
          const auto& RH = atmosphere.relative_humidity(k);
          auto invalid_RH = haero::MaskType((RH < rel_hum_bounds.first) or
                                            (RH > rel_hum_bounds.second));
          violation += (invalid_T.any() or invalid_RH.any());
        },
        violations);
    return (violations == 0);
  }

  /// Returns true if all state data in the given prognostics container falls
  /// within this region of validity, false otherwise.
  KOKKOS_INLINE_FUNCTION
  bool contains(const Prognostics& prognostics) const {
    return (mmrs_are_valid_(prognostics.interstitial_aerosols,
                            int_aero_indices_, int_aero_bounds_) and
            mmrs_are_valid_(prognostics.cloud_aerosols, cld_aero_indices_,
                            cld_aero_bounds_) and
            num_concs_are_valid_(prognostics.interstitial_num_concs,
                                 int_n_indices_, int_n_bounds_) and
            num_concs_are_valid_(prognostics.cloudborne_num_concs,
                                 cld_n_indices_, cld_n_bounds_) and
            mmrs_are_valid_(prognostics.gases, gas_indices_, gas_bounds_));
  }

 private:
  using IndexArray = kokkos_device_type::view_1d<int>;
  using BoundsArray = kokkos_device_type::view_1d<Bounds>;

  // This helper inserts a new entry into an IndexArray/BoundsArray pair at
  // the appropriate (sorted) position.
  void insert_bounds_at_index_(int index, Real min, Real max,
                               IndexArray& indices, BoundsArray& bounds) {
    EKAT_ASSERT(index >= 0);
    EKAT_ASSERT(min < max);

    int* begin = indices.data();
    int* end = begin + indices.extent(0);
    int* iter = std::lower_bound(begin, end, index);
    int pos = (iter == end) ? 0 : *iter;
    if ((iter == end) || (*iter != index)) {  // index must be inserted
      Kokkos::resize(indices, indices.extent(0) + 1);
      Kokkos::resize(bounds, bounds.extent(0) + 1);
      for (int p = indices.extent(0) - 1; p > pos; --p) {
        indices(p) = indices(p - 1);
        bounds(p) = bounds(p - 1);
      }
      indices(pos) = index;
    }
    bounds(pos).first = min;
    bounds(pos).second = max;
  }

  // This helper returns the bounds found in the given array at the given
  // index, or the default bounds if the index is not found.
  Bounds get_bounds_(const IndexArray& indices, const BoundsArray& bounds,
                     int index, Real default_min, Real default_max) const {
    EKAT_ASSERT(index >= 0);
    int* begin = indices.data();
    int* end = begin + indices.extent(0);
    int* iter = std::lower_bound(begin, end, index);
    if ((iter == end) || (*iter != index)) {
      return Bounds({default_min, default_max});
    } else {
      return bounds(*iter);
    }
  }

  // This helper returns true if the mass mixing ratios in the given
  // array fall within the bounds in the given indices/bounds arrays.
  KOKKOS_INLINE_FUNCTION
  bool mmrs_are_valid_(const SpeciesColumnView& mmrs, const IndexArray& indices,
                       const BoundsArray& bounds) const {
    int violations = 0;
    for (int p = 0; p < indices.extent(0); ++p) {
      int index = indices(p);
      Real min = bounds(index).first;
      Real max = bounds(index).second;
      int num_mmrs = mmrs.extent(0);
      int num_levels = mmrs.extent(1);
      if (index < num_mmrs) {
        Kokkos::parallel_reduce(
            "RegionOfValidity::contains(aero)", num_levels,
            KOKKOS_LAMBDA(const int k, int& violation) {
              const auto& q = mmrs(index, k);
              auto invalid_mmr = haero::MaskType((q < min) or (q > max));
              violation += invalid_mmr.any();
            },
            violations);
      }
    }
    return (violations == 0);
  }

  // This helper returns true if the number concentrations in the given array
  // fall within the bounds in the given indices/bounds arrays.
  KOKKOS_INLINE_FUNCTION
  bool num_concs_are_valid_(const ModalColumnView& concs,
                            const IndexArray& indices,
                            const BoundsArray& bounds) const {
    int violations = 0;
    for (int p = 0; p < indices.extent(0); ++p) {
      int index = indices(p);
      Real min = bounds(index).first;
      Real max = bounds(index).second;
      int num_concs = concs.extent(0);
      int num_levels = concs.extent(1);
      if (index < num_concs) {
        Kokkos::parallel_reduce(
            "RegionOfValidity::contains(aero)", num_levels,
            KOKKOS_LAMBDA(const int k, int& violation) {
              const auto& q = concs(index, k);
              auto invalid_n = haero::MaskType((q < min) or (q > max));
              violation += invalid_n.any();
            },
            violations);
      }
    }
    return (violations == 0);
  }

  /// Minimum and maximum bounds on specific gas species, indexed by (case-
  /// insensitive) symbols. Arrays are sorted in ascending index order.
  IndexArray int_aero_indices_, cld_aero_indices_;
  BoundsArray int_aero_bounds_, cld_aero_bounds_;
  IndexArray int_n_indices_, cld_n_indices_;
  BoundsArray int_n_bounds_, cld_n_bounds_;
  IndexArray gas_indices_;
  BoundsArray gas_bounds_;
};

}  // namespace haero

#endif
