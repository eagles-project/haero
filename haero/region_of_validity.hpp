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
  using Bounds = std::pair<Real, Real>;

  /// Constructor
  RegionOfValidity() : temp_bounds({0, 500}), rel_hum_bounds({0, 1}) {}

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
    EKAT_ASSERT(pop_index >= 0);
    EKAT_ASSERT(min_mmr < max_mmr);

    int* begin = int_aero_indices_.data();
    int* end = begin + int_aero_indices_.extent(0);
    int* iter = std::lower_bound(begin, end, pop_index);
    if ((iter == end) || (*iter != pop_index)) {  // index must be inserted
      auto aero_indices =
          IndexArray("aero indices", int_aero_indices_.extent(0) + 1);
      auto aero_bounds =
          BoundsArray("aero bounds", int_aero_bounds_.extent(0) + 1);
      int pos = *iter;
      for (int p = aero_indices.extent(0)-1; p > pos; --p) {
        aero_indices(p) = int_aero_indices_(p - 1);
        aero_bounds(p) = int_aero_bounds_(p - 1);
      }
      aero_indices(pos) = pop_index;
      int_aero_indices_ = aero_indices;
      int_aero_bounds_ = aero_bounds;
    }
    int_aero_bounds_(pop_index).first = min_mmr;
    int_aero_bounds_(pop_index).second = max_mmr;
  }

  /// Adds a set of bounds for the cloudborne aerosol species with the given
  /// population index. Callable from host only.
  /// @param [in] pop_index The index of the aerosol species within a mode as
  ///                       indexed by a population index in a prognostics
  ///                       container
  /// @param [in] min_mmr The minimum valid mass mixing ratio [kg aero/kg air]
  /// @param [in] max_mmr The maximum valid mass mixing ratio [kg aero/kg air]
  void set_cloud_aerosol_mmr_bounds(int pop_index, Real min_mmr, Real max_mmr) {
    EKAT_ASSERT(pop_index >= 0);
    EKAT_ASSERT(min_mmr < max_mmr);

    int* begin = cld_aero_indices_.data();
    int* end = begin + cld_aero_indices_.extent(0);
    int* iter = std::lower_bound(begin, end, pop_index);
    if ((iter == end) || (*iter != pop_index)) {  // index must be inserted
      auto aero_indices =
          IndexArray("aero indices", int_aero_indices_.extent(0) + 1);
      auto aero_bounds =
          BoundsArray("aero bounds", int_aero_bounds_.extent(0) + 1);
      int pos = *iter;
      for (int p = aero_indices.extent(0)-1; p > pos; --p) {
        aero_indices(p) = cld_aero_indices_(p - 1);
        aero_bounds(p) = cld_aero_bounds_(p - 1);
      }
      aero_indices(pos) = pop_index;
      cld_aero_indices_ = aero_indices;
      cld_aero_bounds_ = aero_bounds;
    }
    cld_aero_bounds_(pop_index).first = min_mmr;
    cld_aero_bounds_(pop_index).second = max_mmr;
  }

  /// Adds a set of bounds for the gas species with the given index. Callable
  /// from host only.
  /// @param [in] gas_index The index of the gas as it appears in the set of
  ///                       gases in a prognostics container
  /// @param [in] min_mmr The minimum valid mass mixing ratio [kg gas/kg air]
  /// @param [in] max_mmr The maximum valid mass mixing ratio [kg gas/kg air]
  void set_gas_mmr_bounds(int gas_index, Real min_mmr, Real max_mmr) {
    EKAT_ASSERT(gas_index >= 0);
    EKAT_ASSERT(min_mmr < max_mmr);

    int* begin = gas_indices_.data();
    int* end = begin + gas_indices_.extent(0);
    int* iter = std::lower_bound(begin, end, gas_index);
    if ((iter == end) || (*iter != gas_index)) {  // index must be inserted
      auto gas_indices = IndexArray("gas indices", gas_indices_.extent(0) + 1);
      auto gas_bounds = BoundsArray("gas bounds", gas_bounds_.extent(0) + 1);
      int pos = *iter;
      for (int g = gas_indices.extent(0)-1; g > pos; --g) {
        gas_indices(g) = gas_indices_(g - 1);
        gas_bounds(g) = gas_bounds_(g - 1);
      }
      gas_indices(pos) = gas_index;
      gas_indices_ = gas_indices;
      gas_bounds_ = gas_bounds;
    }
    gas_bounds_(gas_index).first = min_mmr;
    gas_bounds_(gas_index).second = max_mmr;
  }

  /// Minimum and maximum bounds on the mass mixing ratio for the cloudborne
  /// aerosol corresponding to the given population index.
  KOKKOS_INLINE_FUNCTION
  const Bounds& interstitial_aerosol_bounds(int pop_index) const {
    EKAT_KERNEL_ASSERT(pop_index >= 0);
    EKAT_KERNEL_ASSERT(pop_index < int_aero_bounds_.extent(0));
    return int_aero_bounds_(pop_index);
  }

  /// Minimum and maximum bounds on the mass mixing ratio for the cloudborne
  /// aerosol corresponding to the given population index.
  KOKKOS_INLINE_FUNCTION
  const Bounds& cloud_aerosol_bounds(int pop_index) const {
    EKAT_KERNEL_ASSERT(pop_index >= 0);
    EKAT_KERNEL_ASSERT(pop_index < cld_aero_bounds_.extent(0));
    return cld_aero_bounds_(pop_index);
  }

  /// Minimum and maximum bounds on the mass mixing ratio for the gas species
  /// corresponding to the given index.
  KOKKOS_INLINE_FUNCTION
  const Bounds& gas_bounds(int gas_index) const {
    EKAT_KERNEL_ASSERT(gas_index >= 0);
    EKAT_KERNEL_ASSERT(gas_index < gas_bounds_.extent(0));
    return gas_bounds_(gas_index);
  }

  /// Returns true if all state data in the given atmosphere container falls
  /// within this region of validity, false otherwise.
  KOKKOS_INLINE_FUNCTION
  bool contains(const Atmosphere& atmosphere) const {
    int violations = 0;
    Kokkos::parallel_reduce(
        "RegionOfValidity::contains(atm)", atmosphere.num_levels(),
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
    int violations = 0;

    // Interstitial aerosol MMRs
    for (int p = 0; p < int_aero_indices_.extent(0); ++p) {
      int pop_index = int_aero_indices_(p);
      const auto& bounds = int_aero_bounds_(pop_index);
      if (pop_index < prognostics.num_aerosol_populations()) {
        Kokkos::parallel_reduce(
            "RegionOfValidity::contains(aero)", prognostics.num_levels(),
            KOKKOS_LAMBDA(const int k, int& violation) {
              const auto& q = prognostics.interstitial_aerosols(pop_index, k);
              auto invalid_mmr =
                  haero::MaskType((q < bounds.first) or (q > bounds.second));
              violation += invalid_mmr.any();
            },
            violations);
      }
    }

    // Cloudborne aerosol MMRs
    for (int p = 0; p < cld_aero_indices_.extent(0); ++p) {
      int pop_index = cld_aero_indices_(p);
      const auto& bounds = cld_aero_bounds_(pop_index);
      if (pop_index < prognostics.num_aerosol_populations()) {
        Kokkos::parallel_reduce(
            "RegionOfValidity::contains(aero)", prognostics.num_levels(),
            KOKKOS_LAMBDA(const int k, int& violation) {
              const auto& q = prognostics.cloud_aerosols(pop_index, k);
              auto invalid_mmr =
                  haero::MaskType((q < bounds.first) or (q > bounds.second));
              violation += invalid_mmr.any();
            },
            violations);
      }
    }

    // Gas MMRs
    for (int g = 0; g < gas_indices_.extent(0); ++g) {
      int gas_index = gas_indices_(g);
      const auto& bounds = gas_bounds_(gas_index);
      if (gas_index < prognostics.num_gases()) {
        Kokkos::parallel_reduce(
            "RegionOfValidity::contains(aero)", prognostics.num_levels(),
            KOKKOS_LAMBDA(const int k, int& violation) {
              const auto& q = prognostics.gases(gas_index, k);
              auto invalid_mmr =
                  haero::MaskType((q < bounds.first) or (q > bounds.second));
              violation += invalid_mmr.any();
            },
            violations);
      }
    }
    return (violations == 0);
  }

 protected:
  /// Minimum and maximum bounds on specific gas species, indexed by (case-
  /// insensitive) symbols. Arrays are sorted in ascending index order.
  using IndexArray = kokkos_device_type::view_1d<int>;
  using BoundsArray = kokkos_device_type::view_1d<Bounds>;
  IndexArray int_aero_indices_, cld_aero_indices_;
  BoundsArray int_aero_bounds_, cld_aero_bounds_;
  IndexArray gas_indices_;
  BoundsArray gas_bounds_;
};

}  // namespace haero

#endif
