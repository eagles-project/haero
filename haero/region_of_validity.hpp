#ifndef HAERO_REGION_OF_VALIDITY_HPP
#define HAERO_REGION_OF_VALIDITY_HPP

#include <algorithm>

#include "haero/atmosphere.hpp"
#include "haero/conversions.hpp"
#include "haero/haero.hpp"
#include "haero/ideal_gas.hpp"
#include "haero/modal_aerosol_config.hpp"
#include "haero/prognostics.hpp"

namespace haero {

/// @class RegionOfValidity
/// This type expresses the "region of validity" for a given aerosol process
/// in terms of bounds on atmospheric and aerosol prognostic states. A region
/// of validity is initialized on the host and subsequently accessed from the
/// device via an aerosol process.
class RegionOfValidity final {
 public:
  using Bounds = Kokkos::pair<Real, Real>;

  /// Constructor
  RegionOfValidity()
      : temp_bounds({0, 500}),
        rel_hum_bounds({0, 1}),
        num_modes_(0),
        num_aero_species_in_mode_("Number of aerosol species per mode", 0),
        aero_species_("Aerosol species data", 0),
        gas_species_("Gas species data", 0),
        int_aero_indices_("Interstitial aerosol indices", 0),
        cld_aero_indices_("Cloudborne aerosol indices", 0),
        int_aero_bounds_("Interstitial aerosol bounds", 0),
        cld_aero_bounds_("Cloudborne aerosol bounds", 0),
        gas_indices_("Gas indices", 0),
        gas_bounds_("Gas bounds", 0) {}

  /// Destructor.
  ~RegionOfValidity() {}

  /// Minimum and maximum bounds on atmospheric temperature [K]
  Bounds temp_bounds;
  /// Minimum and maximum bounds on relative humidity [-]
  Bounds rel_hum_bounds;

  /// On host: initializes a RegionOfValidity with aerosol configuration
  /// metadata.
  void init(const ModalAerosolConfig& config) {
    // Gather information from the config.
    num_modes_ = config.num_modes();
    int num_pops = config.num_aerosol_populations;
    DeviceType::view_1d<int> host_num_aero_species_in_mode("", num_modes_);
    DeviceType::view_1d<AerosolSpecies> host_aero_species("", num_pops);
    for (int m = 0; m < num_modes_; ++m) {
      auto aero_species = config.aerosol_species_for_mode(m);
      host_num_aero_species_in_mode(m) = aero_species.size();
      for (int s = 0; s < aero_species.size(); ++s) {
        int pop_index = config.population_index(m, s);
        host_aero_species(pop_index) = aero_species[s];
      }
    }

    int num_gases = config.num_gases();
    DeviceType::view_1d<GasSpecies> host_gas_species("", num_gases);
    for (int g = 0; g < num_gases; ++g) {
      host_gas_species(g) = config.gas_species[g];
    }

    // Copy the information over to the device.
    Kokkos::resize(num_aero_species_in_mode_, host_num_aero_species_in_mode.extent(0));
    Kokkos::deep_copy(num_aero_species_in_mode_, host_num_aero_species_in_mode);
    Kokkos::resize(aero_species_, host_aero_species.extent(0));
    Kokkos::deep_copy(aero_species_, host_aero_species);
    Kokkos::resize(gas_species_, host_gas_species.extent(0));
    Kokkos::deep_copy(gas_species_, host_gas_species);
  }

  /// On host: adds a set of bounds for the number concentration of the
  /// interstitial aerosol species with the given population index.
  /// @param [in] pop_index The index of the aerosol species within a mode as
  ///                       indexed by a population index in a prognostics
  ///                       container
  /// @param [in] min_n The minimum valid number concentration [m-3]
  /// @param [in] max_n The maximum valid number concentration [m-3]
  void set_interstitial_aerosol_bounds(int pop_index, Real min_n, Real max_n) {
    insert_bounds_at_index_(pop_index, min_n, max_n, int_aero_indices_,
                            int_aero_bounds_);
  }

  /// On host: Adds a set of bounds for the number concentration of the
  /// cloudborne aerosol species with the given population index.
  /// @param [in] pop_index The index of the aerosol species within a mode as
  ///                       indexed by a population index in a prognostics
  ///                       container
  /// @param [in] min_n The minimum valid number concentration [m-3]
  /// @param [in] max_n The maximum valid number concentration [m-3]
  void set_cloud_aerosol_bounds(int pop_index, Real min_n, Real max_n) {
    insert_bounds_at_index_(pop_index, min_n, max_n, int_aero_indices_,
                            int_aero_bounds_);
  }

  /// On host: adds a set of bounds for the number concentration of the gas
  /// species with the given index.
  /// @param [in] gas_index The index of the gas as it appears in the set of
  ///                       gases in a prognostics container
  /// @param [in] min_n The minimum valid number concentration [m-3]
  /// @param [in] max_n The maximum valid number concentration [m-3]
  void set_gas_bounds(int gas_index, Real min_n, Real max_n) {
    insert_bounds_at_index_(gas_index, min_n, max_n, gas_indices_, gas_bounds_);
  }

  /// On host: Minimum and maximum bounds on the number concentration for the
  /// interstitial aerosol corresponding to the given population index.
  /// @param [in] pop_index The index of the aerosol species within a mode as
  ///                       indexed by a population index in a prognostics
  ///                       container
  const Bounds interstitial_aerosol_bounds(int pop_index) const {
    return get_bounds_(int_aero_indices_, int_aero_bounds_, pop_index, 0.0,
                       1.0);
  }

  /// On host: Minimum and maximum bounds on the number concentration for the
  /// cloudborne aerosol corresponding to the given population index.
  /// @param [in] pop_index The index of the aerosol species within a mode as
  ///                       indexed by a population index in a prognostics
  ///                       container
  const Bounds cloud_aerosol_bounds(int pop_index) const {
    return get_bounds_(cld_aero_indices_, cld_aero_bounds_, pop_index, 0.0,
                       1.0);
  }

  /// Minimum and maximum bounds on the number concentration for the gas
  /// species corresponding to the given index. Callable from host only.
  /// @param [in] gas_index The index of the gas as it appears in the set of
  ///                       gases in a prognostics container
  const Bounds gas_bounds(int gas_index) const {
    return get_bounds_(gas_indices_, gas_bounds_, gas_index, 0.0, 1.0);
  }

  /// On device: returns true if all state data in the given atmospheric and
  /// prognostic state fall within this region of validity, false otherwise.
  /// @param [in] atmosphere The atmospheric state under consideration
  /// @param [in] prognostics The prognostic aerosol state under consideration
  KOKKOS_INLINE_FUNCTION
  bool contains(const Atmosphere& atmosphere,
                const Prognostics& prognostics) const {
    using namespace haero::conversions;
    // Check atmospheric thresholds.
    int violations = 0;
    int num_levels = atmosphere.temperature.extent(0);
    Kokkos::parallel_reduce(
        "RegionOfValidity::contains(atm)", num_levels,
        KOKKOS_LAMBDA(const int k, int& violation) {
          const auto& T = atmosphere.temperature(k);
          auto invalid_T = haero::MaskType((T < temp_bounds.first) or
                                           (T > temp_bounds.second));
          const auto& qv = atmosphere.vapor_mixing_ratio(k);
          const auto& p = atmosphere.pressure(k);
          const auto RH =
              conversions::relative_humidity_from_vapor_mixing_ratio(qv, p, T);
          auto invalid_RH = haero::MaskType((RH < rel_hum_bounds.first) or
                                            (RH > rel_hum_bounds.second));
          violation += (invalid_T.any() or invalid_RH.any());
        },
        violations);
    if (violations == 0) {
      // Check the prognostic state.
      return (aero_n_valid_(prognostics.interstitial_aerosols,
                            atmosphere.vapor_mixing_ratio, atmosphere.pressure,
                            atmosphere.temperature, int_aero_indices_,
                            int_aero_bounds_) and
              aero_n_valid_(prognostics.cloud_aerosols,
                            atmosphere.vapor_mixing_ratio, atmosphere.pressure,
                            atmosphere.temperature, cld_aero_indices_,
                            cld_aero_bounds_) and
              gas_n_valid_(prognostics.gases,
                           atmosphere.vapor_mixing_ratio, atmosphere.pressure,
                           atmosphere.temperature, gas_indices_, gas_bounds_));
    }
    return (violations == 0);
  }

  /// On host: returns the intersection of two regions of validity r1 and r2.
  static RegionOfValidity intersection(const RegionOfValidity& r1,
                                       const RegionOfValidity& r2) {
    RegionOfValidity int_rov;
    int_rov.temp_bounds.first =
        std::max(r1.temp_bounds.first, r2.temp_bounds.first);
    int_rov.temp_bounds.second =
        std::min(r1.temp_bounds.second, r2.temp_bounds.second);
    int_rov.rel_hum_bounds.first =
        std::max(r1.rel_hum_bounds.first, r2.rel_hum_bounds.first);
    int_rov.rel_hum_bounds.second =
        std::min(r1.rel_hum_bounds.second, r2.rel_hum_bounds.second);
    intersect_bounds_(r1.int_aero_indices_, r1.int_aero_bounds_,
                      r2.int_aero_indices_, r2.int_aero_bounds_,
                      int_rov.int_aero_indices_, int_rov.int_aero_bounds_);
    intersect_bounds_(r1.cld_aero_indices_, r1.cld_aero_bounds_,
                      r2.cld_aero_indices_, r2.cld_aero_bounds_,
                      int_rov.cld_aero_indices_, int_rov.cld_aero_bounds_);
    intersect_bounds_(r1.gas_indices_, r1.gas_bounds_, r2.gas_indices_,
                      r2.gas_bounds_, int_rov.gas_indices_,
                      int_rov.gas_bounds_);
    return int_rov;
  }

 private:
  using IndexArray = kokkos_device_type::view_1d<int>;
  using BoundsArray = kokkos_device_type::view_1d<Bounds>;

  // This helper inserts a new entry into an IndexArray/BoundsArray pair at
  // the appropriate (sorted) position.
  static void insert_bounds_at_index_(int index, Real min, Real max,
                                      IndexArray& indices,
                                      BoundsArray& bounds) {
    EKAT_ASSERT(index >= 0);
    EKAT_ASSERT(min < max);

    int* begin = indices.data();
    int* end = begin + indices.extent(0);
    int* iter = std::lower_bound(begin, end, index);
    int pos = (iter == end) ? indices.extent(0) : *iter;
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

  // This helper returns true if the aerosol mass mixing ratios in the given
  // array correspond to number concentrations that fall within the bounds in
  // the given indices/bounds arrays.
  KOKKOS_INLINE_FUNCTION
  bool aero_n_valid_(const SpeciesColumnView& mmrs, const ColumnView& qv,
                     const ColumnView& p, const ColumnView& T,
                     const IndexArray& indices,
                     const BoundsArray& bounds) const {
    using namespace haero::conversions;
    int violations = 0;
    for (int i = 0; i < indices.extent(0); ++i) {
      int pop_index = indices(i);
      Real n_min = bounds(pop_index).first;
      Real n_max = bounds(pop_index).second;
      int num_pop = mmrs.extent(0);
      int num_levels = mmrs.extent(1);
      if (pop_index < num_pop) {
        const Real mu = aero_species_(pop_index).molecular_weight;
        Kokkos::parallel_reduce(
            "RegionOfValidity::aero_n_valid_", num_levels,
            KOKKOS_LAMBDA(const int k, int& violation) {
              const auto& q = mmrs(pop_index, k);
              auto SH = specific_humidity_from_vapor_mixing_ratio(qv(k));
              auto rho = ideal_gas::mass_density(p(k), T(k), qv(k));
              auto rho_d = (1 - SH) * rho;
              auto n = number_conc_from_mmr(q, mu, rho_d);
              auto invalid_n = haero::MaskType((n < n_min) or (n > n_max));
              violation += invalid_n.any();
            },
            violations);
      }
    }
    return (violations == 0);
  }

  // This helper returns true if the gas mass mixing ratios in the given
  // array correspond to number concentrations that fall within the bounds in
  // the given indices/bounds arrays.
  KOKKOS_INLINE_FUNCTION
  bool gas_n_valid_(const SpeciesColumnView& mmrs, const ColumnView& qv,
                    const ColumnView& p, const ColumnView& T,
                    const IndexArray& indices,
                    const BoundsArray& bounds) const {
    using namespace haero::conversions;
    int violations = 0;
    for (int i = 0; i < indices.extent(0); ++i) {
      int index = indices(i);
      Real n_min = bounds(index).first;
      Real n_max = bounds(index).second;
      int num_gases = mmrs.extent(0);
      int num_levels = mmrs.extent(1);
      if (index < num_gases) {
        const auto& species = gas_species_(index);
        const Real mu = species.molecular_weight;
        Kokkos::parallel_reduce(
            "RegionOfValidity::gas_n_valid_", num_levels,
            KOKKOS_LAMBDA(const int k, int& violation) {
              const auto& q = mmrs(index, k);
              auto SH = specific_humidity_from_vapor_mixing_ratio(qv(k));
              auto rho = ideal_gas::mass_density(p(k), T(k), qv(k));
              auto rho_d = (1 - SH) * rho;
              auto n = number_conc_from_mmr(q, mu, rho_d);
              auto invalid_n = haero::MaskType((n < n_min) or (n > n_max));
              violation += invalid_n.any();
            },
            violations);
      }
    }
    return (violations == 0);
  }

  // This helper generates the intersection of the bounds between two regions
  // of validity.
  static void intersect_bounds_(const IndexArray& i1, const BoundsArray& b1,
                                const IndexArray& i2, const BoundsArray& b2,
                                IndexArray& int_i, BoundsArray& int_b) {
    int_i = i1;
    int_b = b1;

    for (int i = 0; i < i2.extent(0); ++i) {
      auto index = i2(i);
      const Bounds& bounds2 = b2(i);
      auto* begin = int_i.data();
      auto* end = begin + int_i.extent(0);
      auto* iter = std::lower_bound(begin, end, index);
      if ((iter == end) || (*iter != index)) {  // bounds not found--insert
        insert_bounds_at_index_(index, bounds2.first, bounds2.second, int_i,
                                int_b);
      } else {  // overwrite bounds
        Bounds& int_bounds = int_b(index);
        int_bounds.first = std::max(bounds2.first, int_bounds.first);
        int_bounds.second = std::min(bounds2.second, int_bounds.second);
      }
    }
  }

  // Aerosol and gas species in the aerosol configuration.
  int num_modes_;
  DeviceType::view_1d<int> num_aero_species_in_mode_;
  DeviceType::view_1d<AerosolSpecies> aero_species_;
  DeviceType::view_1d<GasSpecies> gas_species_;

  /// Minimum and maximum bounds on specific aerosol and gas species, indexed
  /// by (case-insensitive) symbols. Elements in arrays are sorted in ascending
  /// index order.
  IndexArray int_aero_indices_, cld_aero_indices_;
  BoundsArray int_aero_bounds_, cld_aero_bounds_;
  IndexArray gas_indices_;
  BoundsArray gas_bounds_;
};

}  // namespace haero

#endif
