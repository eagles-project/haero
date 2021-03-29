#ifndef HAERO_MODAL_AEROSOL_CONFIG_HPP
#define HAERO_MODAL_AEROSOL_CONFIG_HPP

#include "haero/mode.hpp"
#include "haero/aerosol_species.hpp"
#include "haero/view_pack_helpers.hpp"
#include <map>
#include <algorithm>
#include <numeric>

namespace haero {

/// @struct ModalAerosolConfig

/// This type represents the configuration of mode-based aerosols and gas
/// species that define a physical system of aerosols in the atmosphere. It's
/// used to convey information to aerosol processes that run within a
/// \ref Model instance.
class ModalAerosolConfig final {
  public:

  /// Default constructor -- use with caution.
  ModalAerosolConfig() = default;

  /// Creates a model configuration that defines the represenation of a set of
  /// aerosols and gases in the atmosphere. This constructor only works on the
  /// host.
  /// @param [in] aerosol_modes a list of aerosol modes supported by the Context
  /// @param [in] aerosol_species a list of aerosol species supported by the
  ///                             Context
  /// @param [in] mode_species a map that defines the association of aerosol
  ///                          species with an aerosol mode. The keys in this
  ///                          map are names of aerosol modes (corresponding to
  ///                          those found in `modes`), and the values are lists
  ///                          of symbolic names (symbols) of aerosol species
  ///                          (supplied in `aerosol_species`) that belong to
  ///                          those modes.
  /// @param [in] gas_species a list of gas species supported by the Context
  ModalAerosolConfig(const std::vector<Mode>& aerosol_modes,
                     const std::vector<AerosolSpecies>& aerosol_species,
                     const std::map<std::string, std::vector<std::string> >& mode_species,
                     const std::vector<AerosolSpecies>& gas_species):
    d_aerosol_modes  (vector_to_1dview(aerosol_modes,   "aerosol_modes")),
    d_aerosol_species(vector_to_1dview(aerosol_species, "aerosol_species")),

    // Sum the length of all vectors in the mode_species map.
    // Is the use of std::accumulate here too obscure?
    num_aerosol_populations(std::accumulate(mode_species.begin(), mode_species.end(), 0,
      [](int i, std::pair<std::string, std::vector<std::string> > m){return i + m.second.size();})),
    d_gas_species    (vector_to_1dview(gas_species, "gas_species"))
  {
    h_aerosol_modes = Kokkos::create_mirror_view(d_aerosol_modes);
    Kokkos::deep_copy(h_aerosol_modes, d_aerosol_modes);

    h_aerosol_species = Kokkos::create_mirror_view(d_aerosol_species);
    Kokkos::deep_copy(h_aerosol_species, d_aerosol_species);

    h_gas_species = Kokkos::create_mirror_view(d_gas_species);
    Kokkos::deep_copy(h_gas_species, d_gas_species);

    EKAT_REQUIRE_MSG(h_aerosol_modes.size(),   "ModalAerosolConfig: No modes were defined!");
    EKAT_REQUIRE_MSG(h_aerosol_species.size(), "ModalAerosolConfig: No aerosol species were given!");
    EKAT_REQUIRE_MSG(h_gas_species.size(),     "ModalAerosolConfig: No gas species were given!");
    index_modal_species(mode_species);
  }

  /// Copy constructor.
  ModalAerosolConfig(const ModalAerosolConfig&) = default;

  /// Destructor.
  ~ModalAerosolConfig() = default;

  /// Assignment operator.
  ModalAerosolConfig& operator=(const ModalAerosolConfig&) = default;

  /// The list of aerosol modes associated with this aerosol model.
  DeviceType::view_1d<Mode> d_aerosol_modes;
  HostType::view_1d<Mode>   h_aerosol_modes;

  /// The list of all aerosol species associated with this aerosol
  /// model.
  DeviceType::view_1d<AerosolSpecies> d_aerosol_species;
  HostType::view_1d<AerosolSpecies>   h_aerosol_species;

  /// The total number of distinct aerosol species populations in the
  /// model, counting appearances of one species in different modes separately.
  int num_aerosol_populations;

  /// The list of gas species associated with this aerosol model.
  DeviceType::view_1d<AerosolSpecies> d_gas_species;
  HostType::view_1d<AerosolSpecies>   h_gas_species;

  /// Returns the list of aerosol species associated with the model with the
  /// given mode index.
  /// @param [in]mode_index An integer index identifying the mode in question. This
  ///                       This index goes from 0 to num_modes-1.
  std::vector<AerosolSpecies> aerosol_species_for_mode(const int mode_index) const {
    EKAT_ASSERT(mode_index >= 0);
    EKAT_ASSERT(mode_index < h_species_for_mode.extent(0));
    // Construct this vector from our association data.
    std::vector<AerosolSpecies> species;
    for (int s = 0; s < h_species_for_mode.extent(1) && 0 <= h_species_for_mode(mode_index,s); ++s) {
      species.push_back(h_aerosol_species[h_species_for_mode(mode_index,s)]);
    }
    return species;
  }

  /// On host: returns the index of a specific aerosol mode, or -1 if the
  /// desired mode is not found.
  /// @param [in] mode_name The name of the mode for which the index is retrieved
  int aerosol_mode_index(const std::string& mode_name) const {
    for (int m = 0; m < h_aerosol_modes.size(); ++m) {
      if (h_aerosol_modes[m].name() == mode_name) {
        return m;
        break;
      }
    }
    return -1;
  }

  /// On host: returns the index of a specific aerosol species within the
  /// aerosol mode with the given index, or -1 if the desired mode is not found.
  /// @param [in] mode_index The index of the mode in which the aerosol is sought.
  /// @param [in] aerosol_symbol The symbolic name of the aerosol species for
  ///                            which the index is retrieved within the given
  ///                            mode.
  int aerosol_species_index(int mode_index, const std::string& aerosol_name) const {
    for (int s = 0; s < h_species_for_mode.extent(1); ++s) {
      int species_index = h_species_for_mode(mode_index, s);
      if ((species_index >= 0) &&
          (h_aerosol_species[species_index].name() == aerosol_name)) {
        return s;
        break;
      }
    }
    return -1;
  }

  /// On host: returns the population index corresponding to the given aerosol
  /// mode and species indices.
  /// @param [in] mode_index The index of the aerosol mode
  /// @param [in] species_index The index of the aerosol species
  int population_index(int mode_index, int species_index) const {
    bool found = false;
    int p = 0;
    for (int m = 0; m < h_species_for_mode.extent(0); ++m) {
      for (int s = 0; s < h_species_for_mode.extent(1) && 0 <= h_species_for_mode(m, s); ++s, ++p) {
        if ((m == mode_index) && (s == species_index)) {
          found = true;
          break;
        }
      }
      if (found) {
        break;
      }
    }
    return (found) ? p : -1;
  }

  /// On host: returns the index of a specific gas species, or -1 if the
  /// desired species is not found.
  /// @param [in] gas_symbol The symbolic name of the gas for which the index is
  ///                        retrieved
  int gas_index(const std::string& gas_symbol) const {
    for (int g = 0; g < h_gas_species.size(); ++g) {
      if (h_gas_species[g].symbol() == gas_symbol) {
        return g;
        break;
      }
    }
    return -1;
  }

  /// On Device: Returns the list of aerosol species associated with the model with the
  /// given mode index.
  /// @param [in]mode_index An integer index identifying the mode in question. This
  ///                       This index goes from 0 to num_modes-1.
  /// @param [out]aerosol_species A list of f aerosol species associated with the model with the
  ///                       given mode index.  Must have been sized large enough to hold as many
  ///                       species as may be found.
  KOKKOS_INLINE_FUNCTION
  void aerosol_species_for_mode(const int mode_index, DeviceType::view_1d<AerosolSpecies> aerosol_species) const {
    EKAT_KERNEL_ASSERT(mode_index >= 0);
    EKAT_KERNEL_ASSERT(mode_index < d_species_for_mode.extent(0));
    // Construct this vector from our association data.
    for (int s = 0; s < d_species_for_mode.extent(1) && 0 <= d_species_for_mode(mode_index,s); ++s) {
      EKAT_KERNEL_ASSERT(s < aerosol_species.extent(0));
      const int species_for_mode = d_species_for_mode(mode_index,s);
      const AerosolSpecies species = d_aerosol_species[species_for_mode];
      aerosol_species[s] = species;
    }
  }

  private:
  // The association of aerosol species with modes.
  // species_for_modes_[mode_name] = vector of species names
  DeviceType::view_2d<int> d_species_for_mode;
  HostType::view_2d<int>   h_species_for_mode;

  // This sets mode->species indexing. Throws an exception if the mode_species
  // mapping produces an inconsistent configuration.
  void index_modal_species(const std::map<std::string, std::vector<std::string> >&);
};

inline ModalAerosolConfig create_mam4_modal_aerosol_config() {
  return ModalAerosolConfig(create_mam4_modes(), create_mam4_aerosol_species(),
                            create_mam4_mode_species(), create_mam4_gas_species());
}

}

#endif
