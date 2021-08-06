#ifndef HAERO_MODAL_AEROSOL_CONFIG_HPP
#define HAERO_MODAL_AEROSOL_CONFIG_HPP

#include <algorithm>
#include <map>
#include <numeric>
#include <string>

#include "haero/aerosol_species.hpp"
#include "haero/gas_species.hpp"
#include "haero/mode.hpp"
#include "haero/view_pack_helpers.hpp"

namespace haero {

/// @class ModalAerosolConfig
/// This type represents the configuration of mode-based aerosols and gas
/// species that define a physical system of aerosols in the atmosphere. It's
/// used to convey information to aerosol processes that run within a
/// \ref Model instance. Objects of this type can only be created on the host--
/// no ModalAerosolConfig object may be created on the device.
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
  ModalAerosolConfig(
      const std::vector<Mode>& aerosol_modes,
      const std::vector<AerosolSpecies>& aerosol_species,
      const std::map<std::string, std::vector<std::string>>& mode_species,
      const std::vector<GasSpecies>& gas_species)
      : aerosol_modes(aerosol_modes),
        aerosol_modes_size(aerosol_modes.size()),
        aerosol_species(aerosol_species),
        aerosol_species_size(aerosol_species.size()),
        gas_species(gas_species),
        gas_species_size(gas_species.size()),
        // Sum the length of all vectors in the mode_species map.
        // Is the use of std::accumulate here too obscure?
        num_aerosol_populations(std::accumulate(
            mode_species.begin(), mode_species.end(), 0,
            [](int i, std::pair<std::string, std::vector<std::string>> m) {
              return i + m.second.size();
            })) {
    EKAT_REQUIRE_MSG(aerosol_modes.size(),
                     "ModalAerosolConfig: No modes were defined!");
    EKAT_REQUIRE_MSG(aerosol_species.size(),
                     "ModalAerosolConfig: No aerosol species were given!");
    EKAT_REQUIRE_MSG(gas_species.size(),
                     "ModalAerosolConfig: No gas species were given!");
    index_mode_species(mode_species);
  }

  /// Copy constructor.
  ModalAerosolConfig(const ModalAerosolConfig&) = default;

  /// Destructor.
  ~ModalAerosolConfig() = default;

  /// Assignment operator.
  ModalAerosolConfig& operator=(const ModalAerosolConfig&) = default;

  /// The list of aerosol modes associated with this aerosol model.
  std::vector<Mode> aerosol_modes;
  int aerosol_modes_size;

  /// The list of all aerosol species associated with this aerosol model.
  std::vector<AerosolSpecies> aerosol_species;
  int aerosol_species_size;

  /// The list of gas species associated with this aerosol model.
  std::vector<GasSpecies> gas_species;
  int gas_species_size;

  /// The total number of distinct aerosol species populations in the
  /// model, counting appearances of one species in different modes separately.
  int num_aerosol_populations;

  KOKKOS_INLINE_FUNCTION
  int num_modes() const { return aerosol_modes_size; }

  KOKKOS_INLINE_FUNCTION
  int num_gases() const { return gas_species_size; }

  inline int max_species_per_mode() const {
    int result = 0;
    for (int m = 0; m < num_modes(); ++m) {
      const auto mode_species = aerosol_species_for_mode(m);
      result = std::max(result, int(mode_species.size()));
    }
    return result;
  }

  /// Returns the list of aerosol species associated with the model with the
  /// given mode index.
  /// @param [in] mode_index An integer index identifying the mode in question.
  ///                        This index goes from 0 to num_modes-1.
  std::vector<AerosolSpecies> aerosol_species_for_mode(
      const int mode_index) const {
    EKAT_ASSERT(mode_index >= 0);
    EKAT_ASSERT(mode_index < species_for_mode_.size());
    std::vector<AerosolSpecies> species;
    for (int s = 0; s < species_for_mode_[mode_index].size(); ++s) {
      species.push_back(aerosol_species[species_for_mode_[mode_index][s]]);
    }
    return species;
  }

  /// Returns the index of a specific aerosol mode, or -1 if the desired mode
  /// is not found.
  /// @param [in] mode_name The name of the mode for which the index is
  ///                       retrieved
  /// @param [in] case_sensitive If true, mode_name must exactly match the name
  ///                            of the aerosol mode. Otherwise, a case-
  ///                            insensitive comparison is made.
  int aerosol_mode_index(const std::string& mode_name,
                         bool case_sensitive = true) const {
    for (int m = 0; m < aerosol_modes.size(); ++m) {
      const auto& mode = aerosol_modes[m];
      if ((mode.name() == mode_name) or
          (not case_sensitive and
           (strcasecmp(mode.name().c_str(), mode_name.c_str()) == 0))) {
        return m;
      }
    }
    return -1;
  }

  /// Returns the index of a specific aerosol species within the
  /// aerosol mode with the given index, or -1 if the desired mode is not found.
  /// @param [in] mode_index The index of the mode in which the aerosol is
  /// sought.
  /// @param [in] aerosol_symbol The symbolic name of the aerosol species for
  ///                            which the index is retrieved within the given
  ///                            mode.
  /// @param [in] case_sensitive If true, aerosol_symbol must exactly match the
  ///                            aerosol's symbol. Otherwise, a case-insensitive
  ///                            comparison is made.
  int aerosol_species_index(int mode_index, const std::string& aerosol_symbol,
                            bool case_sensitive = true) const {
    EKAT_REQUIRE(mode_index >= 0);
    EKAT_REQUIRE(mode_index < species_for_mode_.size());
    for (int s = 0; s < species_for_mode_[mode_index].size(); ++s) {
      const auto& species = aerosol_species[species_for_mode_[mode_index][s]];
      if ((species.symbol() == aerosol_symbol) or
          (not case_sensitive and (strcasecmp(species.symbol().c_str(),
                                              aerosol_symbol.c_str()) == 0))) {
        return s;
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
    for (int m = 0; m < species_for_mode_.size(); ++m) {
      for (int s = 0; s < species_for_mode_[m].size(); ++s, ++p) {
        if ((m == mode_index) and (s == species_index)) {
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

  /// Returns the index of a specific gas species, or -1 if the desired species
  /// is not found.
  /// @param [in] gas_symbol The symbolic name of the gas for which the index is
  ///                        retrieved
  int gas_index(const std::string& gas_symbol,
                bool case_sensitive = true) const {
    for (int g = 0; g < gas_species.size(); ++g) {
      if ((gas_species[g].symbol() == gas_symbol) ||
          (not case_sensitive and (strcasecmp(gas_species[g].symbol().c_str(),
                                              gas_symbol.c_str()) == 0))) {
        return g;
        break;
      }
    }
    return -1;
  }

  std::string info_string(const int tab_level = 0) const;

 private:
  // The association of aerosol species with modes.
  // species_for_mode_[mode_index] = vector of species indices
  std::vector<std::vector<int>> species_for_mode_;

  // This sets mode->species indexing. Throws an exception if the mode_species
  // mapping produces an inconsistent configuration.
  void index_mode_species(
      const std::map<std::string, std::vector<std::string>>& mode_species);
};

inline ModalAerosolConfig create_mam4_modal_aerosol_config() {
  return ModalAerosolConfig(create_mam4_modes(), create_mam4_aerosol_species(),
                            create_mam4_mode_species(),
                            create_mam4_gas_species());
}

ModalAerosolConfig create_simple_test_config();

}  // namespace haero

#endif
