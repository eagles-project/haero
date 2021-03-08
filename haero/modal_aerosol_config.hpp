#ifndef HAERO_MODAL_AEROSOL_CONFIG_HPP
#define HAERO_MODAL_AEROSOL_CONFIG_HPP

#include "haero/mode.hpp"
#include "haero/species.hpp"
#include <map>

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
                     const std::vector<Species>& aerosol_species,
                     const std::map<std::string, std::vector<std::string> >& mode_species,
                     const std::vector<Species>& gas_species):
    aerosol_modes(aerosol_modes), aerosol_species(aerosol_species),
    num_aerosol_populations(0), gas_species(gas_species)
  {
    EKAT_REQUIRE_MSG(not this->aerosol_modes.empty(),
      "ModalAerosolConfig: No modes were defined!");
    EKAT_REQUIRE_MSG(not this->aerosol_species.empty(),
      "ModalAerosolConfig: No aerosol species were given!");
    EKAT_REQUIRE_MSG(not this->gas_species.empty(),
      "ModalAerosolConfig: No gas species were given!");
    index_modal_species(mode_species);
  }

  /// Copy constructor.
  ModalAerosolConfig(const ModalAerosolConfig&) = default;

  /// Destructor.
  ~ModalAerosolConfig() = default;

  /// Assignment operator.
  ModalAerosolConfig& operator=(const ModalAerosolConfig&) = default;

  /// The list of aerosol modes associated with this aerosol model.
  std::vector<Mode> aerosol_modes;

  /// The list of all aerosol species associated with this aerosol
  /// model.
  std::vector<Species> aerosol_species;

  /// The total number of distinct aerosol species populations in the
  /// model, counting appearances of one species in different modes separately.
  int num_aerosol_populations;

  /// The list of gas species associated with this aerosol model.
  std::vector<Species> gas_species;

  /// Returns the list of aerosol species associated with the model with the
  /// given mode index.
  /// @param [in]mode_index An integer index identifying the mode in question. This
  ///                       This index goes from 0 to num_modes-1.
  std::vector<Species> aerosol_species_for_mode(int mode_index) const {
    EKAT_ASSERT(mode_index >= 0);
    EKAT_ASSERT(mode_index < species_for_mode_.size());
    // Construct this vector from our association data.
    std::vector<Species> species;
    for (int s = 0; s < species_for_mode_[mode_index].size(); ++s) {
      species.push_back(aerosol_species[species_for_mode_[mode_index][s]]);
    }
    return species;
  }

  private:

  // This sets mode->species indexing. Throws an exception if the mode_species
  // mapping produces an inconsistent configuration.
  void index_modal_species(const std::map<std::string, std::vector<std::string> >& mode_species) {
    species_for_mode_.resize(aerosol_modes.size());
    num_aerosol_populations = 0;
    for (auto iter = mode_species.begin(); iter != mode_species.end(); ++iter) {
      const auto& mode_name = iter->first;
      const auto& aero_species = iter->second;
      num_aerosol_populations += aero_species.size();

      auto m_iter = std::find_if(aerosol_modes.begin(), aerosol_modes.end(),
          [&] (const Mode& mode) { return mode.name == mode_name; });
      int mode_index = m_iter - aerosol_modes.begin();

      for (int s = 0; s < aero_species.size(); ++s) {
        auto s_iter = std::find_if(aerosol_species.begin(), aerosol_species.end(),
            [&] (const Species& species) {
              return species.symbol == aerosol_species[s].symbol;
            });
        int species_index = s_iter - aerosol_species.begin();
        species_for_mode_[mode_index].push_back(species_index);
      }
    }

    // Make sure each mode contains at least one species.
    for (int m = 0; m < species_for_mode_.size(); ++m) {
      EKAT_REQUIRE_MSG(not species_for_mode_[m].empty(),
        aerosol_modes[m].name.c_str() << " mode contains no aerosol species!");
    }
  }

  // The association of aerosol species with modes.
  // species_for_modes_[mode_name] = vector of species names
  std::vector<std::vector<int> > species_for_mode_;
};

}

#endif
