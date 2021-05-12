#include "haero/modal_aerosol_config.hpp"

namespace haero {

// This sets mode->species indexing. Throws an exception if the mode_species
// mapping produces an inconsistent configuration.
void ModalAerosolConfig::index_modal_species(
    const std::map<std::string, std::vector<std::string>>& mode_species) {
  HostType::view_1d<const Mode> aerosol_modes = h_aerosol_modes;
  std::vector<std::vector<int>> species_for_mode(aerosol_modes.size());
  for (auto iter = mode_species.begin(); iter != mode_species.end(); ++iter) {
    const std::string& mode_name = iter->first;
    const std::vector<std::string>& aero_species = iter->second;

    int mode_index = -1;
    for (unsigned i = 0; i < aerosol_modes.size() && mode_index == -1; ++i)
      if (aerosol_modes[i].name() == mode_name) mode_index = i;

    for (int s = 0; s < aero_species.size(); ++s) {
      int species_index = -1;
      for (unsigned i = 0; i < h_aerosol_species.size() && species_index == -1;
           ++i)
        if (h_aerosol_species[i].symbol() == aero_species[s]) species_index = i;
      species_for_mode[mode_index].push_back(species_index);
    }
  }
  // Make sure each mode contains at least one species.
  for (int m = 0; m < species_for_mode.size(); ++m) {
    EKAT_REQUIRE_MSG(not species_for_mode[m].empty(),
                     aerosol_modes[m].name().c_str()
                         << " mode contains no aerosol species!");
  }

  // Move to a 2D View by making rectangular and filling ragged arrays with -1.
  unsigned long max_num_species = 0;
  for (int m = 0; m < species_for_mode.size(); ++m)
    max_num_species = std::max(max_num_species, species_for_mode[m].size());
  for (int m = 0; m < species_for_mode.size(); ++m)
    for (int s = species_for_mode[m].size(); s < max_num_species; ++s)
      species_for_mode[m].push_back(-1);
  d_species_for_mode = vector_to_2dview(species_for_mode, "species_for_mode");

  h_species_for_mode = Kokkos::create_mirror_view(d_species_for_mode);
  Kokkos::deep_copy(h_species_for_mode, d_species_for_mode);
}

}  // namespace haero
