#include "haero/modal_aerosol_config.hpp"

#include <iomanip>
#include <sstream>

#include "haero/utils.hpp"

namespace haero {

// This sets mode->species indexing. Throws an exception if the mode_species
// mapping produces an inconsistent configuration.
void ModalAerosolConfig::index_modal_species(
    const std::map<std::string, std::vector<std::string>>& mode_species) {
  HostType::view_1d<const Mode> aerosol_modes = h_aerosol_modes;
  std::vector<std::vector<int>> species_for_mode(aerosol_modes.size());

  h_n_species_per_mode = Kokkos::create_mirror_view(d_n_species_per_mode);

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
    h_n_species_per_mode(m) = species_for_mode[m].size();
  }
  Kokkos::deep_copy(d_n_species_per_mode, h_n_species_per_mode);

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

std::string ModalAerosolConfig::info_string(const int tab_level) const {
  auto tabstr = indent_string(tab_level);
  std::ostringstream ss;
  ss << tabstr << "ModalAerosolConfig info\n";
  tabstr += "\t";
  ss << tabstr << "        Mode    "
     << "nspec  "
     << "species"
     << "\n";
  for (int m = 0; m < num_modes(); ++m) {
    ss << std::setw(20) << h_aerosol_modes(m).name() << std::setw(7)
       << h_n_species_per_mode(m) << "  ";
    const auto species_for_mode = aerosol_species_for_mode(m);
    for (int s = 0; s < h_n_species_per_mode(m); ++s) {
      ss << species_for_mode[s].symbol()
         << (s < h_n_species_per_mode(m) - 1 ? ", " : "\n");
    }
    ss << "\n";
  }
  return ss.str();
}

ModalAerosolConfig create_simple_test_config() {
  const int nmodes = 2;
  const std::vector<std::string> mode_names = {"test_mode0", "test_mode1"};
  const std::vector<std::pair<Real, Real>> mode_min_max_diams = {
      std::make_pair(1e-8, 1e-7), std::make_pair(1e-7, 1e-6)};
  const std::vector<Real> mode_sigmas = {1, 1.5};
  const Real rh_deliq = 0.8;
  const Real rh_cryst = 0.35;

  std::map<std::string, std::vector<std::string>> mode_spec_map;
  std::vector<Mode> test_modes(nmodes);
  for (int m = 0; m < nmodes; ++m) {
    test_modes[m] =
        Mode(mode_names[m], mode_min_max_diams[m].first,
             mode_min_max_diams[m].second, mode_sigmas[m], rh_cryst, rh_deliq);
    mode_spec_map.emplace(mode_names[m], std::vector<std::string>());
  }

  const int nspec = 2;
  const std::vector<std::string> spec_names = {"test_spec0", "test_spec1"};
  const std::vector<std::string> spec_symbs = {"TS0", "TS1"};
  const Real g_to_kg = 1e-3;
  const std::vector<Real> spec_molec_weights = {g_to_kg * 10, g_to_kg * 100};
  const std::vector<Real> spec_dry_radius = {5e-8, 5e-7};
  const std::vector<Real> spec_hygro = {0.5, 1};
  const std::vector<Real> spec_dens = {1e3, 2e3};
  std::vector<AerosolSpecies> aeros(2);
  for (int s = 0; s < nspec; ++s) {
    aeros[s] = AerosolSpecies(spec_names[s], spec_symbs[s], "test_aerosol",
                              spec_molec_weights[s],
                              spec_dens[s], spec_hygro[s]);
    mode_spec_map[mode_names[1]].push_back(spec_symbs[s]);
  }
  mode_spec_map[mode_names[0]].push_back(spec_symbs[0]);

  const std::vector<GasSpecies> test_gases = {
      GasSpecies("test_gas0", "TG0", "test_gas", g_to_kg * 40)};

  return ModalAerosolConfig(test_modes, aeros, mode_spec_map, test_gases);
}

}  // namespace haero
