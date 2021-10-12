#include "haero/modal_aerosol_config.hpp"

#include <iomanip>
#include <sstream>

#include "haero/utils.hpp"

namespace haero {

// This sets mode->species indexing. Throws an exception if the mode_species
// mapping produces an inconsistent configuration.
void ModalAerosolConfig::index_mode_species(
    const std::map<std::string, std::vector<std::string>>& mode_species) {
  species_for_mode_.clear();
  for (const auto& kv_pair : mode_species) {
    const std::string& mode_name = kv_pair.first;
    const std::vector<std::string>& aero_species = kv_pair.second;

    // Make sure an array of species exists for the mode with the given name.
    auto m_iter =
        std::find_if(aerosol_modes.begin(), aerosol_modes.end(),
                     [&](auto mode) { return mode.name() == mode_name; });
    if (m_iter == aerosol_modes.end()) continue;
    auto mode_index = std::distance(aerosol_modes.begin(), m_iter);
    species_for_mode_.resize(
        std::max(species_for_mode_.size(), size_t(mode_index + 1)));

    // Place the appropriate species into the mode array.
    for (const auto& species_name : aero_species) {
      auto s_iter = std::find_if(
          aerosol_species.begin(), aerosol_species.end(),
          [&](auto species) { return species.symbol() == species_name; });
      if (s_iter == aerosol_species.end()) continue;
      auto species_index = std::distance(aerosol_species.begin(), s_iter);
      species_for_mode_[mode_index].push_back(species_index);
    }
  }

  // Make sure each mode contains at least one species.
  for (int m = 0; m < species_for_mode_.size(); ++m) {
    EKAT_REQUIRE_MSG(not species_for_mode_[m].empty(),
                     aerosol_modes[m].name().c_str()
                         << " mode contains no aerosol species!");
  }
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
    ss << std::setw(20) << aerosol_modes[m].name() << std::setw(7)
       << species_for_mode_[m].size() << "  ";
    const auto species_for_mode = aerosol_species_for_mode(m);
    for (int s = 0; s < species_for_mode.size(); ++s) {
      ss << species_for_mode[s].symbol()
         << (s < species_for_mode.size() - 1 ? ", " : "\n");
    }
    ss << "\n";
  }
  return ss.str();
}

ModalAerosolConfig create_simple_test_config() {
  const int nmodes = 2;
  const std::vector<std::string> mode_names = {"test_mode0", "test_mode1"};
  // Mode's minimum, nominal and maximum diameters [m]
  const auto mode_diams = std::vector<std::tuple<Real, Real, Real>>{
      std::make_tuple(1e-7, 5e-7, 1e-6), std::make_tuple(1e-6, 4e-6, 1e-5)};
  const std::vector<Real> mode_sigmas = {1, 1.5};
  const Real rh_deliq = 0.8;
  const Real rh_cryst = 0.35;

  std::map<std::string, std::vector<std::string>> mode_spec_map;
  std::vector<Mode> test_modes(nmodes);
  for (int m = 0; m < nmodes; ++m) {
    test_modes[m] = Mode(mode_names[m], std::get<0>(mode_diams[m]),
                         std::get<1>(mode_diams[m]), std::get<2>(mode_diams[m]),
                         mode_sigmas[m], rh_cryst, rh_deliq);
    mode_spec_map.emplace(mode_names[m], std::vector<std::string>());
  }

  const int nspec = 2;
  const std::vector<std::string> spec_names = {"test_spec0", "test_spec1"};
  const std::vector<std::string> spec_symbs = {"TS0", "TS1"};
  const Real g_to_kg = 1e-3;
  const std::vector<Real> spec_molec_weights = {g_to_kg * 10, g_to_kg * 100};
  const std::vector<Real> spec_hygro = {0.5, 1};
  const std::vector<Real> spec_dens = {1e3, 2e3};
  std::vector<AerosolSpecies> aeros(2);
  for (int s = 0; s < nspec; ++s) {
    aeros[s] =
        AerosolSpecies(spec_names[s], spec_symbs[s], "test_aerosol",
                       spec_molec_weights[s], spec_dens[s], spec_hygro[s]);
    mode_spec_map[mode_names[1]].push_back(spec_symbs[s]);
  }
  mode_spec_map[mode_names[0]].push_back(spec_symbs[0]);

  const std::vector<GasSpecies> test_gases = {
      GasSpecies("test_gas0", "TG0", "test_gas", g_to_kg * 40)};

  return ModalAerosolConfig(test_modes, aeros, mode_spec_map, test_gases);
}

}  // namespace haero
