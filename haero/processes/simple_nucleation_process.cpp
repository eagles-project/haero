#include "haero/processes/simple_nucleation_process.hpp"

#include <cmath>

#include "haero/conversions.hpp"
#include "haero/processes/merikanto2007.hpp"
#include "haero/processes/vehkamaki2002.hpp"

namespace haero {

SimpleNucleationProcess::SimpleNucleationProcess()
    : DeviceAerosolProcess<SimpleNucleationProcess>(NucleationProcess, "SimpleNucleationProcess"),
      nucleation_rate_factor(1),
      pbl_factor(1),
      tendency_factor(1),
      nucleation_method(2),
      pbl_method(0),
      igas_h2so4(-1),
      igas_nh3(-1),
      num_modes(),
      iaer_so4(),
      ipop_so4(),
      iaer_nh4(),
      ipop_nh4(),
      d_mean_aer(),
      d_min_aer(),
      d_max_aer() {}

void SimpleNucleationProcess::init_(const ModalAerosolConfig &config) {
  num_modes = config.num_modes();

  // Set indices for species and gases.
  iaer_so4.resize(num_modes);
  ipop_so4.resize(num_modes);
  iaer_nh4.resize(num_modes);
  ipop_nh4.resize(num_modes);
  Kokkos::parallel_for(num_modes, KOKKOS_LAMBDA(int m) {
    iaer_so4[m] = config.aerosol_species_index(m, "SO4", false);
    ipop_so4[m] = config.population_index(m, iaer_so4[m]);
    iaer_nh4[m] = config.aerosol_species_index(m, "NH4", false);
    ipop_nh4[m] = config.population_index(m, iaer_nh4[m]);
  });
  igas_h2so4 = config.gas_index("H2SO4", false);
  igas_nh3 = config.gas_index("NH3", false);

  // Set mode diameters.
  d_mean_aer.resize(num_modes);
  d_min_aer.resize(num_modes);
  d_max_aer.resize(num_modes);
  Kokkos::parallel_for(num_modes, KOKKOS_LAMBDA(int m) {
    d_mean_aer[m] = config.aerosol_modes[m].mean_std_dev;
    d_min_aer[m] = config.aerosol_modes[m].min_diameter;
    d_max_aer[m] = config.aerosol_modes[m].max_diameter;
  });

  // Jot down the molecular weights for our gases.
  if (igas_h2so4 != -1) {
    for (int g = 0; g < config.num_gases(); ++g) {
      const auto &species = config.gas_species[g];
      if (species.symbol() == "H2SO4") {
        mu_h2so4 = species.molecular_weight;
        break;
      }
    }
  }
  if (igas_nh3 != -1) {
    for (int g = 0; g < config.num_gases(); ++g) {
      const auto &species = config.gas_species[g];
      if (species.symbol() == "NH3") {
        mu_nh3 = species.molecular_weight;
        break;
      }
    }
  }

  // Do the same for our nucleating aerosols.
  auto iter = std::find_if(config.aerosol_species.begin(),
                           config.aerosol_species.end(),
                           [&](auto species) { return species.symbol() == "SO4"; });
  if (iter != config.aerosol_species.end()) {
    const AerosolSpecies& species = *iter;
    mu_so4 = species.molecular_weight;
  }
  iter = std::find_if(config.aerosol_species.begin(),
                      config.aerosol_species.end(),
                      [&](auto species) { return species.symbol() == "NH4"; });
  if (iter != config.aerosol_species.end()) {
    const AerosolSpecies& species = *iter;
    mu_nh4 = species.molecular_weight;
  }

  // Set our region of validity based on our nucleation parameterization.
  auto &rov = region_of_validity();
  if (nucleation_method == 2) {
    rov.temp_bounds = vehkamaki2002::valid_temp_range();
    rov.rel_hum_bounds = vehkamaki2002::valid_rel_hum_range();
    if (igas_h2so4 != -1) {
      auto c_h2so4_range = vehkamaki2002::valid_c_h2so4_range();
      rov.set_gas_bounds(igas_h2so4, c_h2so4_range.first, c_h2so4_range.second);
    }
  } else {
    EKAT_ASSERT(nucleation_method == 3);
    rov.temp_bounds = merikanto2007::valid_temp_range();
    rov.rel_hum_bounds = merikanto2007::valid_rel_hum_range();
    if (igas_h2so4 != -1) {
      auto c_h2so4_range = merikanto2007::valid_c_h2so4_range();
      rov.set_gas_bounds(igas_h2so4, c_h2so4_range.first, c_h2so4_range.second);
    }
    if (igas_nh3 != -1) {
      auto xi_nh3_range = merikanto2007::valid_xi_nh3_range();  // [ppt]
      Real c_nh3_min =
          1e-12 * conversions::mmr_from_vmr(xi_nh3_range.first, mu_nh3);
      Real c_nh3_max =
          1e-12 * conversions::mmr_from_vmr(xi_nh3_range.second, mu_nh3);
      rov.set_gas_bounds(igas_nh3, c_nh3_min, c_nh3_max);
    }
  }
}

void SimpleNucleationProcess::set_param_(const std::string &name, Real value) {
  if ("pbl_factor" == name) {
    if (value > 0) {
      pbl_factor = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else if ("nucleation_rate_factor" == name) {
    if (value > 0) {
      nucleation_rate_factor = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else if ("tendency_factor" == name) {
    if (value > 0) {
      tendency_factor = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else {
    EKAT_REQUIRE_MSG(false, "Invalid parameter: " << name);
  }
}

void SimpleNucleationProcess::set_param_(const std::string &name, int value) {
  if ("nucleation_method" == name) {
    if ((value == 2) or (value == 3)) {
      nucleation_method = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else if ("pbl_method" == name) {
    if ((value == 0) or (value == 1) or (value == 2)) {
      pbl_method = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else {
    EKAT_REQUIRE_MSG(false, "Invalid parameter: " << name);
  }
}

}  // namespace haero
