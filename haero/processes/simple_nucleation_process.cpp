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
      imode(-1),
      igas_h2so4(-1),
      igas_nh3(-1),
      iaer_so4(-1),
      ipop_so4(-1),
      iaer_nh4(-1),
      ipop_nh4(-1),
      d_mean_aer(0),
      d_min_aer(0),
      d_max_aer(0) {}

void SimpleNucleationProcess::init_(const ModalAerosolConfig &config) {
  // Set indices for species and gases.
  // TODO: How to do we do this without assuming a specific nucleation mode?
  iaer_so4 = config.aerosol_species_index(imode, "SO4", false);
  ipop_so4 = config.population_index(imode, iaer_so4);
  iaer_nh4 = config.aerosol_species_index(imode, "NH4", false);
  ipop_nh4 = config.population_index(imode, iaer_nh4);
  igas_h2so4 = config.gas_index("H2SO4", false);
  igas_nh3 = config.gas_index("NH3", false);

  // Set mode diameters.
  d_mean_aer.resize(config.num_modes());
  d_min_aer.resize(config.num_modes());
  d_max_aer.resize(config.num_modes());
  Kokkos::parallel_for(config.num_modes(), KOKKOS_LAMBDA(int m) {
    d_mean_aer.assign(m, config.aerosol_modes[m].mean_std_dev);
    d_min_aer.assign(m, config.aerosol_modes[m].min_diameter);
    d_max_aer.assign(m, config.aerosol_modes[m].max_diameter);
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
  if (iaer_so4 != -1) {
    for (const auto &species : config.aerosol_species_for_mode(imode)) {
      if (species.symbol() == "SO4") {
        mu_so4 = species.molecular_weight;
        break;
      }
    }
  }
  if (igas_nh3 != -1) {
    for (const auto &species : config.aerosol_species_for_mode(imode)) {
      if (species.symbol() == "NH4") {
        mu_nh4 = species.molecular_weight;
        break;
      }
    }
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
