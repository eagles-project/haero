#include "haero/processes/simple_nucleation_process.hpp"

#include <cmath>

#include "haero/conversions.hpp"
#include "haero/processes/merikanto2007.hpp"
#include "haero/processes/vehkamaki2002.hpp"

namespace haero {

SimpleNucleationProcess::SimpleNucleationProcess()
    : DeviceAerosolProcess<SimpleNucleationProcess>("SimpleNucleationProcess"),
      nucleation_rate_factor_(1),
      pbl_factor_(1),
      tendency_factor_(1),
      nucleation_method_(2),
      pbl_method_(0),
      accom_coeff_h2so4_(0.65),
      density_so4_(1770.0), // TODO: verify dry mass density
      mw_h2so4_(Constants::molec_weight_h2so4),
      mw_nh3_(Constants::molec_weight_nh3),
      mw_so4_(Constants::molec_weight_so4),
      mw_nh4_(Constants::molec_weight_nh4),
      igas_h2so4_(-1),
      igas_nh3_(-1),
      num_modes_(),
      inuc_mode_(-1),
      iaer_so4_(),
      ipop_so4_(),
      iaer_nh4_(),
      ipop_nh4_(),
      d_mean_aer_(),
      d_min_aer_(),
      d_max_aer_() {}

void SimpleNucleationProcess::init_(const ModalAerosolConfig &config) {
  num_modes_ = config.num_modes();

  // Find the Aitken mode, into which we place nucleated particles.
  inuc_mode_ = config.aerosol_mode_index("aitken", false);
  EKAT_REQUIRE_MSG(inuc_mode_ >= 0, "Aitken mode not found in aerosol config!");

  // Set indices for species and gases.
  Kokkos::resize(iaer_so4_, num_modes_);
  Kokkos::resize(ipop_so4_, num_modes_);
  Kokkos::resize(iaer_nh4_, num_modes_);
  Kokkos::resize(ipop_nh4_, num_modes_);
  auto h_iaer_so4 = Kokkos::create_mirror_view(iaer_so4_);
  auto h_ipop_so4 = Kokkos::create_mirror_view(ipop_so4_);
  auto h_iaer_nh4 = Kokkos::create_mirror_view(iaer_nh4_);
  auto h_ipop_nh4 = Kokkos::create_mirror_view(ipop_nh4_);
  for (int m = 0; m < num_modes_; ++m) {
    h_iaer_so4[m] = config.aerosol_species_index(m, "SO4", false);
    h_ipop_so4[m] = config.population_index(m, h_iaer_so4[m]);
    h_iaer_nh4[m] = config.aerosol_species_index(m, "NH4", false);
    h_ipop_nh4[m] = config.population_index(m, h_iaer_nh4[m]);
  }
  Kokkos::deep_copy(iaer_so4_, h_iaer_so4);
  Kokkos::deep_copy(ipop_so4_, h_ipop_so4);
  Kokkos::deep_copy(iaer_nh4_, h_iaer_nh4);
  Kokkos::deep_copy(ipop_nh4_, h_ipop_nh4);

  igas_h2so4_ = config.gas_index("H2SO4", false);
  igas_nh3_ = config.gas_index("NH3", false);

  // Set mode diameters.
  Kokkos::resize(d_mean_aer_, num_modes_);
  Kokkos::resize(d_min_aer_, num_modes_);
  Kokkos::resize(d_max_aer_, num_modes_);
  auto h_mean_aer = Kokkos::create_mirror_view(d_mean_aer_);
  auto h_min_aer = Kokkos::create_mirror_view(d_min_aer_);
  auto h_max_aer = Kokkos::create_mirror_view(d_max_aer_);
  for (int m = 0; m < num_modes_; ++m) {
    h_mean_aer[m] = config.aerosol_modes[m].mean_std_dev;
    h_min_aer[m] = config.aerosol_modes[m].min_diameter;
    h_max_aer[m] = config.aerosol_modes[m].max_diameter;
  }
  Kokkos::deep_copy(d_mean_aer_, h_mean_aer);
  Kokkos::deep_copy(d_min_aer_, h_min_aer);
  Kokkos::deep_copy(d_max_aer_, h_max_aer);

  // Jot down the molecular weights for our gases.
  if (igas_h2so4_ != -1) {
    for (int g = 0; g < config.num_gases(); ++g) {
      const auto &species = config.gas_species[g];
      if (species.symbol() == "H2SO4") {
        mw_h2so4_ = species.molecular_weight;
        break;
      }
    }
  }
  if (igas_nh3_ != -1) {
    for (int g = 0; g < config.num_gases(); ++g) {
      const auto &species = config.gas_species[g];
      if (species.symbol() == "NH3") {
        mw_nh3_ = species.molecular_weight;
        break;
      }
    }
  }

  // Do the same for our nucleating aerosols.
  auto iter =
      std::find_if(config.aerosol_species.begin(), config.aerosol_species.end(),
                   [&](auto species) { return species.symbol() == "SO4"; });
  if (iter != config.aerosol_species.end()) {
    const AerosolSpecies &species = *iter;
    mw_so4_ = species.molecular_weight;
  }
  iter =
      std::find_if(config.aerosol_species.begin(), config.aerosol_species.end(),
                   [&](auto species) { return species.symbol() == "NH4"; });
  if (iter != config.aerosol_species.end()) {
    const AerosolSpecies &species = *iter;
    mw_nh4_ = species.molecular_weight;
  }

  // Set our region of validity based on our nucleation parameterization.
  auto &rov = region_of_validity();
  if (nucleation_method_ == binary_nucleation) {
    rov.temp_bounds = vehkamaki2002::valid_temp_range();
    rov.rel_hum_bounds = vehkamaki2002::valid_rel_hum_range();
    if (igas_h2so4_ != -1) {
      auto c_h2so4_range = vehkamaki2002::valid_c_h2so4_range();
      rov.set_gas_bounds(igas_h2so4_, c_h2so4_range.first,
                         c_h2so4_range.second);
    }
  } else {
    EKAT_ASSERT(nucleation_method_ == ternary_nucleation);
    rov.temp_bounds = merikanto2007::valid_temp_range();
    rov.rel_hum_bounds = merikanto2007::valid_rel_hum_range();
    if (igas_h2so4_ != -1) {
      auto c_h2so4_range = merikanto2007::valid_c_h2so4_range();
      rov.set_gas_bounds(igas_h2so4_, c_h2so4_range.first,
                         c_h2so4_range.second);
    }
    if (igas_nh3_ != -1) {
      auto xi_nh3_range = merikanto2007::valid_xi_nh3_range();  // [ppt]
      Real c_nh3_min =
          1e-12 * conversions::mmr_from_vmr(xi_nh3_range.first, mw_nh3_);
      Real c_nh3_max =
          1e-12 * conversions::mmr_from_vmr(xi_nh3_range.second, mw_nh3_);
      rov.set_gas_bounds(igas_nh3_, c_nh3_min, c_nh3_max);
    }
  }
}

void SimpleNucleationProcess::set_param_(const std::string &name, Real value) {
  if ("pbl_factor" == name) {
    if (value > 0) {
      pbl_factor_ = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else if ("nucleation_rate_factor" == name) {
    if (value > 0) {
      nucleation_rate_factor_ = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else if ("tendency_factor" == name) {
    if (value > 0) {
      tendency_factor_ = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else if ("accom_coeff_h2so4" == name) {
    if (value > 0) {
      accom_coeff_h2so4_ = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else if ("mw_h2so4" == name) {
    if (value > 0) {
      mw_h2so4_ = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else if ("mw_nh3" == name) {
    if (value > 0) {
      mw_nh3_ = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else if ("mw_so4" == name) {
    if (value > 0) {
      mw_so4_ = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else if ("mw_nh4" == name) {
    if (value > 0) {
      mw_nh4_ = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else {
    EKAT_REQUIRE_MSG(false, "Invalid parameter: " << name);
  }
}

void SimpleNucleationProcess::set_param_(const std::string &name, int value) {
  if ("nucleation_method" == name) {
    if ((value == binary_nucleation) or (value == ternary_nucleation)) {
      nucleation_method_ = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else if ("pbl_method" == name) {
    if ((value == no_pbl_correction) or (value == first_order_pbl_correction) or
        (value == second_order_pbl_correction)) {
      pbl_method_ = value;
    } else {
      EKAT_REQUIRE_MSG(false, "Invalid " << name << ": " << value);
    }
  } else {
    EKAT_REQUIRE_MSG(false, "Invalid parameter: " << name);
  }
}

}  // namespace haero
