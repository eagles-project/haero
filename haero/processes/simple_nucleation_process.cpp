#include "haero/processes/simple_nucleation_process.hpp"

#include <cmath>

namespace haero {

SimpleNucleationProcess::SimpleNucleationProcess()
    : AerosolProcess(NucleationProcess, "SimpleNucleationProcess"),
      nucleation_rate_factor(1),
      pbl_factor(1),
      tendency_factor(1),
      nucleation_method(2),
      pbl_method(0) {}

//------------------------------------------------------------------------
//                                Accessors
//------------------------------------------------------------------------

void SimpleNucleationProcess::init(const ModalAerosolConfig &config) {
  // Set indices for modes, species, and gases.
  imode = config.aerosol_mode_index(nucleation_mode, false);
  iaer_so4 = config.aerosol_species_index(imode, "SO4", false);
  igas_h2so4 = config.gas_index("H2SO4", false);
  igas_nh3 = config.gas_index("NH3", false);

  // Set mode diameters.
  Kokkos::resize(d_mean_aer, config.num_modes());
  Kokkos::resize(d_min_aer, config.num_modes());
  Kokkos::resize(d_max_aer, config.num_modes());
  {
    auto d = Kokkos::create_mirror_view(d_mean_aer);
    for (int m = 0; m < config.num_modes(); ++m)
      d(m) = config.h_aerosol_modes(m).mean_std_dev;
    Kokkos::deep_copy(d_mean_aer, d);
  }
  {
    auto d = Kokkos::create_mirror_view(d_min_aer);
    for (int m = 0; m < config.num_modes(); ++m)
      d(m) = config.h_aerosol_modes(m).min_diameter;
    Kokkos::deep_copy(d_min_aer, d);
  }
  {
    auto d = Kokkos::create_mirror_view(d_max_aer);
    for (int m = 0; m < config.num_modes(); ++m)
      d(m) = config.h_aerosol_modes(m).max_diameter;
    Kokkos::deep_copy(d_max_aer, d);
  }
}

void SimpleNucleationProcess::set_param(const std::string &name, Real value) {
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

void SimpleNucleationProcess::set_param(const std::string &name, int value) {
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

void SimpleNucleationProcess::set_param(const std::string &name,
                                        const std::string &value) {
  if ("nucleation_mode" == name) {
    if (value == "") {
      EKAT_REQUIRE_MSG(false, "Invalid mode name: '" << value << "'");
    } else {
      nucleation_mode = value;
    }
  } else {
    EKAT_REQUIRE_MSG(false, "Invalid parameter: " << name);
  }
}

}  // namespace haero
