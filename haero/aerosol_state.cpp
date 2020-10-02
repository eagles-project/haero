#include "haero/aerosol_state.hpp"

namespace haero {

AerosolState::AerosolState(int num_columns, int num_levels):
  num_columns_(num_columns),
  num_levels_(num_levels),
  finalized_(false) {
}

AerosolState::~AerosolState() {
}

int AerosolState::add_aerosol_mode(const Mode& mode,
                                   const std::vector<Species>& aero_species) {
}

int AerosolState::add_aerosol_mode(const Mode& mode,
                                   const std::vector<Species>& aero_species,
                                   ColumnSpeciesView& int_aero_data,
                                   ColumnSpeciesView& cld_aero_data) {
}

int AerosolState::add_gas_species(const std::vector<Species>& gas_species) {
}

int AerosolState::add_gas_species(const std::vector<Species>& gas_species,
                                  ColumnSpeciesView& int_aero_data) {
}

void AerosolState::finalize() {
  finalized_ = true;
}

bool AerosolState::is_finalized() const {
  return finalized_;
}

int AerosolState::num_aerosol_modes() const {
  return num_modes_;
}

int AerosolState::num_aerosol_species(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < num_modes_);
}

int AerosolState::num_gas_species() const {
}

AerosolState::ColumnSpeciesView&
AerosolState::interstitial_aerosols(int mode_index) {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < num_modes_);
  return int_aero_species_[mode_index];
}

const AerosolState::ColumnSpeciesView&
AerosolState::interstitial_aerosols(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < num_modes_);
  return int_aero_species_[mode_index];
}

AerosolState::ColumnSpeciesView&
AerosolState::cloudborne_aerosols(int mode_index) {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < num_modes_);
  return cld_aero_species_[mode_index];
}

const AerosolState::ColumnSpeciesView&
AerosolState::cloudborne_aerosols(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < num_modes_);
  return cld_aero_species_[mode_index];
}

AerosolState::ColumnSpeciesView&
AerosolState::gas_mole_fractions() {
  return gas_mole_fractions_;
}

const AerosolState::ColumnSpeciesView&
AerosolState::gas_mole_fractions() const {
  return gas_mole_fractions_;
}

AerosolState::ColumnView& AerosolState::modal_densities() {
  return modal_num_densities_;
}

const AerosolState::ColumnView& AerosolState::modal_densities() const {
  return modal_num_densities_;
}

}

