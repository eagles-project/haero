#include "haero/aero_state.hpp"

namespace haero {

AeroState::AeroState(int num_columns, int num_levels):
  num_columns_(num_columns),
  num_levels_(num_levels),
  finalized_(false) {
}

AeroState::~AeroState() {
}

int AeroState::add_aerosol_mode(const Mode& mode,
                                const std::vector<Species>& aero_species) {
  EKAT_ASSERT_MSG(not finalized_,
                  "Cannot add an aerosol mode to a finalized AeroState!");
  num_aero_species_.push_back(aero_species.size());
  return static_cast<int>(num_aero_species_.size()-1);
}

int AeroState::add_aerosol_mode(const Mode& mode,
                                const std::vector<Species>& aero_species,
                                ColumnSpeciesView& int_aero_data,
                                ColumnSpeciesView& cld_aero_data) {
  EKAT_ASSERT_MSG(not finalized_,
                  "Cannot add an aerosol mode to a finalized AeroState!");
  num_aero_species_.push_back(aero_species.size());
  return static_cast<int>(num_aero_species_.size()-1);
}

int AeroState::add_gas_species(const std::vector<Species>& gas_species) {
  EKAT_ASSERT_MSG(not finalized_,
                  "Cannot add gas species to a finalized AeroState!");
  num_gas_species_ += gas_species.size();
  return static_cast<int>(num_aero_species_.size()-gas_species.size()-1);
}

int AeroState::add_gas_species(const std::vector<Species>& gas_species,
                               ColumnSpeciesView& int_aero_data) {
  EKAT_ASSERT_MSG(not finalized_,
                  "Cannot add gas species to a finalized AeroState!");
  num_gas_species_ += gas_species.size();
  return static_cast<int>(num_aero_species_.size()-gas_species.size()-1);
}

void AeroState::finalize() {
  finalized_ = true;
}

bool AeroState::is_finalized() const {
  return finalized_;
}

int AeroState::num_aerosol_modes() const {
  return num_aero_species_.size();
}

int AeroState::num_aerosol_species(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < num_aero_species_.size());
  return num_aero_species_[mode_index];
}

int AeroState::num_gas_species() const {
  return num_gas_species_;
}

int AeroState::num_columns() const {
  return num_columns_;
}

int AeroState::num_levels() const {
  return num_levels_;
}

AeroState::ColumnSpeciesView&
AeroState::interstitial_aerosols(int mode_index) {
  EKAT_ASSERT_MSG(finalized_,
                  "Cannot access data in a non-finalized AeroState!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < num_aero_species_.size());
  return int_aero_species_[mode_index];
}

const AeroState::ColumnSpeciesView&
AeroState::interstitial_aerosols(int mode_index) const {
  EKAT_ASSERT_MSG(finalized_,
                  "Cannot access data in a non-finalized AeroState!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < num_aero_species_.size());
  return int_aero_species_[mode_index];
}

AeroState::ColumnSpeciesView& AeroState::cloudborne_aerosols(int mode_index) {
  EKAT_ASSERT_MSG(finalized_,
                  "Cannot access data in a non-finalized AeroState!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < num_aero_species_.size());
  return cld_aero_species_[mode_index];
}

const AeroState::ColumnSpeciesView&
AeroState::cloudborne_aerosols(int mode_index) const {
  EKAT_ASSERT_MSG(finalized_,
                  "Cannot access data in a non-finalized AeroState!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < num_aero_species_.size());
  return cld_aero_species_[mode_index];
}

AeroState::ColumnSpeciesView& AeroState::gas_mole_fractions() {
  EKAT_ASSERT_MSG(finalized_,
                  "Cannot access data in a non-finalized AeroState!");
  return gas_mole_fractions_;
}

const AeroState::ColumnSpeciesView& AeroState::gas_mole_fractions() const {
  EKAT_ASSERT_MSG(finalized_,
                  "Cannot access data in a non-finalized AeroState!");
  return gas_mole_fractions_;
}

AeroState::ColumnView& AeroState::modal_num_density(int mode_index) {
  EKAT_ASSERT_MSG(finalized_,
                  "Cannot access data in a non-finalized AeroState!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < modal_num_densities_.size());
  return modal_num_densities_[mode_index];
}

const AeroState::ColumnView& AeroState::modal_num_density(int mode_index) const {
  EKAT_ASSERT_MSG(finalized_,
                  "Cannot access data in a non-finalized AeroState!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < modal_num_densities_.size());
  return modal_num_densities_[mode_index];
}

}

