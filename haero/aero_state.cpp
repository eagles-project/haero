#include "haero/aero_state.hpp"

namespace haero {

AeroState::AeroState(int num_columns, int num_levels):
  num_columns_(num_columns),
  num_levels_(num_levels),
  assembled_(false) {
}

AeroState::~AeroState() {
}

int AeroState::add_aerosol_mode(const Mode& mode,
                                const std::vector<Species>& aero_species) {
  auto int_view_name = std::string("interstitial_aerosols(") + mode.name +
                       std::string(")");
  auto int_aero_data = ManagedColumnSpeciesView(int_view_name,
                                                num_columns_,
                                                aero_species.size(),
                                                num_levels_);
  managed_column_species_views_.push_back(int_aero_data);
  auto cld_view_name = std::string("cloudborne_aerosols(") + mode.name +
                       std::string(")");
  auto cld_aero_data = ManagedColumnSpeciesView(cld_view_name,
                                                num_columns_,
                                                aero_species.size(),
                                                num_levels_);
  auto modal_data_name = std::string("num_density(") + mode.name +
                         std::string(")");
  auto modal_data = ManagedColumnView(modal_data_name,
                                      num_columns_,
                                      num_levels_);
  managed_column_species_views_.push_back(cld_aero_data);
  return add_aerosol_mode(mode,
                          aero_species,
                          ColumnSpeciesView(int_aero_data),
                          ColumnSpeciesView(cld_aero_data),
                          ColumnView(modal_data));
}

int AeroState::add_aerosol_mode(const Mode& mode,
                                const std::vector<Species>& aero_species,
                                ColumnSpeciesView int_aero_data,
                                ColumnSpeciesView cld_aero_data,
                                ColumnView modal_data) {
  EKAT_ASSERT_MSG(not assembled_,
                  "Cannot add an aerosol mode to an assembled AeroState!");
  aero_species_names_.push_back(std::vector<std::string>());
  for (int s = 0; s < aero_species.size(); ++s) {
    aero_species_names_.back().push_back(aero_species[s].symbol);
  }
  int_aero_species_.push_back(int_aero_data);
  cld_aero_species_.push_back(cld_aero_data);
  return static_cast<int>(aero_species_names_.size()-1);
}

void AeroState::add_gas_species(const std::vector<Species>& gas_species) {
  int num_species = gas_species.size();
  auto gas_mole_fracs = ManagedColumnSpeciesView("gas_mole_fractions",
                                                 num_columns_,
                                                 num_species,
                                                 num_levels_);
  managed_column_species_views_.push_back(gas_mole_fracs);
  add_gas_species(gas_species, ColumnSpeciesView(gas_mole_fracs));
}

void AeroState::add_gas_species(const std::vector<Species>& gas_species,
                                ColumnSpeciesView gas_data) {
  EKAT_ASSERT_MSG(not assembled_,
                  "Cannot add gas species to an assembled AeroState!");
  EKAT_ASSERT_MSG(gas_species_names_.empty(),
                  "Cannot add more than one set of gas species to an AeroState!");
  for (int s = 0; s < gas_species.size(); ++s) {
    gas_species_names_.push_back(gas_species[s].symbol);
  }
  gas_mole_fractions_ = gas_data;
}

void AeroState::assemble() {
  assembled_ = true;
}

bool AeroState::is_assembled() const {
  return assembled_;
}

int AeroState::num_aerosol_modes() const {
  return aero_species_names_.size();
}

int AeroState::num_aerosol_species(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < aero_species_names_.size());
  return aero_species_names_[mode_index].size();
}

int AeroState::num_gas_species() const {
  return static_cast<int>(gas_species_names_.size());
}

int AeroState::num_columns() const {
  return num_columns_;
}

int AeroState::num_levels() const {
  return num_levels_;
}

AeroState::ColumnSpeciesView&
AeroState::interstitial_aerosols(int mode_index) {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in a unassembled AeroState!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < aero_species_names_.size());
  return int_aero_species_[mode_index];
}

const AeroState::ColumnSpeciesView&
AeroState::interstitial_aerosols(int mode_index) const {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled AeroState!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < aero_species_names_.size());
  return int_aero_species_[mode_index];
}

AeroState::ColumnSpeciesView& AeroState::cloudborne_aerosols(int mode_index) {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled AeroState!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < aero_species_names_.size());
  return cld_aero_species_[mode_index];
}

const AeroState::ColumnSpeciesView&
AeroState::cloudborne_aerosols(int mode_index) const {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled AeroState!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < aero_species_names_.size());
  return cld_aero_species_[mode_index];
}

AeroState::ColumnSpeciesView& AeroState::gas_mole_fractions() {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled AeroState!");
  return gas_mole_fractions_;
}

const AeroState::ColumnSpeciesView& AeroState::gas_mole_fractions() const {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled AeroState!");
  return gas_mole_fractions_;
}

AeroState::ColumnView& AeroState::modal_num_density(int mode_index) {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled AeroState!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < modal_num_densities_.size());
  return modal_num_densities_[mode_index];
}

const AeroState::ColumnView& AeroState::modal_num_density(int mode_index) const {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled AeroState!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < modal_num_densities_.size());
  return modal_num_densities_[mode_index];
}

}

