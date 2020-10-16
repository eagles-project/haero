#include "haero/prognostics.hpp"

namespace haero {

Prognostics::Prognostics(int num_columns, int num_levels):
  num_columns_(num_columns),
  num_levels_(num_levels),
  assembled_(false) {
}

Prognostics::~Prognostics() {
}

int Prognostics::add_aerosol_mode(const Mode& mode,
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
                          ColumnSpeciesView(cld_aero_data));
}

int Prognostics::add_aerosol_mode(const Mode& mode,
                                  const std::vector<Species>& aero_species,
                                  ColumnSpeciesView int_aero_data,
                                  ColumnSpeciesView cld_aero_data) {
  EKAT_ASSERT_MSG(not assembled_,
                  "Cannot add an aerosol mode to an assembled Prognostics!");
  aero_species_names_.push_back(std::vector<std::string>());
  for (int s = 0; s < aero_species.size(); ++s) {
    aero_species_names_.back().push_back(aero_species[s].symbol);
  }
  int_aero_species_.push_back(int_aero_data);
  cld_aero_species_.push_back(cld_aero_data);
  return static_cast<int>(aero_species_names_.size()-1);
}

void Prognostics::add_gas_species(const std::vector<Species>& gas_species) {
  int num_species = gas_species.size();
  auto gas_mole_fracs = ManagedColumnSpeciesView("gas_mole_fractions",
                                                 num_columns_,
                                                 num_species,
                                                 num_levels_);
  managed_column_species_views_.push_back(gas_mole_fracs);
  add_gas_species(gas_species, ColumnSpeciesView(gas_mole_fracs));
}

void Prognostics::add_gas_species(const std::vector<Species>& gas_species,
                                  ColumnSpeciesView gas_data) {
  EKAT_ASSERT_MSG(not assembled_,
                  "Cannot add gas species to an assembled Prognostics!");
  EKAT_ASSERT_MSG(gas_species_names_.empty(),
                  "Cannot add more than one set of gas species to an Prognostics!");
  for (int s = 0; s < gas_species.size(); ++s) {
    gas_species_names_.push_back(gas_species[s].symbol);
  }
  gas_mole_fractions_ = gas_data;
}

void Prognostics::set_modal_number_densities(ModalColumnView modal_num_densities) {
  EKAT_ASSERT_MSG(not assembled_,
                  "Cannot set modal number densities in an assembled Prognostics!");
  modal_num_densities_ = modal_num_densities;
}

void Prognostics::assemble() {
  assembled_ = true;
  // TODO: Allocate modal number densities, etc.
}

bool Prognostics::is_assembled() const {
  return assembled_;
}

int Prognostics::num_aerosol_modes() const {
  return aero_species_names_.size();
}

int Prognostics::num_aerosol_species(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < aero_species_names_.size());
  return aero_species_names_[mode_index].size();
}

int Prognostics::num_gas_species() const {
  return static_cast<int>(gas_species_names_.size());
}

int Prognostics::num_columns() const {
  return num_columns_;
}

int Prognostics::num_levels() const {
  return num_levels_;
}

Prognostics::ColumnSpeciesView&
Prognostics::interstitial_aerosols(int mode_index) {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in a unassembled Prognostics!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < aero_species_names_.size());
  return int_aero_species_[mode_index];
}

const Prognostics::ColumnSpeciesView&
Prognostics::interstitial_aerosols(int mode_index) const {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled Prognostics!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < aero_species_names_.size());
  return int_aero_species_[mode_index];
}

Prognostics::ColumnSpeciesView& Prognostics::cloudborne_aerosols(int mode_index) {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled Prognostics!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < aero_species_names_.size());
  return cld_aero_species_[mode_index];
}

const Prognostics::ColumnSpeciesView&
Prognostics::cloudborne_aerosols(int mode_index) const {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled Prognostics!");
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < aero_species_names_.size());
  return cld_aero_species_[mode_index];
}

Prognostics::ColumnSpeciesView& Prognostics::gas_mole_fractions() {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled Prognostics!");
  return gas_mole_fractions_;
}

const Prognostics::ColumnSpeciesView& Prognostics::gas_mole_fractions() const {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled Prognostics!");
  return gas_mole_fractions_;
}

Prognostics::ModalColumnView& Prognostics::modal_num_densities() {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled Prognostics!");
  return modal_num_densities_;
}

const Prognostics::ModalColumnView& Prognostics::modal_num_densities() const {
  EKAT_ASSERT_MSG(assembled_,
                  "Cannot access data in an unassembled Prognostics!");
  return modal_num_densities_;
}

void Prognostics::scale_and_add(Real scale_factor,
                                const Tendencies& tendencies) {
}

}

