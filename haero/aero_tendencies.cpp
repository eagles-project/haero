#include "haero/aero_tendencies.hpp"

namespace haero {

AeroTendencies::AeroTendencies(const AeroState& state) {
  EKAT_ASSERT_MSG(state.is_assembled(),
                  "Cannot construct AeroTendencies from an unassembled AeroState!");
  int num_modes = state.num_aerosol_modes();
  int num_columns = state.num_columns();
  int num_levels = state.num_levels();
  for (int m = 0; m < num_modes; ++m) {
    int num_aero_species = state.num_aerosol_species(m);
    auto int_view_name = std::string("d/dt[") +
                         state.interstitial_aerosols(m).label() +
                         std::string(")]");
    int_aero_species_.push_back(ColumnSpeciesView("dqi/dt",
                                                  num_columns,
                                                  num_aero_species,
                                                  num_levels));
    auto cld_view_name = std::string("d/dt[") +
                         state.cloudborne_aerosols(m).label() +
                         std::string(")]");
    cld_aero_species_.push_back(ColumnSpeciesView(cld_view_name,
                                                  num_columns,
                                                  num_aero_species,
                                                  num_levels));
    auto n_view_name = std::string("d/dt[") +
                       state.modal_num_density(m).label() +
                       std::string(")]");
    modal_num_densities_.push_back(ColumnView(n_view_name, num_columns, num_levels));
  }
  int num_gas_species = state.num_gas_species();
  auto gas_view_name = std::string("d/dt[") +
                       state.gas_mole_fractions().label() +
                       std::string(")]");
  gas_mole_fractions_ = ColumnSpeciesView(gas_view_name,
                                          num_columns,
                                          num_gas_species,
                                          num_levels);
}

AeroTendencies::~AeroTendencies() {
}

int AeroTendencies::num_aerosol_modes() const {
  return static_cast<int>(modal_num_densities_.size());
}

int AeroTendencies::num_aerosol_species(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < int_aero_species_.size());
  return int_aero_species_[mode_index].extent(2);
}

int AeroTendencies::num_gas_species() const {
  return gas_mole_fractions_.extent(1);
}

int AeroTendencies::num_columns() const {
  return gas_mole_fractions_.extent(0);
}

int AeroTendencies::num_levels() const {
  return gas_mole_fractions_.extent(2);
}

AeroTendencies::ColumnSpeciesView&
AeroTendencies::interstitial_aerosols(int mode_index) {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < int_aero_species_.size());
  return int_aero_species_[mode_index];
}

const AeroTendencies::ColumnSpeciesView&
AeroTendencies::interstitial_aerosols(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < int_aero_species_.size());
  return int_aero_species_[mode_index];
}

AeroTendencies::ColumnSpeciesView&
AeroTendencies::cloudborne_aerosols(int mode_index) {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < cld_aero_species_.size());
  return cld_aero_species_[mode_index];
}

const AeroTendencies::ColumnSpeciesView&
AeroTendencies::cloudborne_aerosols(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < cld_aero_species_.size());
  return cld_aero_species_[mode_index];
}

AeroTendencies::ColumnSpeciesView& AeroTendencies::gas_mole_fractions() {
  return gas_mole_fractions_;
}

const AeroTendencies::ColumnSpeciesView& AeroTendencies::gas_mole_fractions() const {
  return gas_mole_fractions_;
}

AeroTendencies::ColumnView&
AeroTendencies::modal_num_density(int mode_index) {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < modal_num_densities_.size());
  return modal_num_densities_[mode_index];
}

const AeroTendencies::ColumnView&
AeroTendencies::modal_num_density(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < modal_num_densities_.size());
  return modal_num_densities_[mode_index];
}

AeroTendencies& AeroTendencies::scale(Real factor) {
  return *this;
}

void AeroTendencies::add(const AeroTendencies& tendencies) {
}

}

