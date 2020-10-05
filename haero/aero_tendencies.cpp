#include "haero/aero_tendencies.hpp"

namespace haero {

AeroTendencies::AeroTendencies(const AeroState& state) {
  EKAT_ASSERT_MSG(state.is_finalized(),
                  "Cannot construct AeroTendencies from an AeroState that isn't finalized!");
  int num_modes = state.num_aerosol_modes();
  int num_columns = state.num_columns();
  int num_levels = state.num_levels();
  for (int m = 0; m < num_modes; ++m) {
    int num_aero_species = state.num_aerosol_species(m);
    int_aero_species_.push_back(ColumnSpeciesView("dqi/dt",
                                                  num_columns,
                                                  num_aero_species,
                                                  num_levels));
    cld_aero_species_.push_back(ColumnSpeciesView("dqc/dt",
                                                  num_columns,
                                                  num_aero_species,
                                                  num_levels));
    modal_num_densities_.push_back(ColumnView("dn/dt", num_columns, num_levels));
  }
  int num_gas_species = state.num_gas_species();
  gas_mole_fractions_ = ColumnSpeciesView("d(gas_mole_fractions)/dt",
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
  return gas_mole_fractions_.extent(2);
}

int AeroTendencies::num_columns() const {
  return gas_mole_fractions_.extent(1);
}

int AeroTendencies::num_levels() const {
  return gas_mole_fractions_.extent(3);
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

