#include "haero/aero_tendencies.hpp"

namespace haero {

AeroTendencies::AeroTendencies(const AeroState& state)
{
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

