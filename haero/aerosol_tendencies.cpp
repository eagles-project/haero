#include "haero/aerosol_tendencies.hpp"

namespace haero {

AerosolTendencies::AerosolTendencies(const AerosolState& state) {
}

AerosolTendencies::~AerosolTendencies() {
}

int AerosolTendencies::num_aerosol_modes() const {
  return num_modes_;
}

int AerosolTendencies::num_aerosol_species(int mode_index) const {
}

int AerosolTendencies::num_gas_species() const {
}

AerosolTendencies::ColumnSpeciesView&
AerosolTendencies::interstitial_aerosols(int mode_index) {
}

const AerosolTendencies::ColumnSpeciesView&
AerosolTendencies::interstitial_aerosols(int mode_index) const {
}

AerosolTendencies::ColumnSpeciesView& AerosolTendencies::cloudborne_aerosols(int mode_index) {
}

const AerosolTendencies::ColumnSpeciesView& AerosolTendencies::cloudborne_aerosols(int mode_index) const {
}

AerosolTendencies::ColumnSpeciesView& AerosolTendencies::gas_mole_fractions() {
}

const AerosolTendencies::ColumnSpeciesView& AerosolTendencies::gas_mole_fractions() const {
}

AerosolTendencies::ColumnView& AerosolTendencies::modal_densities() {
}

const AerosolTendencies::ColumnView& AerosolTendencies::modal_densities() const {
}

AerosolTendencies& AerosolTendencies::scale(Real factor) {
  return *this;
}

void AerosolTendencies::add(const AerosolTendencies& tendencies) {
}

}

