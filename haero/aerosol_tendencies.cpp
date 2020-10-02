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

ColumnSpeciesView& AerosolTendencies::interstitial_aerosols(int mode_index) {
}

const ColumnSpeciesView& AerosolTendencies::interstitial_aerosols(int mode_index) const {
}

ColumnSpeciesView& AerosolTendencies::cloudborne_aerosols(int mode_index) {
}

const ColumnSpeciesView& AerosolTendencies::cloudborne_aerosols(int mode_index) const {
}

ColumnSpeciesView& AerosolTendencies::gas_mole_fractions() {
}

const ColumnSpeciesView& AerosolTendencies::gas_mole_fractions() const {
}

ColumnView& AerosolTendencies::modal_densities() {
}

const ColumnView& AerosolTendencies::modal_densities() const {
}

AerosolTendencies& AerosolTendencies::scale(Real factor) {
  return *this;
}

void AerosolTendencies::add(const AerosolTendencies& tendencies) {
}

}

