#include "haero/tendencies.hpp"

namespace haero {

Tendencies::Tendencies(const Prognostics& prognostics) {
  EKAT_ASSERT_MSG(prognostics.is_assembled(),
                  "Cannot construct Tendencies from unassembled Prognostics!");
  int num_modes = prognostics.num_aerosol_modes();
  int num_columns = prognostics.num_columns();
  int num_levels = prognostics.num_levels();
  for (int m = 0; m < num_modes; ++m) {
    int num_aero_species = prognostics.num_aerosol_species(m);
    auto int_view_name = std::string("d/dt[") +
                         prognostics.interstitial_aerosols(m).label() +
                         std::string(")]");
    int_aero_species_.push_back(ColumnSpeciesView("dqi/dt",
                                                  num_columns,
                                                  num_aero_species,
                                                  num_levels));
    auto cld_view_name = std::string("d/dt[") +
                         prognostics.cloudborne_aerosols(m).label() +
                         std::string(")]");
    cld_aero_species_.push_back(ColumnSpeciesView(cld_view_name,
                                                  num_columns,
                                                  num_aero_species,
                                                  num_levels));
  }
  auto n_view_name = std::string("d/dt[") +
                     prognostics.modal_num_densities().label() +
                     std::string(")]");
  modal_num_densities_ = ModalColumnView(n_view_name, num_modes, num_columns,
                                         num_levels);
  int num_gas_species = prognostics.num_gas_species();
  auto gas_view_name = std::string("d/dt[") +
                       prognostics.gas_mole_fractions().label() +
                       std::string(")]");
  gas_mole_fractions_ = ColumnSpeciesView(gas_view_name,
                                          num_columns,
                                          num_gas_species,
                                          num_levels);
}

Tendencies::~Tendencies() {
}

int Tendencies::num_aerosol_modes() const {
  return modal_num_densities_.extent(0);
}

int Tendencies::num_aerosol_species(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < int_aero_species_.size());
  return int_aero_species_[mode_index].extent(1);
}

int Tendencies::num_gas_species() const {
  return gas_mole_fractions_.extent(1);
}

int Tendencies::num_columns() const {
  return gas_mole_fractions_.extent(0);
}

int Tendencies::num_levels() const {
  return gas_mole_fractions_.extent(2);
}

Tendencies::ColumnSpeciesView&
Tendencies::interstitial_aerosols(int mode_index) {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < int_aero_species_.size());
  return int_aero_species_[mode_index];
}

const Tendencies::ColumnSpeciesView&
Tendencies::interstitial_aerosols(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < int_aero_species_.size());
  return int_aero_species_[mode_index];
}

Tendencies::ColumnSpeciesView&
Tendencies::cloudborne_aerosols(int mode_index) {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < cld_aero_species_.size());
  return cld_aero_species_[mode_index];
}

const Tendencies::ColumnSpeciesView&
Tendencies::cloudborne_aerosols(int mode_index) const {
  EKAT_ASSERT(mode_index >= 0);
  EKAT_ASSERT(mode_index < cld_aero_species_.size());
  return cld_aero_species_[mode_index];
}

Tendencies::ColumnSpeciesView& Tendencies::gas_mole_fractions() {
  return gas_mole_fractions_;
}

const Tendencies::ColumnSpeciesView& Tendencies::gas_mole_fractions() const {
  return gas_mole_fractions_;
}

Tendencies::ModalColumnView&
Tendencies::modal_num_densities() {
  return modal_num_densities_;
}

const Tendencies::ModalColumnView&
Tendencies::modal_num_densities() const {
  return modal_num_densities_;
}

Tendencies& Tendencies::scale(Real factor) {
  return *this;
}

void Tendencies::accumulate(const Tendencies& tendencies) {
}

}

