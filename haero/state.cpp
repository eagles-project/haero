#include "haero/state.hpp"

namespace haero {

State::State(int num_modes,
             const std::vector<int>& num_modal_species,
             int num_gas_species,
             int num_columns,
             int num_levels):
  num_modes_(num_modes),
  modal_species_offsets_(num_modal_species.size()),
  num_gas_species_(num_gas_species),
  num_columns_(num_columns),
  num_levels_(num_levels),
  interstitial_aerosols_(num_modes),
  cloud_borne_aerosols_(num_modes),
  gases_(num_gas_species),
  modal_n_(num_modes) {
}

State::~State() {
}

}

