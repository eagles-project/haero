#ifndef HAERO_CONVERSIONS_HPP
#define HAERO_CONVERSIONS_HPP

#include "haero/haero.hpp"
#include "haero/physical_constants.hpp"

namespace haero {

namespace conversions {

/// Given a number concentration for a species or mixture [m-3], computes and
/// returns a mass mixing ratio [kg species/kg air] based on its molecular
/// weight and on the density of dry air in the vicinity.
/// @param [in] number_conc The number concentration of the species/mixture
/// [m-3]
/// @param [in] molecular_wt The molecular weight of the species/mixture
/// [kmol/kg]
/// @param [in] air_density The mass density of dry air [kg/m3]
KOKKOS_INLINE_FUNCTION
PackType mmr_from_number_conc(const PackType& number_conc, Real molecular_wt,
                              Real air_density) {
  return number_conc * molecular_wt / (air_density * constants::avogadro);
}

/// Given a mass mixing ratio (mmr) for a species or mixture [kg species/kg
/// air], computes and returns a number density [m-3] based on its molecular
/// weight and on the density of dry air in the vicinity.
/// @param [in] mmr The number concentration of the species/mixture [m3]
/// @param [in] molecular_wt The molecular weight of the species/mixture
/// [kmol/kg]
/// @param [in] air_density The mass density of dry air [kg/m3]
KOKKOS_INLINE_FUNCTION
PackType number_conc_from_mmr(const PackType& mmr, Real molecular_wt,
                              Real air_density) {
  return mmr * (air_density * constants::avogadro) / molecular_wt;
}

/// Given a molar mixing ratio (vmr) for a species or mixture
/// [kmol species/kmol air], computes and returns a mass mixing ratio
/// [kg species/kg air] based on its molecular weight.
/// @param [in] vmr The molar mixing ratio of the species/mixture [kmol/kmol
/// air]
/// @param [in] molecular_wt The molecular weight of the species/mixture
/// [kmol/kg]
KOKKOS_INLINE_FUNCTION
PackType mmr_from_vmr(const PackType& vmr, Real molecular_wt) {
  return vmr * molecular_wt / constants::molec_weight_dry_air;
}

/// Given a mass mixing ratio (mmr) for a species or mixture [kg species/kg
/// air], computes and returns a molar mixing ratio [kmol species/k air] based
/// on its molecular weight.
/// @param [in] mmr The molar mixing ratio of the species/mixture [kmol/kmol
/// air]
/// @param [in] molecular_wt The molecular weight of the species/mixture
/// [kmol/kg]
KOKKOS_INLINE_FUNCTION
PackType vmr_from_mmr(const PackType& mmr, Real molecular_wt) {
  return mmr * constants::molec_weight_dry_air / molecular_wt;
}

}  // namespace conversions

}  // namespace haero

#endif
