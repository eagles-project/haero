#include "haero/processes/mam_calcsize_process.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>

namespace haero {

MAMCalcsizeProcess::MAMCalcsizeProcess()
    : DeviceAerosolProcess<MAMCalcsizeProcess>("MAMCalcsizeProcess") {}
//------------------------------------------------------------------------
//                                Accessors
//------------------------------------------------------------------------

void MAMCalcsizeProcess::init_(const ModalAerosolConfig &modal_aerosol_config) {
  nmodes = modal_aerosol_config.num_aerosol_modes();

  num_populations = modal_aerosol_config.num_aerosol_populations;

  Kokkos::resize(num_mode_species, nmodes);
  {
    auto nm_spec = Kokkos::create_mirror_view(num_mode_species);
    for (int i = 0; i < nmodes; i++)
      nm_spec[i] = modal_aerosol_config.aerosol_species_for_mode(i).size();

    max_nspec = nm_spec(0);
    Kokkos::resize(density, max_nspec);

    for (int i = 0; i < nm_spec.extent(0); i++)
      max_nspec = std::max(max_nspec, static_cast<std::size_t>(nm_spec(i)));

    Kokkos::deep_copy(num_mode_species, nm_spec);
  }

  Kokkos::resize(spec_density, max_nspec * nmodes);
  {
    const auto h_spec_density = Kokkos::create_mirror_view(spec_density);
    const auto all_species = modal_aerosol_config.aerosol_species;
    for (int i = 0; i < all_species.size(); i++)
      h_spec_density(i) = all_species[i].density;
    Kokkos::deep_copy(spec_density, h_spec_density);
  }

  Kokkos::resize(population_offsets, nmodes);
  {
    const auto pop_offsets = Kokkos::create_mirror_view(population_offsets);
    for (int i = 0; i < nmodes; i++)
      pop_offsets[i] = modal_aerosol_config.population_index(i, 0);
    Kokkos::deep_copy(population_offsets, pop_offsets);
  }

  Kokkos::resize(do_aitacc_transfer_allowed, npair_csizxf);
  {
    const auto h_do_xfer_allowed =
        Kokkos::create_mirror_view(do_aitacc_transfer_allowed);
    for (std::size_t i = 0; i < npair_csizxf; i++)
      h_do_xfer_allowed(i) = 0 /*false*/;
    Kokkos::deep_copy(do_aitacc_transfer_allowed, h_do_xfer_allowed);
  }

  Kokkos::resize(v2nmin_nmodes, nmodes);
  Kokkos::resize(v2nmax_nmodes, nmodes);
  Kokkos::resize(v2nnom_nmodes, nmodes);
  Kokkos::resize(dgnmin_nmodes, nmodes);
  Kokkos::resize(dgnmax_nmodes, nmodes);
  Kokkos::resize(dgnnom_nmodes, nmodes);
  Kokkos::resize(common_factor_nmodes, nmodes);

  Kokkos::resize(srcmode_csizxf, maxpair_csizxf);
  Kokkos::resize(destmode_csizxf, maxpair_csizxf);
  {
    const auto h_srcmode = Kokkos::create_mirror_view(srcmode_csizxf);
    const auto h_destmode = Kokkos::create_mirror_view(destmode_csizxf);
    for (std::size_t i = 0; i < maxpair_csizxf; i++) {
      h_srcmode(i) = -1;
      h_destmode(i) = -1;
    }
    Kokkos::deep_copy(srcmode_csizxf, h_srcmode);
    Kokkos::deep_copy(destmode_csizxf, h_destmode);
  }

  Kokkos::resize(dgncur_a, num_levels_upper_bound, nmodes);
  Kokkos::resize(dgncur_c, num_levels_upper_bound, nmodes);
  Kokkos::resize(v2ncur_a, num_levels_upper_bound, nmodes);
  Kokkos::resize(v2ncur_c, num_levels_upper_bound, nmodes);

  Kokkos::resize(dryvol_a, num_levels_upper_bound);
  Kokkos::resize(dryvol_c, num_levels_upper_bound);

  Kokkos::resize(ait_mode_inter, maxpair_csizxf);
  Kokkos::resize(acc_mode_inter, maxpair_csizxf);
  Kokkos::resize(ait_mode_cldbrn, maxpair_csizxf);
  Kokkos::resize(acc_mode_cldbrn, maxpair_csizxf);
  Kokkos::resize(nspec_common, maxpair_csizxf);

  Kokkos::resize(ait_spec_in_acc_inter, maxpair_csizxf, index_mapping_extent_1);
  Kokkos::resize(acc_spec_in_ait_inter, maxpair_csizxf, index_mapping_extent_1);
  Kokkos::resize(ait_spec_in_acc_cldbrn, maxpair_csizxf,
                 index_mapping_extent_1);
  Kokkos::resize(acc_spec_in_ait_cldbrn, maxpair_csizxf,
                 index_mapping_extent_1);

  {
    const auto h_v2nmin_nmodes = Kokkos::create_mirror_view(v2nmin_nmodes);
    const auto h_v2nmax_nmodes = Kokkos::create_mirror_view(v2nmax_nmodes);
    const auto h_v2nnom_nmodes = Kokkos::create_mirror_view(v2nnom_nmodes);
    const auto h_dgnmin_nmodes = Kokkos::create_mirror_view(dgnmin_nmodes);
    const auto h_dgnmax_nmodes = Kokkos::create_mirror_view(dgnmax_nmodes);
    const auto h_dgnnom_nmodes = Kokkos::create_mirror_view(dgnnom_nmodes);
    const auto h_common_factor_nmodes =
        Kokkos::create_mirror_view(common_factor_nmodes);

    for (int i = 0; i < nmodes; i++) {
      const auto &mode = modal_aerosol_config.aerosol_modes[i];
      using T = decltype(v2nmin_nmodes)::value_type;
      h_v2nmin_nmodes[i] = mode.min_vol_to_num_ratio<T>();
      h_v2nmax_nmodes[i] = mode.max_vol_to_num_ratio<T>();
      h_v2nnom_nmodes[i] = mode.nom_vol_to_num_ratio<T>();
      h_dgnmin_nmodes[i] = mode.min_diameter;
      h_dgnmax_nmodes[i] = mode.max_diameter;
      h_dgnnom_nmodes[i] = mode.nom_diameter;
      h_common_factor_nmodes[i] =
          std::exp(4.5 * std::log(std::pow(mode.mean_std_dev, 2.0))) *
          Constants::pi_sixth;
    }

    Kokkos::deep_copy(v2nmin_nmodes, h_v2nmin_nmodes);
    Kokkos::deep_copy(v2nmax_nmodes, h_v2nmax_nmodes);
    Kokkos::deep_copy(v2nnom_nmodes, h_v2nnom_nmodes);
    Kokkos::deep_copy(dgnmin_nmodes, h_dgnmin_nmodes);
    Kokkos::deep_copy(dgnmax_nmodes, h_dgnmax_nmodes);
    Kokkos::deep_copy(dgnnom_nmodes, h_dgnnom_nmodes);
    Kokkos::deep_copy(common_factor_nmodes, h_common_factor_nmodes);
  }

  Kokkos::resize(drv_a_aitsv, num_levels_upper_bound);
  Kokkos::resize(drv_a_accsv, num_levels_upper_bound);
  Kokkos::resize(drv_c_aitsv, num_levels_upper_bound);
  Kokkos::resize(drv_c_accsv, num_levels_upper_bound);
  Kokkos::resize(num_a_aitsv, num_levels_upper_bound);
  Kokkos::resize(num_a_accsv, num_levels_upper_bound);
  Kokkos::resize(num_c_aitsv, num_levels_upper_bound);
  Kokkos::resize(num_c_accsv, num_levels_upper_bound);
  Kokkos::resize(drv_a_sv, num_levels_upper_bound, nmodes);
  Kokkos::resize(drv_c_sv, num_levels_upper_bound, nmodes);
  Kokkos::resize(num_a_sv, num_levels_upper_bound, nmodes);
  Kokkos::resize(num_c_sv, num_levels_upper_bound, nmodes);

  Kokkos::resize(no_transfer_acc2ait, max_nspec);

  aitken_idx = modal_aerosol_config.aerosol_mode_index("aitken");
  accum_idx = modal_aerosol_config.aerosol_mode_index("accum");

  find_species_mapping(modal_aerosol_config);
}

void MAMCalcsizeProcess::find_species_mapping(
    const ModalAerosolConfig &modal_aerosol_config) {
  static constexpr auto do_adjust_allowed = true;

  // ------------------------------------------------------------------------------------------
  // Original comment in the code:
  // "do_aitacc_transfer_allowed" allows aitken <--> accum mode transfer to be
  // turned on/off NOTE: it can only be true when aitken & accum modes are both
  // present
  //       and have prognosed number and diagnosed surface/sigmag
  // ------------------------------------------------------------------------------------------

  bool any_xfer_allowed = false;

  // find out mode index for accum and aitken modes in the prognostic
  // radiation list (rad_climate)
  const auto nait = aitken_idx;  // mode number of aitken mode
  const auto nacc = accum_idx;   // mode number of accumulation mode

  for (std::size_t i = 0; i < npair_csizxf; i++) {
    // find out accumulation and aitken modes in the radiation list
    // FIXME: This should get aitken and accumulation mode indices for
    // diagnostic lists
    const auto iacc = nacc;  // FIXME:Hardwired!! we should get this index
                             // frm the diagnostic lists
    const auto iait = nait;  // FIXME:Hardwired!! we should get this index
                             // frm the diagnostic lists

    // find out if aitken or accumulation modes exist in the radiation list
    // (a positive value means that the mode exists)
    const auto accum_exists = iacc > 0;
    const auto aitken_exists = iait > 0;

    // if both aitken and accumulation modes exist, make it True and assign
    // source and destination mode variables
    const auto assign_mapping =
        accum_exists and aitken_exists and (iacc != iait);

    if (assign_mapping) {
      any_xfer_allowed = true;
      do_aitacc_transfer_allowed(i) = 1 /*true*/;
      srcmode_csizxf(i) = iait;
      destmode_csizxf(i) = iacc;
    }
  }

  if (any_xfer_allowed) {
    for (std::size_t ilist = 0; ilist < npair_csizxf; ilist++) {
      if (do_aitacc_transfer_allowed(ilist)) {
        // aitken  mode of this list
        const auto imode_ait = srcmode_csizxf(ilist);
        // accumulation mode of this list
        const auto imode_acc = destmode_csizxf(ilist);

        //----------------------------------------------------------------------------------------
        // Aerosol *number* indices mapping between aitken and accumulation
        // modes in this list
        //----------------------------------------------------------------------------------------

        ait_mode_inter(ilist) = imode_ait;  // aitken mode
        acc_mode_inter(ilist) = imode_acc;  // accumulation mode

        // indices for cloudborne species array are same as interstitial species

        // aitken mode (cloud borne)
        ait_mode_cldbrn(ilist) = imode_ait;

        // accumulation mode (cloud borne)
        acc_mode_cldbrn(ilist) = imode_acc;

        //--------------------------------------------------------------------------------------
        // find aerosol *mass* indices mapping between aitken and accumulation
        // modes in this list
        //--------------------------------------------------------------------------------------

        // find number of species in the aitken mode of this ilist
        const auto nspec_ait =
            modal_aerosol_config.aerosol_species_for_mode(imode_ait).size();

        // find number of species in the accumulation mode of this list
        const auto nspec_acc =
            modal_aerosol_config.aerosol_species_for_mode(imode_acc).size();

        for (std::size_t i = 0; i < maxpair_csizxf; i++)
          for (std::size_t j = 0; j < index_mapping_extent_1; j++) {
            ait_spec_in_acc_inter(i, j) = -1;
            acc_spec_in_ait_inter(i, j) = -1;
            ait_spec_in_acc_cldbrn(i, j) = -1;
            acc_spec_in_ait_cldbrn(i, j) = -1;
          }

        int icnt = 0;
        for (std::size_t ispec_ait = 1; ispec_ait <= nspec_ait; ispec_ait++) {
          const auto spec_name_ait =
              modal_aerosol_config
                  .aerosol_species_for_mode(imode_ait)[ispec_ait]
                  .symbol();
          // Now find this specie index in the destination mode
          for (std::size_t ispec_acc = 1; ispec_acc <= nspec_acc; ispec_acc++) {
            // find if specie in acc mode is same as ait or not
            const auto spec_name_acc =
                modal_aerosol_config
                    .aerosol_species_for_mode(imode_acc)[ispec_acc]
                    .symbol();
            if (spec_name_ait == spec_name_acc) {
              // if there is a match, find indices of species in cnst array
              icnt++;
              ait_spec_in_acc_inter(ilist, icnt) =
                  modal_aerosol_config.population_index(nait, ispec_ait);
              acc_spec_in_ait_inter(ilist, icnt) =
                  modal_aerosol_config.population_index(nacc, ispec_acc);
              ait_spec_in_acc_cldbrn(ilist, icnt) =
                  modal_aerosol_config.population_index(nait, ispec_ait);
              acc_spec_in_ait_cldbrn(ilist, icnt) =
                  modal_aerosol_config.population_index(nacc, ispec_acc);
            }
          }
        }
        nspec_common(ilist) = icnt;
      }
    }
  }
}

}  // namespace haero
