#include "haero/processes/mam_rename_process.hpp"

#include <Kokkos_Pair.hpp>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>

#include "ekat_pack_math.hpp"
#include "haero/constants.hpp"
#include "haero/mode.hpp"

namespace haero {

namespace {

using HostExec = Kokkos::HostSpace;

KOKKOS_INLINE_FUNCTION
static Real pi_sixth() { return haero::Constants::pi_sixth; }

KOKKOS_INLINE_FUNCTION
static Real frelax() { return 27.0; }

KOKKOS_INLINE_FUNCTION
static Real sqrt_half() { return sqrt((Real)0.5); }

KOKKOS_INLINE_FUNCTION
static Real smallest_dryvol_value() { return 1.0e-25; }

// This will be replaced by another method when the rest of the original
// fortran has been ported. This is why _model_ is an unused parameter - it
// will be used to calculate _dest_mode_of_mode_mapping_ when that calculation
// is in place.
//
// \post _dest_mode_of_mode_mapping_ will have predetermined size of 4 due to
// temporary workaround.
void initialize_dest_mode_of_mode_mapping(
    view_1d_int_type& dest_mode_of_mode_mapping,
    const ModalAerosolConfig& model) {
  (void)model;
  Kokkos::resize(dest_mode_of_mode_mapping, 4);
  Kokkos::parallel_for(
      "delete", 4, KOKKOS_LAMBDA(const int i) {
        if (i == 1)
          dest_mode_of_mode_mapping[i] = 1;
        else
          dest_mode_of_mode_mapping[i] = 0;
      });
}

template <typename PackType>
static inline PackType compute_relaxed_volume_to_num_ratio(
    const Mode& mode, const PackType& diameter_for_current_mode) {
  const auto vol =
      mode.mean_particle_volume_from_diameter(diameter_for_current_mode);
  return 1. / vol;
}

}  // namespace

void MAMRenameProcess::init_(const ModalAerosolConfig& config) {
  num_modes = config.num_modes();
  max_aer = num_modes;
  num_populations = config.num_aerosol_populations;

  // Reserve memory for private fields
  Kokkos::resize(stencil, num_modes);
  Kokkos::resize(dgnumlo, num_modes);
  Kokkos::resize(dgnumhi, num_modes);
  Kokkos::resize(dgnum, num_modes);
  Kokkos::resize(alnsg, num_modes);
  Kokkos::resize(mass_2_vol, num_modes);
  Kokkos::resize(population_offsets, num_modes);

  Kokkos::resize(qi_del_growth, max_aer, num_modes);
  Kokkos::resize(qcld_del_growth, max_aer, num_modes);
  Kokkos::resize(dryvol_a, num_modes);
  Kokkos::resize(dryvol_c, num_modes);
  Kokkos::resize(deldryvol_a, num_modes);
  Kokkos::resize(deldryvol_c, num_modes);

  Kokkos::resize(start_species_for_mode, num_modes);
  Kokkos::resize(end_species_for_mode, num_modes);

  {
    auto population_offsets_host =
        Kokkos::create_mirror_view(population_offsets);
    for (int i = 0; i < num_modes; i++) {
      population_offsets_host(i) = config.population_index(i, 0);
    }
    Kokkos::deep_copy(population_offsets, population_offsets_host);
    Kokkos::deep_copy(start_species_for_mode, population_offsets_host);
  }

  {
    auto mass_2_vol_host = Kokkos::create_mirror_view(mass_2_vol);
    for (int i = 0; i < num_modes; i++) {
      mass_2_vol_host(i) = 0.0;
    }
    Kokkos::deep_copy(mass_2_vol, mass_2_vol_host);
  }

  initialize_dest_mode_of_mode_mapping(dest_mode_of_mode_mapping, config);

  view_1d_pack_type size_factor("size factor", num_modes);
  view_1d_pack_type fmode_dist_tail_fac("fmode dist tail fac", num_modes);
  view_1d_pack_type volume2num_lo_relaxed("volume to num low relaxed",
                                          num_modes);
  view_1d_pack_type volume2num_hi_relaxed("volume to num high relaxed",
                                          num_modes);
  view_1d_pack_type ln_diameter_tail_fac("ln of diameter tail fac", num_modes);
  view_1d_pack_type diameter_cutoff("diameter cutoff", num_modes);
  view_1d_pack_type diameter_belowcutoff("diameter below cutoff", num_modes);
  view_1d_pack_type ln_dia_cutoff("ln of diameter cutoff", num_modes);
  view_1d_pack_type dryvol_smallest("dry volume smallest", num_modes);

  std::size_t num_pairs = 0;

  find_renaming_pairs_(config, dest_mode_of_mode_mapping, num_pairs,
                       size_factor, fmode_dist_tail_fac, volume2num_lo_relaxed,
                       volume2num_hi_relaxed, ln_diameter_tail_fac,
                       diameter_cutoff, ln_dia_cutoff, diameter_belowcutoff,
                       dryvol_smallest);
}

template <typename T>
KOKKOS_INLINE_FUNCTION static auto compute_size_factor(T alnsg) {
  return pi_sixth() * exp(4.5 * pow(alnsg, 2));
}

void MAMRenameProcess::find_renaming_pairs_(
    const ModalAerosolConfig& config,
    view_1d_int_type dest_mode_of_mode_mapping, std::size_t& num_pairs,
    view_1d_pack_type size_factor, view_1d_pack_type fmode_dist_tail_fac,
    view_1d_pack_type volume2num_lo_relaxed,
    view_1d_pack_type volume2num_hi_relaxed,
    view_1d_pack_type ln_diameter_tail_fac, view_1d_pack_type diameter_cutoff,
    view_1d_pack_type ln_dia_cutoff, view_1d_pack_type diameter_belowcutoff,
    view_1d_pack_type dryvol_smallest) const {
  const auto N = dest_mode_of_mode_mapping.extent(0);

  // number of pairs allowed to do inter-mode particle transfer
  // (e.g. if we have a pair "mode_1<-->mode_2", mode_1 and mode_2 can
  // participate in inter-mode aerosol particle transfer where like particles in
  // mode_1 can be transferred to mode_2 and vice-versa)
  {
    auto mapping_host = Kokkos::create_mirror_view(dest_mode_of_mode_mapping);
    Kokkos::deep_copy(mapping_host, dest_mode_of_mode_mapping);
    device_kernels::reduce(
        N,
        KOKKOS_LAMBDA(const int i, std::size_t& acc) {
          acc += static_cast<std::size_t>(mapping_host(i) > 0);
        },
        num_pairs);
  }

  if (num_pairs == 0) return;

  auto dest_mode_of_mode_mapping_h =
      Kokkos::create_mirror_view(dest_mode_of_mode_mapping);
  auto dgnum_h = Kokkos::create_mirror_view(dgnum);
  auto dgnumhi_h = Kokkos::create_mirror_view(dgnumhi);
  auto alnsg_h = Kokkos::create_mirror_view(alnsg);
  auto size_factor_h = Kokkos::create_mirror_view(size_factor);
  auto fmode_dist_tail_fac_h = Kokkos::create_mirror_view(fmode_dist_tail_fac);
  auto volume2num_lo_relaxed_h =
      Kokkos::create_mirror_view(volume2num_lo_relaxed);
  auto volume2num_hi_relaxed_h =
      Kokkos::create_mirror_view(volume2num_hi_relaxed);
  auto ln_diameter_tail_fac_h =
      Kokkos::create_mirror_view(ln_diameter_tail_fac);
  auto diameter_cutoff_h = Kokkos::create_mirror_view(diameter_cutoff);
  auto ln_dia_cutoff_h = Kokkos::create_mirror_view(ln_dia_cutoff);
  auto diameter_belowcutoff_h =
      Kokkos::create_mirror_view(diameter_belowcutoff);
  auto dryvol_smallest_h = Kokkos::create_mirror_view(dryvol_smallest);

  Kokkos::deep_copy(dest_mode_of_mode_mapping_h, dest_mode_of_mode_mapping);
  Kokkos::deep_copy(alnsg_h, alnsg);
  Kokkos::deep_copy(dgnum_h, dgnum);
  Kokkos::deep_copy(dgnumhi_h, dgnumhi);
  Kokkos::deep_copy(size_factor_h, size_factor);
  Kokkos::deep_copy(fmode_dist_tail_fac_h, fmode_dist_tail_fac);
  Kokkos::deep_copy(volume2num_lo_relaxed_h, volume2num_lo_relaxed);
  Kokkos::deep_copy(volume2num_hi_relaxed_h, volume2num_hi_relaxed);
  Kokkos::deep_copy(ln_diameter_tail_fac_h, ln_diameter_tail_fac);
  Kokkos::deep_copy(diameter_cutoff_h, diameter_cutoff);
  Kokkos::deep_copy(ln_dia_cutoff_h, ln_dia_cutoff);
  Kokkos::deep_copy(diameter_belowcutoff_h, diameter_belowcutoff);
  Kokkos::deep_copy(dryvol_smallest_h, dryvol_smallest);

  for (int src_mode = 0; src_mode < N; src_mode++) {
    const auto& dest_mode_of_current_mode =
        dest_mode_of_mode_mapping_h(src_mode);

    const auto& alnsg_for_current_mode = alnsg_h[src_mode];
    const auto& alnsg_for_dest_mode = alnsg_h[dest_mode_of_current_mode];

    // Assign size factors for source and destination nodes.
    // Perform calculation only once, assignment twice.
    {
      const auto current_size_factor =
          compute_size_factor(alnsg_for_current_mode);
      size_factor_h[src_mode] = current_size_factor;
      size_factor_h[dest_mode_of_current_mode] = current_size_factor;
    }

    fmode_dist_tail_fac_h[src_mode] = sqrt_half() / alnsg_for_current_mode;
    dryvol_smallest_h[src_mode] = smallest_dryvol_value();

    // Set relaxed limits of ratios for current source/dest mode pair
    {
      const Mode& mode = config.aerosol_modes[src_mode];
      volume2num_lo_relaxed_h[src_mode] =
          compute_relaxed_volume_to_num_ratio(mode, dgnumhi_h[src_mode]);
      volume2num_lo_relaxed_h[dest_mode_of_current_mode] =
          compute_relaxed_volume_to_num_ratio(
              mode, dgnumhi_h[dest_mode_of_current_mode]);
    }

    // A factor for computing diameter at the tails of the distribution
    ln_diameter_tail_fac_h[src_mode] = 3.0 * pow(alnsg_for_current_mode, 2);

    // Cut-off (based on geometric mean) for making decision to do
    // inter-mode transfers
    //
    // TODO: use dummy values for _dgnum_, or assign to dgnum_low for the
    // moment. Have to figure out how to compute this. We will extract from
    // model at some point.
    {
      const auto sqrt_param_a =
          dgnum_h[src_mode] * exp(1.5 * pow(alnsg_for_current_mode, 2));

      const auto sqrt_param_b = dgnum_h[dest_mode_of_current_mode] *
                                exp(1.5 * pow(alnsg_for_dest_mode, 2));

      const auto diameter_cutoff_for_src_mode =
          sqrt(sqrt_param_a * sqrt_param_b);

      diameter_cutoff_h[src_mode] = diameter_cutoff_for_src_mode;

      ln_dia_cutoff_h[src_mode] = log(diameter_cutoff_for_src_mode);

      diameter_belowcutoff_h[src_mode] = 0.99 * diameter_cutoff_for_src_mode;
    }
  }

  Kokkos::deep_copy(dest_mode_of_mode_mapping, dest_mode_of_mode_mapping_h);
  Kokkos::deep_copy(alnsg, alnsg_h);
  Kokkos::deep_copy(dgnum, dgnum_h);
  Kokkos::deep_copy(dgnumhi, dgnumhi_h);
  Kokkos::deep_copy(size_factor, size_factor_h);
  Kokkos::deep_copy(fmode_dist_tail_fac, fmode_dist_tail_fac_h);
  Kokkos::deep_copy(volume2num_lo_relaxed, volume2num_lo_relaxed_h);
  Kokkos::deep_copy(volume2num_hi_relaxed, volume2num_hi_relaxed_h);
  Kokkos::deep_copy(ln_diameter_tail_fac, ln_diameter_tail_fac_h);
  Kokkos::deep_copy(diameter_cutoff, diameter_cutoff_h);
  Kokkos::deep_copy(ln_dia_cutoff, ln_dia_cutoff_h);
  Kokkos::deep_copy(diameter_belowcutoff, diameter_belowcutoff_h);
  Kokkos::deep_copy(dryvol_smallest, dryvol_smallest_h);
}

}  // namespace haero
