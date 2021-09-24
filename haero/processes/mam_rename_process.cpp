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

using size_type = MAMRenameProcess::size_type;
using integral_type = MAMRenameProcess::integral_type;

static constexpr auto& pi_sixth = haero::Constants::pi_sixth;
static constexpr Real frelax = 27.0;
static const Real sqrt_half = sqrt((Real)0.5);
static constexpr Real smallest_dryvol_value = 1.0e-25;

// This will be replaced by another method when the rest of the original
// fortran has been ported. This is why _model_ is an unused parameter - it
// will be used to calculate _dest_mode_of_mode_mapping_ when that calculation
// is in place.
//
// \post _dest_mode_of_mode_mapping_ will have predetermined size of 4 due to
// temporary workaround.
void initialize_dest_mode_of_mode_mapping(
    view_1d_int_type dest_mode_of_mode_mapping,
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

MAMRenameProcess::MAMRenameProcess()
    : DeviceAerosolProcess<MAMRenameProcess>(RenameProcess,
                                             "MAMRenameProcess")
    , is_cloudy{true} {}

void MAMRenameProcess::init_(const ModalAerosolConfig& config) {
  num_modes = config.num_modes();
  max_aer = num_modes;
  num_populations = config.num_aerosol_populations;

  // Reserve memory for private fields
  Kokkos::resize(dgnumlo, num_modes);
  Kokkos::resize(dgnumhi, num_modes);
  Kokkos::resize(dgnum, num_modes);
  Kokkos::resize(alnsg, num_modes);

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

  size_type num_pairs = 0;

  find_renaming_pairs_(config, dest_mode_of_mode_mapping, num_pairs,
                       size_factor, fmode_dist_tail_fac, volume2num_lo_relaxed,
                       volume2num_hi_relaxed, ln_diameter_tail_fac,
                       diameter_cutoff, ln_dia_cutoff, diameter_belowcutoff,
                       dryvol_smallest);
}

void MAMRenameProcess::find_renaming_pairs_(
    const ModalAerosolConfig& config,
    view_1d_int_type& dest_mode_of_mode_mapping, size_type& num_pairs,
    view_1d_pack_type& size_factor, view_1d_pack_type& fmode_dist_tail_fac,
    view_1d_pack_type& volume2num_lo_relaxed,
    view_1d_pack_type& volume2num_hi_relaxed,
    view_1d_pack_type& ln_diameter_tail_fac, view_1d_pack_type& diameter_cutoff,
    view_1d_pack_type& ln_dia_cutoff, view_1d_pack_type& diameter_belowcutoff,
    view_1d_pack_type& dryvol_smallest) const {
  const auto N = dest_mode_of_mode_mapping.extent(0);

  // number of pairs allowed to do inter-mode particle transfer
  // (e.g. if we have a pair "mode_1<-->mode_2", mode_1 and mode_2 can
  // participate in inter-mode aerosol particle transfer where like particles in
  // mode_1 can be transferred to mode_2 and vice-versa)
  Kokkos::parallel_reduce(
      N,
      KOKKOS_LAMBDA(const int i, decltype(num_pairs)& isum) {
        isum += static_cast<int>(dest_mode_of_mode_mapping(i) > 0);
      },
      num_pairs);

  if (num_pairs == 0) return;

  auto compute_size_factor = [](const auto& alnsg) {
    return pi_sixth * exp(4.5 * pow(alnsg, 2));
  };

  Kokkos::parallel_for(
      N, KOKKOS_LAMBDA(const int src_mode) {
        const auto& dest_mode_of_current_mode =
            dest_mode_of_mode_mapping(src_mode);

        const auto& alnsg_for_current_mode = alnsg[src_mode];
        const auto& alnsg_for_dest_mode = alnsg[dest_mode_of_current_mode];

        // Assign size factors for source and destination nodes.
        // Perform calculation only once, assignment twice.
        {
          const auto current_size_factor =
              compute_size_factor(alnsg_for_current_mode);
          size_factor[src_mode] = current_size_factor;
          size_factor[dest_mode_of_current_mode] = current_size_factor;
        }

        fmode_dist_tail_fac[src_mode] = sqrt_half / alnsg_for_current_mode;
        dryvol_smallest[src_mode] = smallest_dryvol_value;

        // Set relaxed limits of ratios for current source/dest mode pair
        {
          const Mode& mode = config.aerosol_modes[src_mode];
          volume2num_lo_relaxed[src_mode] =
              compute_relaxed_volume_to_num_ratio(mode, dgnumhi[src_mode]);
          volume2num_lo_relaxed[dest_mode_of_current_mode] =
              compute_relaxed_volume_to_num_ratio(
                  mode, dgnumhi[dest_mode_of_current_mode]);
        }

        // A factor for computing diameter at the tails of the distribution
        ln_diameter_tail_fac[src_mode] = 3.0 * pow(alnsg_for_current_mode, 2);

        // Cut-off (based on geometric mean) for making decision to do
        // inter-mode transfers
        //
        // TODO: use dummy values for _dgnum_, or assign to dgnum_low for the
        // moment. Have to figure out how to compute this. We will extract from
        // model at some point.
        {
          // Store in a temporary rather than access element of
          // _diameter_cutoff_ multiple times.
          auto diameter_cutoff_for_src_mode = [&]() -> haero::PackType {
            const auto sqrt_param_a =
                dgnum[src_mode] * exp(1.5 * pow(alnsg_for_current_mode, 2));

            const auto sqrt_param_b = dgnum[dest_mode_of_current_mode] *
                                      exp(1.5 * pow(alnsg_for_dest_mode, 2));

            return sqrt(sqrt_param_a * sqrt_param_b);
          }();

          diameter_cutoff[src_mode] = diameter_cutoff_for_src_mode;

          ln_dia_cutoff[src_mode] = log(diameter_cutoff_for_src_mode);

          diameter_belowcutoff[src_mode] = 0.99 * diameter_cutoff_for_src_mode;
        }
      });
}

}  // namespace haero
