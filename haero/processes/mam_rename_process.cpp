#include "haero/processes/mam_rename_process.hpp"

#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>

#include "haero/physical_constants.hpp"

namespace haero
{

namespace
{
  // Use the same types from the process
  template <typename T>
  using Container = MAMRenameProcess::Container<T>;
  using Integral = MAMRenameProcess::Integral;
  using Size = MAMRenameProcess::Size;

  // Use constants from haero::constants, calculate a few of our own.
  // TODO: Should _smallest_dryvol_value_ be in haero::constants as well?
  using haero::constants::pi_sixth;
  static constexpr Real frelax = 27.0;
  static constexpr Real sqrt_half = std::sqrt((Real)0.5);
  static constexpr Real smallest_dryvol_value = 1.0e-25;

  // Utility function to enumerate a container
  template <typename T>
  Container<std::pair<Size, T>> enumerate(Container<T> enumeratee)
  {
    using Pair = std::pair<Size, T>;
    Container<Pair> enumerated;
    for (int i = 0; i < enumeratee.size(); i++)
      enumerated.push_back(Pair(i, enumeratee[i]));
    return enumerated;
  }

  // This will be replaced by another method when the rest of the original
  // fortran has been ported. This is why _model_ is an unused parameter - it
  // will be used to calculate _dest_mode_of_mode_mapping_ when that calculation is
  // in place.
  //
  // \post _dest_mode_of_mode_mapping_ will have predetermined size of 4 due to
  // temporary workaround.
  void initialize_dest_mode_of_mode_mapping(Container<Integral> dest_mode_of_mode_mapping,
                                    const ModalAerosolConfig& model)
  {
    dest_mode_of_mode_mapping.assign({0, 1, 0, 0});
  }

  // Create a container and reserve memory for it.
  // \remark Why is this not in the stl yet?
  template <typename T = Real>
  Container<T> reserved_container(Size size)
  {
    Container<T> c;
    c.reserve(size);
    return c;
  }

  // TODO: What is a more descriptive name for _diameter_?
  // TODO: Where can I find a formula to reference for this calculation?
  static inline Real compute_relaxed_volume_to_num_ratio(
      const Real& alnsg_for_current_mode,
      const Real& diameter_for_current_mode)
  {
    static const auto unrelaxed_limit =
        1. / pi_sixth * std::pow(diameter_for_current_mode, 3)
        * std::exp(std::pow(4.5 * alnsg_for_current_mode, 2));

    return unrelaxed_limit / frelax;
  }

}  // namespace

MAMRenameProcess::MAMRenameProcess()
    : DeviceAerosolProcess<MAMRenameProcess>(RenameProcess,
                                             "MAMRenameProcess") {}


void MAMRenameProcess::init_(const ModalAerosolConfig& config)
{
  const auto& num_modes = config.num_modes();

  // Reserve memory for private fields
  dgnumlo.reserve(num_modes);
  dgnumhi.reserve(num_modes);
  dgnum.reserve(num_modes);
  alnsg.reserve(num_modes);
}

KOKKOS_FUNCTION
void MAMRenameProcess::run_(const ModalAerosolConfig& modal_aerosol_config,
                           Real t,
                           Real dt,
                           const Prognostics& prognostics,
                           const Atmosphere& atmosphere,
                           const Diagnostics& diagnostics,
                           Tendencies& tendencies) const
{
  const auto& num_modes = modal_aerosol_config.num_modes();

  // Create a vector with `num_modes` elements reserved
  auto reserved_num_modes_real_vector = [&] {
    return reserved_container<Real>(num_modes);
  };
  auto reserved_num_modes_int_vector = [&] {
    return reserved_container<Integral>(num_modes);
  };

  auto dest_mode_of_mode_mapping = reserved_num_modes_int_vector();

  initialize_dest_mode_of_mode_mapping(dest_mode_of_mode_mapping, modal_aerosol_config);

  auto size_factor = reserved_num_modes_real_vector();
  auto fmode_dist_tail_fac = reserved_num_modes_real_vector();
  auto volume2num_lo_relaxed = reserved_num_modes_real_vector();
  auto volume2num_hi_relaxed = reserved_num_modes_real_vector();
  auto ln_diameter_tail_fac = reserved_num_modes_real_vector();
  auto diameter_cutoff = reserved_num_modes_real_vector();
  auto diameter_belowcutoff = reserved_num_modes_real_vector();
  auto ln_dia_cutoff = reserved_num_modes_real_vector();
  auto dryvol_smallest = reserved_num_modes_real_vector();

  Size num_pairs = 0;

  find_renaming_pairs_(num_modes,
                       dest_mode_of_mode_mapping,
                       num_pairs,
                       size_factor,
                       fmode_dist_tail_fac,
                       volume2num_lo_relaxed,
                       volume2num_hi_relaxed,
                       ln_diameter_tail_fac,
                       diameter_cutoff,
                       ln_dia_cutoff,
                       diameter_belowcutoff,
                       dryvol_smallest);
}

void MAMRenameProcess::find_renaming_pairs_(
    const Size nmodes,
    const Container<Integral>& dest_mode_of_mode_mapping,
    Size& num_pairs,
    Container<Real>& size_factor,
    Container<Real>& fmode_dist_tail_fac,
    Container<Real>& volume2num_lo_relaxed,
    Container<Real>& volume2num_hi_relaxed,
    Container<Real>& ln_diameter_tail_fac,
    Container<Real>& diameter_cutoff,
    Container<Real>& ln_dia_cutoff,
    Container<Real>& diameter_belowcutoff,
    Container<Real>& dryvol_smallest) const
{

  // number of pairs allowed to do inter-mode particle transfer
  // (e.g. if we have a pair "mode_1<-->mode_2", mode_1 and mode_2 can
  // participate in inter-mode aerosol particle transfer where like particles in
  // mode_1 can be transferred to mode_2 and vice-versa)
  //
  // Let us assume there are none to start with.
  num_pairs = 0;

  // Filter out values <= 0, as they can't have a pair
  Container<Integral> filtered_mode_mappings;
  std::copy_if(dest_mode_of_mode_mapping.begin(),
               dest_mode_of_mode_mapping.end(),
               std::back_inserter(filtered_mode_mappings),
               std::bind(std::greater<>(), std::placeholders::_1, 0));

  auto compute_size_factor = [](const auto& alnsg) {
    return pi_sixth * std::exp(4.5 * std::pow(alnsg, 2));
  };

  // If we ever bump to c++17, this should use tuple decomposition
  for (const auto& tup : enumerate(filtered_mode_mappings)) {
    const auto src_mode = std::get<0>(tup);
    const auto dest_mode_of_current_mode = std::get<1>(tup);

    const auto& alnsg_for_current_mode = alnsg[src_mode];
    const auto& alnsg_for_dest_mode = alnsg[dest_mode_of_current_mode];

    // Each value in the filtered container of destination modes indicates a
    // pair has been found
    num_pairs++;

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
      volume2num_lo_relaxed[src_mode] =
          compute_relaxed_volume_to_num_ratio(alnsg_for_current_mode,
                                              dgnumhi[src_mode]);
      volume2num_lo_relaxed[dest_mode_of_current_mode] =
          compute_relaxed_volume_to_num_ratio(alnsg_for_current_mode,
                                              dgnumhi[dest_mode_of_current_mode]);
    }

    // A factor for computing diameter at the tails of the distribution
    ln_diameter_tail_fac[src_mode] = 3.0 * std::pow(alnsg_for_current_mode, 2);

    // Cut-off (based on geometric mean) for making decision to do inter-mode
    // transfers
    //
    // TODO: use dummy values for _dgnum_, or assign to dgnum_low for the
    // moment. Have to figure out how to compute this. We will extract from
    // model at some point.
    {
      // Store in a temporary rather than access element of _diameter_cutoff_
      // multiple times.
      auto diameter_cutoff_for_src_mode = [&] () -> Real {
        const auto sqrt_param_a =
            dgnum[src_mode]
            * std::exp(1.5 * std::pow(alnsg_for_current_mode, 2));

        const auto sqrt_param_b =
            dgnum[dest_mode_of_current_mode]
            * std::exp(1.5 * std::pow(alnsg_for_dest_mode, 2));

        return std::sqrt(sqrt_param_a * sqrt_param_b);
      }();

      diameter_cutoff[src_mode] = diameter_cutoff_for_src_mode;

      ln_dia_cutoff[src_mode] = std::log(diameter_cutoff_for_src_mode);

      diameter_belowcutoff[src_mode] = 0.99 * diameter_cutoff_for_src_mode;
    }
  }
}

}  // namespace haero
