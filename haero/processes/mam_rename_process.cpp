#include "haero/processes/mam_rename_process.hpp"

#include <cassert>
#include <iostream>

namespace haero
{

namespace
{
  using Container = std::vector<Real>;

  // This will be replaced by another method when the rest of the original
  // fortran has been ported.
  void initialize_dest_mode_of_mode(Container dest_mode_of_mode,
                                    const ModalAerosolConfig& model)
  {
    assert(dest_mode_of_mode.size() == 4
           && "This test is only configured to use num_modes==4!");
    dest_mode_of_mode.assign({0, 1, 0, 0});
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

  // Initialize vectors used in these methods to their maximum size.
  // Why is this not in the stl yet?
  auto reserved_vector = [](std::size_t sz) {
    Container v;
    v.reserve(sz);
    return v;
  };

  // Create a vector with `num_modes` elements reserved
  auto reserved_num_modes_vector = [&] { return reserved_vector(num_modes); };

  auto dest_mode_of_mode = reserved_num_modes_vector();

  initialize_dest_mode_of_mode(dest_mode_of_mode, modal_aerosol_config);

  auto size_factor = reserved_num_modes_vector();
  auto fmode_dist_tail_fac = reserved_num_modes_vector();
  auto v2n_lo_rlx = reserved_num_modes_vector();
  auto v2n_hi_rlx = reserved_num_modes_vector();
  auto ln_diameter_tail_fac = reserved_num_modes_vector();
  auto diameter_cutoff = reserved_num_modes_vector();
  auto diameter_belowcutoff = reserved_num_modes_vector();
  auto ln_dia_cutoff = reserved_num_modes_vector();
  auto dryvol_smallest = reserved_num_modes_vector();

  std::size_t num_pairs = 0;

  find_renaming_pairs_(num_modes,
                       dest_mode_of_mode,
                       num_pairs,
                       size_factor,
                       fmode_dist_tail_fac,
                       v2n_lo_rlx,
                       v2n_hi_rlx,
                       ln_diameter_tail_fac,
                       diameter_cutoff,
                       ln_dia_cutoff,
                       diameter_belowcutoff,
                       dryvol_smallest);
}

void MAMRenameProcess::find_renaming_pairs_(
    const std::size_t nmodes,
    const Container& dest_mode_of_mode,
    const std::size_t num_pairs,
    const Container& sz_factor,
    const Container& fmode_dist_tail_fac,
    const Container& v2n_lo_rlx,
    const Container& v2n_hi_rlx,
    const Container& ln_diameter_tail_fac,
    const Container& diameter_cutoff,
    const Container& ln_dia_cutoff,
    const Container& diameter_belowcutoff,
    const Container& dryvol_smallest) const
{
  std::cout << "Running find_renaming_pairs_\n";
}

}  // namespace haero
