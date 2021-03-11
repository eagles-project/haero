#ifndef HAERO_HPP
#define HAERO_HPP

#include "haero/haero_config.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"

namespace haero {

/// This is the device on which the data is stored.
using DeviceType = ekat::KokkosTypes<ekat::DefaultDevice>;
using HostType   = ekat::KokkosTypes<ekat::HostDevice>;

/// Pack data structure for SIMD.
using PackType = ekat::Pack<Real, HAERO_PACK_SIZE>;

/// Helpers for pack/array indexing
using PackInfo = ekat::PackInfo<HAERO_PACK_SIZE>;

/// Views of this type store packs in a vertical column, with
/// ceil(num_levels/HAERO_PACK_SIZE) packs spanning the column's vertical
/// levels. If k identifies a vertical level, the corresponding data is
/// accessible at `view(pack_info::pack_idx(k))[pack_info::vec_idx(k)]`.
/// See Chapter 3 of the design document for an explanation of packs.
using ColumnView = DeviceType::view_1d<PackType>;

/// Views of this type store packs defined for each species population (for
/// aerosols or gases), with ceil(num_levels/HAERO_PACK_SIZE) packs spanning
/// the column's vertical levels. If p identifies a species population (aerosol
/// or gas), and k identifies a vertical level,
/// the corresponding data is accessible at
/// `view(p, pack_info::pack_idx(k))[pack_info::vec_idx(k)]`.
/// See Chapter 3 of the design document for an explanation of packs.
using SpeciesColumnView = DeviceType::view_2d<PackType>;

/// Views of this type store packs defined for each mode, with
/// ceil(num_levels/HAERO_PACK_SIZE) packs spanning the column's vertical
/// levels. If m identifies a mode and k identifies a vertical level,
/// the corresponding data is accessible at
/// `view(m, pack_info::pack_idx(k))[pack_info::vec_idx(k)]`.
/// See Chapter 3 of the design document for an explanation of packs.
using ModalColumnView = DeviceType::view_2d<PackType>;

// Returns haero's version string.
const char* version();

// Returns haero's git revision hash, or "unknown" if not found.
const char* revision();

// Returns true iff this build has changes that weren't committed to the git
// repo.
bool has_uncommitted_changes();

} // end namespace haero

#endif
