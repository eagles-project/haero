#ifndef HAERO_HPP
#define HAERO_HPP

#include "haero/haero_config.hpp"

// Cuda "C++" can't handle lambdas consisting of private/protected methods
// This seems awful, but other solutions seem to involve dancing carefully
// around various areas in code, or sacrificing encapsulation in CPU builds.
#ifdef __CUDACC__
#define protected public
#define private public
#endif

namespace haero {

/// MemorySpace refers to the memory space on the device.
#ifdef KOKKOS_ENABLE_CUDA
typedef Kokkos::CudaSpace MemorySpace;
#else
typedef Kokkos::HostSpace MemorySpace;
#endif

/// Helpers for pack/array indexing
using PackInfo = ekat::PackInfo<HAERO_PACK_SIZE>;

/// Pack data structures for SIMD
using PackType = ekat::Pack<Real, HAERO_PACK_SIZE>;
using IntPackType = ekat::Pack<int, HAERO_PACK_SIZE>;
using MaskType = ekat::Mask<HAERO_PACK_SIZE>;

/// Device and host types
using DeviceType = ekat::KokkosTypes<ekat::DefaultDevice>;
using HostType   = ekat::KokkosTypes<ekat::HostDevice>;

/// Execution space (based on DeviceType)
using ExecutionSpace = typename DeviceType::ExeSpace;

/// Kokkos team policy used to perform parallel dispatches of aerosol processes
/// over several columns.
using TeamPolicy = typename DeviceType::TeamPolicy;

/// Kokkos team type for the Kokkos::TeamThreadRange parallel dispatch policy.
/// This type is often called a TeamMember in other Kokkos codes, but it's
/// actually not a member of a team--it's a member of a league: in other words,
/// a team!
using TeamType = typename DeviceType::MemberType;

/// A TracersView іs a rank 3 Kokkos View with the following indices:
/// 1. A tracer index identifying an advected quantity
/// 2. A column index identifying a unique atmospheric column
/// 3. A level index identifying a unique pack containing a set of data
///    associated with adjacent vertical levels.
using TracersView = typename DeviceType::view_3d<PackType>;

/// A DiagnosticsView іs a rank 3 Kokkos View with the following indices:
/// 1. A diagnostic index identifying a diagnostic quantity
/// 2. A column index identifying a unique atmospheric column
/// 3. A level index identifying a unique pack containing a set of data
///    associated with adjacent vertical levels.
using DiagnosticsView = typename DeviceType::view_3d<PackType>;

/// A ColumnView is a rank-1 Kokkos View whose single index identifies a
/// unique pack containing a set of data associated with adjacent vertical
/// levels.
using ColumnView = typename DeviceType::view_1d<PackType>;

// Returns haero's version string.
const char* version();

// Returns haero's git revision hash, or "unknown" if not found.
const char* revision();

// Returns true iff this build has changes that weren't committed to the git
// repo.
bool has_uncommitted_changes();

// Below are min/max functions for Haero's Real type and for mixed Real/Pack
// arguments.

inline Real min(Real a, Real b) {
  return (a < b) ? a : b;
}
inline PackType min(Real a, const PackType& b) {
  PackType result = b;
  result.set(a < b, a);
  return result;
}
inline PackType min(const PackType& a, Real b) {
  PackType result = a;
  result.set(b < a, b);
  return result;
}

inline Real max(Real a, Real b) {
  return (a > b) ? a : b;
}
inline PackType max(Real a, const PackType& b) {
  PackType result = b;
  result.set(a > b, a);
  return result;
}
inline PackType max(const PackType& a, Real b) {
  PackType result = a;
  result.set(b > a, b);
  return result;
}


}  // end namespace haero

#endif
