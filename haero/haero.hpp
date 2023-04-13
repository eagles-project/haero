// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HAERO_HPP
#define HAERO_HPP

#include <haero/haero_config.hpp>

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

/// Device and host types
using DeviceType = ekat::KokkosTypes<ekat::DefaultDevice>;
using HostType = ekat::KokkosTypes<ekat::HostDevice>;

/// Execution space (based on DeviceType)
using ExecutionSpace = typename DeviceType::ExeSpace;

/// Kokkos team policy used to perform parallel dispatches of aerosol processes
/// over several columns.
using ThreadTeamPolicy = typename DeviceType::TeamPolicy;

/// Kokkos team type for the Kokkos::TeamThreadRange parallel dispatch policy.
/// This type is often called a TeamMember in other Kokkos codes, but it's
/// actually not a member of a team--it's a member of a league: in other words,
/// a team!
using ThreadTeam = typename DeviceType::MemberType;

/// A TracersView іs a rank 3 Kokkos View with the following indices:
/// 1. A tracer index identifying an advected quantity
/// 2. A column index identifying a unique atmospheric column
/// 3. A level index identifying a unique vertical level in the column
using TracersView = typename DeviceType::view_3d<Real>;

/// A DiagnosticsView іs a rank 3 Kokkos View with the following indices:
/// 1. A diagnostic index identifying a diagnostic quantity
/// 2. A column index identifying a unique atmospheric column
/// 3. A level index identifying a unique vertical level in the column
using DiagnosticsView = typename DeviceType::view_3d<Real>;

/// A ColumnView is a rank-1 Kokkos View whose single index identifies a
/// unique vertical level. ColumnViews are unmanaged, meaning that they don't
/// own their storage. This allows a host model to manage column data for
/// any Haero aerosol packages.
using ColumnView = ekat::Unmanaged<typename DeviceType::view_1d<Real>>;

// Returns haero's version string.
const char *version();

// Returns haero's git revision hash, or "unknown" if not found.
const char *revision();

// Returns true iff this build has changes that weren't committed to the git
// repo.
bool has_uncommitted_changes();

} // end namespace haero

#endif
