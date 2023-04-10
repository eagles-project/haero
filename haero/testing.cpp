// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#include "testing.hpp"

namespace haero {

namespace {

using Real = haero::Real;

// A simple memory allocation pool for standalone ColumnViews to be used in
// (e.g.) unit tests. A ColumnPool manages a number of ColumnViews with a fixed
// number of vertical levels.
class ColumnPool {
  size_t num_levels_; // number of vertical levels per column (fixed)
  size_t num_cols_;   // number of allocated columns
  std::vector<size_t> col_offsets_; // offsets of columns in memory_
  std::vector<int> col_in_use_;     // an array indicating whether a column is
                                    // being used

  Real *memory_; // memory pool itself (allocated on device)
  size_t size_;  // size of memory pool in bytes
public:
  // constructs a column pool with the given initial number of columns, each
  // with the given number of vertical levels.`
  ColumnPool(size_t num_vertical_levels, size_t initial_num_columns = 64)
      : num_levels_(num_vertical_levels), num_cols_(initial_num_columns),
        col_in_use_(initial_num_columns, 0),
        memory_(reinterpret_cast<Real *>(Kokkos::kokkos_malloc(
            "Column pool",
            sizeof(Real) * num_vertical_levels * initial_num_columns))),
        size_(sizeof(Real) * num_vertical_levels * initial_num_columns) {}

  // no copying of the pool

  // destructor
  ~ColumnPool() { Kokkos::kokkos_free(memory_); }

  // returns a "fresh" (unused) ColumnView from the ColumnPool, marking it as
  // used, and allocating additional memory if needed)
  ColumnView column_view() {
    // find the first unused column
    size_t i;
    for (i = 0; i < num_cols_; ++i) {
      if (!col_in_use_[i])
        break;
    }
    if (i == num_cols_) { // all columns in the pool are in use!
      // double the number of allocated columns in the pool
      num_cols_ *= 2;
      col_in_use_.resize(num_cols_, 0);
      memory_ = reinterpret_cast<Real *>(Kokkos::kokkos_realloc(
          memory_, sizeof(Real) * num_levels_ * num_cols_));
    }

    col_in_use_[i] = 1;
    size_t offset = sizeof(Real) * i * num_levels_;
    return ColumnView(&memory_[offset], num_levels_);
  }
};

// column pools, organized by column resolution
std::map<size_t, std::unique_ptr<ColumnPool>> pools_{};

} // namespace

namespace testing {

Atmosphere create_atmosphere(int num_levels, Real pblh) {
  Atmosphere atm(num_levels, pblh);
  atm.temperature = create_column_view(num_levels);
  atm.pressure = create_column_view(num_levels);
  atm.vapor_mixing_ratio = create_column_view(num_levels);
  atm.height = create_column_view(num_levels);
  atm.hydrostatic_dp = create_column_view(num_levels);
  atm.cloud_fraction = create_column_view(num_levels);
  atm.updraft_vel_ice_nucleation = create_column_view(num_levels);
  return atm;
}

ColumnView create_column_view(int num_levels) {
  // find a column pool for the given number of vertical levels
  auto iter = pools_.find(num_levels);
  if (iter == pools_.end()) {
    auto result = pools_.emplace(
        num_levels, std::unique_ptr<ColumnPool>(new ColumnPool(num_levels)));
    iter = result.first;
  }
  return iter->second->column_view();
}

} // namespace testing

} // namespace haero
