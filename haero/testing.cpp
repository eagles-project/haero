// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#include "testing.hpp"

#include <ekat/ekat_session.hpp>

namespace haero {

namespace {

using Real = haero::Real;

// A simple memory allocation pool for standalone ColumnViews to be used in
// (e.g.) unit tests. A ColumnPool manages a number of ColumnViews with a fixed
// number of vertical levels.
class ColumnPool {
  size_t num_levels_;          // number of vertical levels per column (fixed)
  size_t num_cols_;            // number of allocated columns
  std::vector<int> col_used_;  // columns that are being used already
  std::vector<Real *> memory_; // per-column memory itself (allocated on device)
public:
  // constructs a column pool with the given initial number of columns, each
  // with the given number of vertical levels.`
  ColumnPool(size_t num_vertical_levels, size_t initial_num_columns = 64)
      : num_levels_(num_vertical_levels), num_cols_(initial_num_columns),
        col_used_(initial_num_columns, 0),
        memory_(initial_num_columns, nullptr) {
    for (size_t i = 0; i < num_cols_; ++i) {
      memory_[i] = reinterpret_cast<Real *>(Kokkos::kokkos_malloc(
          "Column pool", sizeof(Real) * num_vertical_levels));
    }
  }

  // no copying of the pool

  // destructor
  ~ColumnPool() {
    for (size_t i = 0; i < num_cols_; ++i) {
      Kokkos::kokkos_free(memory_[i]);
    }
  }

  // returns a "fresh" (unused) ColumnView from the ColumnPool, marking it as
  // used, and allocating additional memory if needed)
  ColumnView column_view() {
    // find the first unused column
    size_t i;
    for (i = 0; i < num_cols_; ++i) {
      if (!col_used_[i])
        break;
    }
    if (i == num_cols_) { // all columns in the pool are in use!
      // double the number of allocated columns in the pool
      size_t new_num_cols = 2 * num_cols_;
      col_used_.resize(new_num_cols, 0);
      memory_.resize(new_num_cols, nullptr);
      for (size_t i = num_cols_; i < new_num_cols; ++i) {
        memory_[i] = reinterpret_cast<Real *>(
            Kokkos::kokkos_malloc(sizeof(Real) * num_levels_));
      }
      num_cols_ = new_num_cols;
    }

    col_used_[i] = 1;
    return ColumnView(memory_[i], num_levels_);
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

void finalize() { pools_.clear(); }

} // namespace testing

} // namespace haero

//------------------------------------------------------------------------
// EKAT test session initialization and finalization overrides
//------------------------------------------------------------------------
// The following functions are called at the beginning and the end of an
// EKAT test session. When calling EkatCreateUnitTest, you must specify the
// option EXCLUDE_TEST_SESSION which prevents EKAT from using its own
// default implementations.
//------------------------------------------------------------------------

// This implementation of ekat_initialize_test_session is identical to the
// default provided by EKAT.
void ekat_initialize_test_session(int argc, char **argv,
                                  const bool print_config) {
  ekat::initialize_ekat_session(argc, argv, print_config);
}

// This implementation of ekat_finalize_test_session calls
// haero::testing::finalize() to deallocate all ColumnView pools.
void ekat_finalize_test_session() {
  haero::testing::finalize();
  ekat::finalize_ekat_session();
}
