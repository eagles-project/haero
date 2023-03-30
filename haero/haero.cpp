// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#include <haero/haero.hpp>

namespace {

// memory allocation pool for standalone ColumnViews
class ColumnPool

ColumnPool *column_pool_ = nullptr;

void delete_column_pool() {
  delete column_pool_;
}

}

namespace haero {

ColumnView create_column_view(int num_levels) {
  if (!column_pool_) {
    column_pool_ = new ColumnPool();
    atexit(delete_column_pool);
  }
}

}

