#include "haero/haero_config.hpp"
#include "haero/view_pack_helpers.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("view_pack_helpers", "") {

  const int mm = 12;
  const int nn = 9;
  std::vector<Real> plus_minus_one(nn);
  for (int j=0; j<nn; ++j) {
    plus_minus_one[j] = std::pow(-1,j);
  }

  std::vector<std::vector<Real>> vectors(mm, plus_minus_one);

  SECTION("basic_views") {

    /// Test vector -> view (1d)
    view_1d_scalar_type v1d = vector_to_basic_1dview(plus_minus_one, "plus_minus_one");
    auto host_v1d = Kokkos::create_mirror_view(v1d);
    Kokkos::deep_copy(host_v1d, v1d);
    int nerr = 0;
    for (int j=0; j<nn; ++j) {
      if (host_v1d(j) != plus_minus_one[j]) ++nerr;
    }
    REQUIRE (nerr == 0);

    /// Test view -> vector (1d)
    std::vector<Real> pm1 = view1d_to_vector(v1d);
    nerr = 0;
    for (int j=0; j<nn; ++j) {
      if (plus_minus_one[j] != pm1[j]) ++nerr;
    }
    REQUIRE (nerr == 0);

    /// Test 2d vectors -> view (2d)
    view_2d_scalar_type v2d = vector_to_basic_2dview(vectors, "vectors");
    auto host_v2d = Kokkos::create_mirror_view(v2d);
    Kokkos::deep_copy(host_v2d, v2d);
    nerr = 0;
    for (int i=0; i<mm; ++i) {
      for (int j=0; j<nn; ++j) {
        if (host_v2d(i,j) != vectors[i][j]) ++nerr;
      }
    }
    REQUIRE (nerr == 0);

  }
  SECTION("packed_views") {

    /// Test vector -> view (packed 1d)
    view_1d_pack_type v1d_pack = vector_to_packed_1dview(plus_minus_one, "plus_minus_one");
    auto host_v1d_pack = Kokkos::create_mirror_view(v1d_pack);
    Kokkos::deep_copy(host_v1d_pack, v1d_pack);
    int nerr = 0;
    for (int j=0; j<nn; ++j) {
      if (host_v1d_pack(pack_info::pack_idx(j))[pack_info::vec_idx(j)]
        != plus_minus_one[j] ) ++nerr;
    }
    REQUIRE (nerr == 0);

    /// Test view -> vector
    std::vector<Real> pm1pack = view1d_to_vector(v1d_pack);
    for (int j=0; j<nn; ++j) {
      if (pm1pack[j] != plus_minus_one[j]) ++nerr;
    }
    REQUIRE (nerr == 0);

    /// Test 2d vectors -> view
    view_2d_pack_type v2d_rowpack = vectors_to_row_packed_2dview(vectors, "vectors");
    auto host_v2d_rowpack = Kokkos::create_mirror_view(v2d_rowpack);
    Kokkos::deep_copy(host_v2d_rowpack, v2d_rowpack);
    for (int i=0; i<mm; ++i) {
      for (int j=0; j<nn; ++j) {
        if (host_v2d_rowpack(i,pack_info::pack_idx(j))[pack_info::vec_idx(j)]
          != vectors[i][j]) ++nerr;
      }
    }
    REQUIRE (nerr == 0);

    view_2d_pack_type v2d_colpack = vectors_to_col_packed_2dview(vectors, "vectors");
    auto host_v2d_colpack = Kokkos::create_mirror_view(v2d_colpack);
    Kokkos::deep_copy(host_v2d_colpack, v2d_colpack);
    for (int i=0; i<mm; ++i) {
      for (int j=0; j<nn; ++j) {
        if (host_v2d_colpack(pack_info::pack_idx(i), j)[pack_info::vec_idx(i)]
          != vectors[i][j]) ++nerr;
      }
    }
  }
}

