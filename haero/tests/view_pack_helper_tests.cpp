#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "catch2/catch.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "haero/haero.hpp"
#include "haero/view_pack_helpers.hpp"

using namespace haero;

TEST_CASE("view_pack_helpers", "") {
  const int mm = 12;
  const int nn = 9;
  std::vector<Real> plus_minus_one(nn);
  for (int j = 0; j < nn; ++j) {
    plus_minus_one[j] = std::pow(-1, j);
  }

  std::vector<std::vector<Real>> vectors(mm, plus_minus_one);

  SECTION("basic_views") {
    /// Test vector -> view (1d)
    RealView1D v1d =
        vector_to_basic_1dview(plus_minus_one, "plus_minus_one");
    auto host_v1d = Kokkos::create_mirror_view(v1d);
    Kokkos::deep_copy(host_v1d, v1d);
    int nerr = 0;
    for (int j = 0; j < nn; ++j) {
      if (host_v1d(j) != plus_minus_one[j]) ++nerr;
    }
    REQUIRE(nerr == 0);

    /// Test view -> vector (1d)
    std::vector<Real> pm1 = view1d_to_vector(v1d);
    nerr = 0;
    for (int j = 0; j < nn; ++j) {
      if (plus_minus_one[j] != pm1[j]) ++nerr;
    }
    REQUIRE(nerr == 0);

    /// Test 2d vectors -> view (2d)
    RealView2D v2d = vector_to_basic_2dview(vectors, "vectors");
    auto host_v2d = Kokkos::create_mirror_view(v2d);
    Kokkos::deep_copy(host_v2d, v2d);
    nerr = 0;
    for (int i = 0; i < mm; ++i) {
      for (int j = 0; j < nn; ++j) {
        if (host_v2d(i, j) != vectors[i][j]) ++nerr;
      }
    }
    REQUIRE(nerr == 0);
  }
  SECTION("packed_views") {
    /// Test vector -> view (packed 1d)
    PackRealView1D v1d_pack =
        vector_to_packed_1dview(plus_minus_one, "plus_minus_one");
    auto host_v1d_pack = Kokkos::create_mirror_view(v1d_pack);
    Kokkos::deep_copy(host_v1d_pack, v1d_pack);
    int nerr = 0;
    for (int j = 0; j < nn; ++j) {
      if (host_v1d_pack(PackInfo::pack_idx(j))[PackInfo::vec_idx(j)] !=
          plus_minus_one[j])
        ++nerr;
    }
    REQUIRE(nerr == 0);

    /// Test view -> vector
    std::vector<Real> pm1pack = view1d_to_vector(v1d_pack, nn);
    for (int j = 0; j < nn; ++j) {
      if (pm1pack[j] != plus_minus_one[j]) ++nerr;
    }
    if (nerr > 0) {
      std::cout << "pm1pack [size " << pm1pack.size() << "] = (";
      for (int j = 0; j < pm1pack.size(); ++j) {
        std::cout << pm1pack[j] << (j < pm1pack.size() - 1 ? " " : ")\n");
      }
    }
    REQUIRE(pm1pack.size() == nn);
    REQUIRE(nerr == 0);

    /// Test 2d vectors -> view
    PackRealView2D v2d_rowpack =
        vectors_to_row_packed_2dview(vectors, "vectors");
    auto host_v2d_rowpack = Kokkos::create_mirror_view(v2d_rowpack);
    Kokkos::deep_copy(host_v2d_rowpack, v2d_rowpack);
    for (int i = 0; i < mm; ++i) {
      for (int j = 0; j < nn; ++j) {
        if (host_v2d_rowpack(i, PackInfo::pack_idx(j))[PackInfo::vec_idx(j)] !=
            vectors[i][j])
          ++nerr;
      }
    }
    REQUIRE(nerr == 0);

    REQUIRE_FALSE(std::is_arithmetic<PackType>::value);
    REQUIRE(
        std::is_arithmetic<ekat::ScalarTraits<PackType>::scalar_type>::value);
    REQUIRE(std::is_arithmetic<ekat::ScalarTraits<Real>::scalar_type>::value);
  }
}

template <int PackSize>
struct PackViewTest {
  typedef ekat::Pack<Real, PackSize> test_pack_type;
  typedef ekat::PackInfo<PackSize> test_pack_info;

  Kokkos::View<test_pack_type*> level_view;
  Kokkos::View<test_pack_type*> intfc_view;
  int nlev;

  PackViewTest(const int nl)
      : level_view("level_view", test_pack_info::num_packs(nl)),
        intfc_view("intfc_view", test_pack_info::num_packs(nl + 1)),
        nlev(nl) {}

  std::string info_string() const;
};

TEST_CASE("scalarized views", "") {
  const int nlev = 72;

  SECTION("pack size = 1") {
    static constexpr int tc_pack_size = 1;
    PackViewTest<tc_pack_size> testcase(nlev);
    std::cout << testcase.info_string();

    REQUIRE(ekat::scalarize(testcase.level_view).extent(0) ==
            tc_pack_size *
                PackViewTest<tc_pack_size>::test_pack_info::num_packs(nlev));
    REQUIRE(
        ekat::scalarize(testcase.intfc_view).extent(0) ==
        tc_pack_size *
            PackViewTest<tc_pack_size>::test_pack_info::num_packs(nlev + 1));
  }
  SECTION("pack size = 2") {
    static constexpr int tc_pack_size = 2;

    PackViewTest<tc_pack_size> testcase(nlev);
    std::cout << testcase.info_string();

    REQUIRE(ekat::scalarize(testcase.level_view).extent(0) ==
            tc_pack_size *
                PackViewTest<tc_pack_size>::test_pack_info::num_packs(nlev));
    REQUIRE(
        ekat::scalarize(testcase.intfc_view).extent(0) ==
        tc_pack_size *
            PackViewTest<tc_pack_size>::test_pack_info::num_packs(nlev + 1));
  }
  SECTION("pack size = 4") {
    static constexpr int tc_pack_size = 4;

    PackViewTest<tc_pack_size> testcase(nlev);
    std::cout << testcase.info_string();

    REQUIRE(ekat::scalarize(testcase.level_view).extent(0) ==
            tc_pack_size *
                PackViewTest<tc_pack_size>::test_pack_info::num_packs(nlev));
    REQUIRE(
        ekat::scalarize(testcase.intfc_view).extent(0) ==
        tc_pack_size *
            PackViewTest<tc_pack_size>::test_pack_info::num_packs(nlev + 1));
  }
}

TEST_CASE("Haero view type basics", "") {
  SECTION("ColumnView") {
    using std::isnan;

    ColumnView empty_view;
    std::cout << "An empty view (default-constructed) has extent "
              << empty_view.extent(0) << "\n";
    REQUIRE(empty_view.extent(0) == 0);

    ColumnView view_18("view_18", PackInfo::num_packs(18));

    auto s_view_18 = ekat::scalarize(view_18);
    auto h_view_18 = Kokkos::create_mirror_view(s_view_18);
    Kokkos::deep_copy(h_view_18, s_view_18);
    for (int i = 0; i < s_view_18.extent(0); ++i) {
      REQUIRE(isnan(h_view_18(i)));
    }

    zero_init(view_18, 18);
    Kokkos::deep_copy(h_view_18, s_view_18);
    for (int i = 0; i < 18; ++i) {
      REQUIRE(h_view_18(i) == 0);
    }
    for (int i=18; i<h_view_18.extent(0); ++i) {
      REQUIRE( isnan(h_view_18(i)) );
    }

  }
}

template <int PackSize>
std::string PackViewTest<PackSize>::info_string() const {
  std::ostringstream ss;
  ss << "PackViewTest info:\n";
  ss << "\t"
     << "nlev = " << nlev << "\n";
  ss << "\t"
     << "pack size = " << PackSize << "\n";
  ss << "\t"
     << "num_packs (levels) = " << test_pack_info::num_packs(nlev) << "\n";
  ss << "\t"
     << "level_view.extent(0) = " << level_view.extent(0) << "\n";
  ss << "\t"
     << "ekat::scalarize(level_view).extent(0) = "
     << ekat::scalarize(level_view).extent(0) << "\n";
  ss << "\t"
     << "num_packs (intfcs) = " << test_pack_info::num_packs(nlev + 1) << "\n";
  ss << "\t"
     << "intfc_view.extent(0) = " << intfc_view.extent(0) << "\n";
  ss << "\t"
     << "ekat::scalarize(intfc_view).extent(0) = "
     << ekat::scalarize(intfc_view).extent(0) << "\n";
  return ss.str();
}
