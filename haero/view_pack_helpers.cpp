#include "view_pack_helpers.hpp"

namespace haero {

RealView1D vector_to_basic_1dview(const std::vector<Real>& vector,
                                           const std::string& view_name) {
  RealView1D result(view_name, vector.size());
  auto hm = Kokkos::create_mirror_view(result);
  for (int i = 0; i < vector.size(); ++i) {
    hm(i) = vector[i];
  }
  Kokkos::deep_copy(result, hm);
  return result;
}

IntView1D vector_to_basic_1dview(const std::vector<int>& vector,
                                        const std::string& view_name) {
  IntView1D result(view_name, vector.size());
  auto hm = Kokkos::create_mirror_view(result);
  for (int i = 0; i < vector.size(); ++i) {
    hm(i) = vector[i];
  }
  Kokkos::deep_copy(result, hm);
  return result;
}

PackRealView1D vector_to_packed_1dview(const std::vector<Real>& vector,
                                          const std::string& view_name) {
  const int nn = vector.size();
  PackRealView1D result(view_name, PackInfo::num_packs(nn));
  auto hm = Kokkos::create_mirror_view(result);
  for (int i = 0; i < nn; ++i) {
    hm(PackInfo::pack_idx(i))[PackInfo::vec_idx(i)] = vector[i];
  }
  Kokkos::deep_copy(result, hm);
  return result;
}

RealView2D vector_to_basic_2dview(
    const std::vector<std::vector<Real>>& vectors,
    const std::string& view_name) {
  const int mm = vectors.size();
  const int nn = vectors[0].size();
  RealView2D result(view_name, mm, nn);
  auto hm = Kokkos::create_mirror_view(result);
  for (int i = 0; i < mm; ++i) {
    for (int j = 0; j < nn; ++j) {
      hm(i, j) = vectors[i][j];
    }
  }
  Kokkos::deep_copy(result, hm);
  return result;
}

PackRealView2D vectors_to_row_packed_2dview(
    const std::vector<std::vector<Real>>& vectors,
    const std::string& view_name) {
  const int mm = vectors.size();
  const int nn = vectors[0].size();
  PackRealView2D result(view_name, mm, PackInfo::num_packs(nn));
  auto hm = Kokkos::create_mirror_view(result);
  for (int i = 0; i < mm; ++i) {
    for (int j = 0; j < nn; ++j) {
      hm(i, PackInfo::pack_idx(j))[PackInfo::vec_idx(j)] = vectors[i][j];
    }
  }
  Kokkos::deep_copy(result, hm);
  return result;
}

}  // namespace haero
