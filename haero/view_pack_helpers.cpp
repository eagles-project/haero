#include "view_pack_helpers.hpp"

namespace haero {

view_1d_scalar_type vector_to_basic_1dview(const std::vector<Real>& vector,
                                           const std::string& view_name) {
  view_1d_scalar_type result(view_name, vector.size());
  auto hm = Kokkos::create_mirror_view(result);
  for (int i = 0; i < vector.size(); ++i) {
    hm(i) = vector[i];
  }
  Kokkos::deep_copy(result, hm);
  return result;
}

view_1d_int_type vector_to_basic_1dview(const std::vector<int>& vector,
                                        const std::string& view_name) {
  view_1d_int_type result(view_name, vector.size());
  auto hm = Kokkos::create_mirror_view(result);
  for (int i = 0; i < vector.size(); ++i) {
    hm(i) = vector[i];
  }
  Kokkos::deep_copy(result, hm);
  return result;
}

view_1d_pack_type vector_to_packed_1dview(const std::vector<Real>& vector,
                                          const std::string& view_name) {
  const int nn = vector.size();
  view_1d_pack_type result(view_name, pack_info::num_packs(nn));
  auto hm = Kokkos::create_mirror_view(result);
  for (int i = 0; i < nn; ++i) {
    hm(pack_info::pack_idx(i))[pack_info::vec_idx(i)] = vector[i];
  }
  Kokkos::deep_copy(result, hm);
  return result;
}

view_2d_scalar_type vector_to_basic_2dview(
    const std::vector<std::vector<Real>>& vectors,
    const std::string& view_name) {
  const int mm = vectors.size();
  const int nn = vectors[0].size();
  view_2d_scalar_type result(view_name, mm, nn);
  auto hm = Kokkos::create_mirror_view(result);
  for (int i = 0; i < mm; ++i) {
    for (int j = 0; j < nn; ++j) {
      hm(i, j) = vectors[i][j];
    }
  }
  Kokkos::deep_copy(result, hm);
  return result;
}

view_2d_pack_type vectors_to_row_packed_2dview(
    const std::vector<std::vector<Real>>& vectors,
    const std::string& view_name) {
  const int mm = vectors.size();
  const int nn = vectors[0].size();
  view_2d_pack_type result(view_name, mm, pack_info::num_packs(nn));
  auto hm = Kokkos::create_mirror_view(result);
  for (int i = 0; i < mm; ++i) {
    for (int j = 0; j < nn; ++j) {
      hm(i, pack_info::pack_idx(j))[pack_info::vec_idx(j)] = vectors[i][j];
    }
  }
  Kokkos::deep_copy(result, hm);
  return result;
}

}  // namespace haero
