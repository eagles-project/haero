#include <cstdio>

#include "catch2/catch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "haero/view_pack_helpers.hpp"
#include "kokkos/Kokkos_Core.hpp"

using namespace haero;

TEST_CASE("kokkos_teams", "simple_loop") {
  typedef ekat::ExeSpaceUtils<>::TeamPolicy::member_type TeamHandleType;
  const auto &teamPolicy = ekat::ExeSpaceUtils<>::get_default_team_policy(2, 0);
  Kokkos::parallel_for(
      teamPolicy, KOKKOS_LAMBDA(const TeamHandleType &team) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, 3u),
                             [&](const int &i) {});
      });
}

// Useful little function to allocate class on device. Probably move to header
// later.
#ifdef KOKKOS_ENABLE_CUDA
typedef Kokkos::CudaSpace MemSpace;
#else
typedef Kokkos::HostSpace MemSpace;
#endif
template <class T>
inline T *allocate_on_device() {
  const std::string debuggingName(typeid(T).name());
  T *t = static_cast<T *>(
      Kokkos::kokkos_malloc<MemSpace>(debuggingName + "_malloc", sizeof(T)));
  Kokkos::parallel_for(
      debuggingName + "_format", 1, KOKKOS_LAMBDA(const int) { new (t) T(); });
  return t;
}

template <class T>
inline T *allocate_on_device(const T &s) {
  const std::string debuggingName(typeid(T).name());
  T *t = static_cast<T *>(
      Kokkos::kokkos_malloc<MemSpace>(debuggingName + "_malloc", sizeof(T)));
  const T &S = s;  // Suck into lambda capture space.
  Kokkos::parallel_for(
      debuggingName + "_format", 1, KOKKOS_LAMBDA(const int) { new (t) T(S); });
  return t;
}

inline void free_on_device(void *t) { Kokkos::kokkos_free<MemSpace>(t); }

// Define some simple classes for the virutal function on device test.
class BaseClass {
 public:
  KOKKOS_INLINE_FUNCTION
  virtual ~BaseClass(){};
  KOKKOS_INLINE_FUNCTION
  virtual int virtual_function() const { return 0; }
};
class DerivedClass : public BaseClass {
 public:
  KOKKOS_INLINE_FUNCTION
  virtual ~DerivedClass(){};
  KOKKOS_INLINE_FUNCTION
  virtual int virtual_function() const { return 1; }
};

TEST_CASE("kokkos_reduce", "virtual_function") {
  BaseClass *base_class = allocate_on_device<DerivedClass>();

  int reduced = 0;
  Kokkos::Sum<int> reducer_int(reduced);
  Kokkos::parallel_reduce(
      "CallOnDevice", 3,
      KOKKOS_LAMBDA(const int, int &value) {
        value += base_class->virtual_function();
      },
      reducer_int);
  REQUIRE(reduced == 3);
  free_on_device(base_class);
}

// Define some simple classes for the virutal classes with Kokkos Views test.
class ViewsBaseClass {
 public:
  KOKKOS_INLINE_FUNCTION
  virtual ~ViewsBaseClass(){};
  KOKKOS_INLINE_FUNCTION
  virtual Real virtual_function(int) const { return 0; }
};
class ViewsDerivedClass : public ViewsBaseClass {
 public:
  ViewsDerivedClass(const std::vector<Real> &V)
      : device_view(vector_to_basic_1dview(V, "CopyVectorToDevice")) {}

  KOKKOS_INLINE_FUNCTION
  ViewsDerivedClass(const ViewsDerivedClass &T) : device_view(T.device_view) {}

  KOKKOS_INLINE_FUNCTION
  virtual ~ViewsDerivedClass(){};
  KOKKOS_INLINE_FUNCTION
  virtual Real virtual_function(int i) const { return device_view[i]; }
  view_1d_scalar_type device_view;
};

TEST_CASE("kokkos_parallel", "virtual_class_with_views") {
  static const int TEN = 10;
  std::vector<Real> V(TEN);
  for (int i = 0; i < TEN; ++i) V[i] = 2 * i + 1;
  ViewsDerivedClass derived_class(V);

  ViewsBaseClass *base_class = allocate_on_device(derived_class);

  Real reduced = 0;
  Kokkos::Sum<Real> reducer_real(reduced);
  Kokkos::parallel_reduce(
      "CallOnDevice", TEN,
      KOKKOS_LAMBDA(const int i, Real &value) {
        value += base_class->virtual_function(i);
      },
      reducer_real);
  REQUIRE(reduced == TEN * TEN);
  free_on_device(base_class);
}
