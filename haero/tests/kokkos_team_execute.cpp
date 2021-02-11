#include "catch2/catch.hpp"
#include <Kokkos_Core.hpp>
#include <cstdio>


TEST_CASE("kokkos_teams", "simple_loop") {
typedef Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type TeamHandleType;
const auto& teamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>(9u, Kokkos::AUTO);
Kokkos::parallel_for(teamPolicy,
                    KOKKOS_LAMBDA(const TeamHandleType& team)
                    {
                      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, 3u), [&] (const int& i)
                      {
                      });
                    });
}





// Define some simple classes for the virutal function on device test.
class BaseClass {
public:
  KOKKOS_INLINE_FUNCTION
  virtual ~BaseClass() {};
  KOKKOS_INLINE_FUNCTION
  virtual int virtual_function() const {
    return 0;
  }
};
class DerivedClass : public BaseClass  {
public :
  KOKKOS_INLINE_FUNCTION
  virtual ~DerivedClass() {};
  KOKKOS_INLINE_FUNCTION
  virtual int virtual_function() const {
    return 1;
  }
};

// Useful little function to allocate class on device. Probably move to header later.
#ifdef KOKKOS_ENABLE_CUDA
  typedef Kokkos::CudaSpace MemSpace;
#else
  typedef Kokkos::HostSpace MemSpace;
#endif
template <class T>
inline T* allocate_on_device() {
  const std::string debuggingName(typeid(T).name());
  T* t = static_cast<T*>(Kokkos::kokkos_malloc<MemSpace>(debuggingName+"_malloc",sizeof(T)));
  Kokkos::parallel_for(debuggingName+"_format", 1, KOKKOS_LAMBDA (const int) {
    new (t) T();
  });
  return t;
}
inline void free_on_device(void * t) {
  Kokkos::kokkos_free<MemSpace>(t);
}


TEST_CASE("kokkos_reduce", "virtual_function") {
  BaseClass* base_class = allocate_on_device<DerivedClass>();

  int reduced=0;
  Kokkos::Sum<int> reducer_int(reduced);
  Kokkos::parallel_reduce("CallOnDevice", 3, KOKKOS_LAMBDA(const int, int &value)
  {
    value += base_class->virtual_function();
  },
  reducer_int);
  REQUIRE(reduced==3);
  free_on_device(base_class);
}

