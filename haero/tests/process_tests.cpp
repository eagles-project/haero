#include "catch2/catch.hpp"
#include <cstdio>

#include "haero/model.hpp"
#include "haero/prognostics.hpp"
#include "haero/atmosphere.hpp"
#include "haero/diagnostics.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "haero/view_pack_helpers.hpp"
#include "haero/process.hpp"
#include "kokkos/Kokkos_Core.hpp"

using namespace haero;


// Useful little function to allocate class on device. Probably move to header
// later.
#ifdef KOKKOS_ENABLE_CUDA
typedef Kokkos::CudaSpace MemSpace;
#else
typedef Kokkos::HostSpace MemSpace;
#endif
#if 0
template <class T> inline T *allocate_on_device() {
  const std::string debuggingName(typeid(T).name());
  T *t = static_cast<T *>(
      Kokkos::kokkos_malloc<MemSpace>(debuggingName + "_malloc", sizeof(T)));
  Kokkos::parallel_for(debuggingName + "_format", 1,
                       KOKKOS_LAMBDA(const int) { new (t) T(); });
  return t;
}

template <class T> inline T *allocate_on_device(const T &s) {
  const std::string debuggingName(typeid(T).name());
  T *t = static_cast<T *>(
      Kokkos::kokkos_malloc<MemSpace>(debuggingName + "_malloc", sizeof(T)));
  const T &S = s; // Suck into lambda capture space.
  Kokkos::parallel_for(debuggingName + "_format", 1,
                       KOKKOS_LAMBDA(const int) { new (t) T(S); });
  return t;
}

inline void free_on_device(void *t) { Kokkos::kokkos_free<MemSpace>(t); }
#endif



class MyPrognosticProcess : public PrognosticProcess {
public :
  MyPrognosticProcess(ProcessType type, const std::string& name) : PrognosticProcess(type, name) {}

  KOKKOS_INLINE_FUNCTION
  virtual ~MyPrognosticProcess() {}

  KOKKOS_INLINE_FUNCTION
  MyPrognosticProcess(const MyPrognosticProcess& pp) : PrognosticProcess(pp) {}

  KOKKOS_FUNCTION
  virtual void run(const Model& model,
                   Real t, Real dt,
                   const Prognostics& prognostics,
                   const Atmosphere& atmosphere,
                   const Diagnostics& diagnostics,
                   Tendencies& tendencies) const {
    Prognostics::SpeciesColumnView aero_species = prognostics.interstitial_aerosols();
    Real v = aero_species(0,1)[0];
    printf ("%s:%d aero_species:%lf\n", __FILE__,__LINE__,v); 
  };

  MyPrognosticProcess *copy_to_device() {
    const std::string debuggingName(name());
    MyPrognosticProcess *pp = static_cast<MyPrognosticProcess *>(
      Kokkos::kokkos_malloc<MemSpace>(debuggingName + "_malloc", sizeof(MyPrognosticProcess)));
    const MyPrognosticProcess &this_pp = *this; // Suck into lambda capture space.
    Kokkos::parallel_for(debuggingName + "_format", 1,
                         KOKKOS_LAMBDA(const int) { new (pp) MyPrognosticProcess(this_pp); });
    return pp;
  }
};


TEST_CASE("process_tests", "prognostic_process") {

  
  const int num_levels = 72;
  const int num_gases  = 1;
  const int num_modes  = 1;
  const Real t         = 2.3;
  const Real dt        = 0.15;
  int num_vert_packs = num_levels / HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels) {
    num_vert_packs++;
  }
  using kokkos_device_type = ekat::KokkosTypes<ekat::DefaultDevice>;
  using SpeciesColumnView  = kokkos_device_type::view_2d<PackType>;
  using ModalColumnView    = kokkos_device_type::view_2d<PackType>;

  SpeciesColumnView dev_gases;
  {
    // example of filling a device view from a std::vector
    std::vector<std::vector<Real>> host_gases(num_gases, std::vector<Real>(num_vert_packs));
    for (int i=0; i<num_vert_packs; ++i) {
      for (int j=0; j<num_gases; ++j) {
        host_gases[j][i] = i+j;
      }
    }
    dev_gases = vectors_to_row_packed_2dview(host_gases, "gases");
  }

  SpeciesColumnView dev_int_aerosols("interstitial aerosols", 1, num_vert_packs);
  {
    // example of filling a device view from a host view
    auto host_int_aerosols  =  Kokkos::create_mirror_view(dev_int_aerosols);
    for (int i=0; i<num_vert_packs; ++i) {
      host_int_aerosols(0,i) = i;
    }
    Kokkos::deep_copy(dev_int_aerosols, host_int_aerosols);
  }

  SpeciesColumnView dev_cld_aerosols("cloudborne aerosols",   1, num_vert_packs);
  ModalColumnView   dev_modal_concs("modal number concs", num_modes, num_vert_packs);
  auto host_cld_aerosols  =  Kokkos::create_mirror_view(dev_cld_aerosols);
  auto host_modal_concs   =  Kokkos::create_mirror_view(dev_modal_concs);
  for (int i=0; i<num_vert_packs; ++i) {
    host_cld_aerosols(0,i) = i;
    host_modal_concs(0,i) = i;
  }
  Kokkos::deep_copy(dev_cld_aerosols, host_cld_aerosols);
  Kokkos::deep_copy(dev_modal_concs,  host_modal_concs);

  Prognostics progs(num_modes, {1}, num_gases, num_levels, 
                    dev_int_aerosols, 
                    dev_cld_aerosols, 
                    dev_gases, 
                    dev_modal_concs);

  Kokkos::View<PackType*> temp("temperature", num_levels);
  Kokkos::View<PackType*> press("pressure", num_levels);
  Kokkos::View<PackType*> rel_hum("relative humidity", num_levels);
  Kokkos::View<PackType*> ht("height", num_levels+1);
  Atmosphere atmos(num_levels, temp, press, rel_hum, ht);

  std::vector<int> num_aero_species(num_modes);
  std::vector<Mode> modes = create_mam4_modes();
  std::map<std::string, std::vector<std::string>> mode_species = create_mam4_mode_species();
  for (int m = 0; m < num_modes; ++m) {
    num_aero_species[m] = mode_species[modes[m].name].size();
  }
  Diagnostics diags(num_modes, num_aero_species, num_gases, num_levels);

  Tendencies tends(progs);

  ProcessType type = CloudBorneWetRemovalProcess;
  const std::string name = "CloudPrognosticProcess";
  MyPrognosticProcess pp(type, name);
  MyPrognosticProcess *device_pp = pp.copy_to_device();

  std::vector<Species> aero_species = create_mam4_aerosol_species();
  std::vector<Species> gas_species  = create_mam4_gas_species();
  Model* model = Model::ForUnitTests(modes, aero_species, mode_species, gas_species, num_levels);

  typedef Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type
      TeamHandleType;
  const auto &teamPolicy =
      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>(9u, Kokkos::AUTO);
  Kokkos::parallel_for(teamPolicy, KOKKOS_LAMBDA(const TeamHandleType &team) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, 3u),
                         [&](const int &i) {
     // Const cast because everything in lambda is const. Need to google how to fix.
     Tendencies* tendency = const_cast<Tendencies*>(&tends);
     device_pp->run(*model, t, dt, progs, atmos, diags, *tendency);
    });
  });

  Kokkos::kokkos_free<MemSpace>(device_pp);
}


