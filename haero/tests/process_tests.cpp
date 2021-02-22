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
  MyPrognosticProcess(ProcessType type, 
                      const std::string& name,
                      const Diagnostics::TOKEN aer_0,
                      const Diagnostics::TOKEN aer_1,
                      const Diagnostics::TOKEN gen_0) : 
   PrognosticProcess(type, name),
   aersol_0 (aer_0),
   aersol_1 (aer_1),
   generic_0 (gen_0)
 {}

  KOKKOS_INLINE_FUNCTION
  virtual ~MyPrognosticProcess() {}

  KOKKOS_INLINE_FUNCTION
  MyPrognosticProcess(const MyPrognosticProcess& pp) : 
   PrognosticProcess(pp),
   aersol_0 (pp.aersol_0),
   aersol_1 (pp.aersol_1),
   generic_0 (pp.generic_0)
 {}

  KOKKOS_FUNCTION
  virtual void run(const Model& model,
                   Real t, Real dt,
                   const Prognostics& prognostics,
                   const Atmosphere& atmosphere,
                   const Diagnostics& diagnostics,
                   Tendencies& tendencies) const {

    const Prognostics::SpeciesColumnView int_aerosols = prognostics.interstitial_aerosols();
    const Atmosphere::ColumnView temp = atmosphere.temperature();
    const Diagnostics::SpeciesColumnView first_aersol = diagnostics.aerosol_var(aersol_0);
    const Diagnostics::SpeciesColumnView second_aersol = diagnostics.aerosol_var(aersol_1);
    const Diagnostics::ColumnView        generic_var   = diagnostics.var(generic_0);
    Tendencies::SpeciesColumnView aero_tend = tendencies.interstitial_aerosols();

    const int num_vert_packs  = int_aerosols.extent(1);
    const int num_levels      = temp.extent(0);
    const int num_populations = first_aersol.extent(0);
    const int num_aerosol_populations = aero_tend.extent(0);
    // int_aerosols 1 x num_vert_packs
    // temp  num_levels
    // first_aersol num_populations x num_vert_packs
    // second_aersol num_populations x num_vert_packs
    // generic_var num_levels
    // aero_tend num_aerosol_populations x num_vert_packs 
    for (int k=0; k<num_levels; ++k) generic_var(k) = 0;
    for (int i=0; i<num_vert_packs; ++i) {
      for (int j=0; j<num_aerosol_populations; ++j) {
        for (int k=0; k<num_levels; ++k) {
          aero_tend(j,i) += int_aerosols(0,i) * temp(k);
        }
      }
      for (int k=0; k<num_levels; ++k) {
        generic_var(k) += temp(k);
      }
     for (int j=0; j<num_populations; ++j) {
        first_aersol(j,i) = 0;
        for (int k=0; k<num_levels; ++k) {
          first_aersol(j,i) += temp(k) * aero_tend(j,i);
        }
        second_aersol(j,i) = aero_tend(j,i);
      }
    }
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
private:
  const Diagnostics::TOKEN aersol_0;
  const Diagnostics::TOKEN aersol_1;
  const Diagnostics::TOKEN generic_0;
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
  {
    auto host_temp  =  Kokkos::create_mirror_view(temp);
    for (int i=0; i<num_levels; ++i) {
      host_temp(i) = 100*i;
    }
    Kokkos::deep_copy(temp, host_temp);
  }
  Atmosphere atmos(num_levels, temp, press, rel_hum, ht);

  std::vector<int> num_aero_species(num_modes);
  std::vector<Mode> modes = create_mam4_modes();
  std::map<std::string, std::vector<std::string>> mode_species = create_mam4_mode_species();
  for (int m = 0; m < num_modes; ++m) {
    num_aero_species[m] = mode_species[modes[m].name].size();
  }

  Diagnostics diags(num_modes, num_aero_species, num_gases, num_levels);
  const Diagnostics::TOKEN aersol_0 = diags.create_aerosol_var("First Aerosol");
  const Diagnostics::TOKEN aersol_1 = diags.create_aerosol_var("Second Aerosol");
  const Diagnostics::TOKEN generic_0 = diags.create_var("Generic Aerosol");

  Tendencies tends(progs);
  {
    const int num_populations = progs.num_aerosol_populations();
    Tendencies::SpeciesColumnView aero_tend = tends.interstitial_aerosols();
    auto host_aero_tend  =  Kokkos::create_mirror_view(aero_tend);
    for (int i=0; i<num_vert_packs; ++i) {
      for (int j=0; j<num_populations; ++j) {
        host_aero_tend(j,i) = .1*i + .1*j;
      }
    }
    Kokkos::deep_copy(aero_tend, host_aero_tend);
  }

  ProcessType type = CloudBorneWetRemovalProcess;
  const std::string name = "CloudPrognosticProcess";
  MyPrognosticProcess pp(type, name, aersol_0, aersol_1, generic_0);
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

