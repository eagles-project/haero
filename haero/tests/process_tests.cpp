#include <cstdio>

#include "catch2/catch.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "haero/aerosol_process.hpp"
#include "haero/atmosphere.hpp"
#include "haero/diagnostics.hpp"
#include "haero/floating_point.hpp"
#include "haero/model.hpp"
#include "haero/prognostics.hpp"
#include "haero/view_pack_helpers.hpp"
#include "kokkos/Kokkos_Core.hpp"

using namespace haero;

class MyAerosolProcess final : public DeviceAerosolProcess<MyAerosolProcess> {
 public:
  MyAerosolProcess(AerosolProcessType type, const std::string &name,
                   const int num_lev, const Diagnostics::Token aer_0,
                   const Diagnostics::Token aer_1,
                   const Diagnostics::Token gen_0)
      : DeviceAerosolProcess<MyAerosolProcess>(type, name),
        num_levels(num_lev),
        aersol_0(aer_0),
        aersol_1(aer_1),
        generic_0(gen_0) {}

  KOKKOS_INLINE_FUNCTION
  virtual ~MyAerosolProcess() {}

  KOKKOS_INLINE_FUNCTION
  MyAerosolProcess(const MyAerosolProcess &pp)
      : DeviceAerosolProcess<MyAerosolProcess>(pp),
        num_levels(pp.num_levels),
        aersol_0(pp.aersol_0),
        aersol_1(pp.aersol_1),
        generic_0(pp.generic_0) {}

 protected:
  //------------------------------------------------------------------------
  //                                Overrides
  //------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void run_(const TeamType &team, Real t, Real dt,
            const Prognostics &prognostics, const Atmosphere &atmosphere,
            const Diagnostics &diagnostics,
            const Tendencies &tendencies) const {
    const SpeciesColumnView int_aerosols = prognostics.interstitial_aerosols;
    const ColumnView temp = atmosphere.temperature;
    const SpeciesColumnView first_aersol = diagnostics.aerosol_var(aersol_0);
    const SpeciesColumnView second_aersol = diagnostics.aerosol_var(aersol_1);
    const ColumnView generic_var = diagnostics.var(generic_0);
    SpeciesColumnView aero_tend = tendencies.interstitial_aerosols;

    const int num_populations = first_aersol.extent(0);
    const int num_aerosol_populations = aero_tend.extent(0);
    for (int k = 0; k < num_levels; ++k) {
      generic_var(pack_info::pack_idx(k))[pack_info::vec_idx(k)] = 0;
    }
    for (int i = 0; i < num_levels; ++i) {
      for (int k = 0; k < num_levels; ++k) {
        generic_var(pack_info::pack_idx(k))[pack_info::vec_idx(k)] +=
            temp(pack_info::pack_idx(k))[pack_info::vec_idx(k)];
      }
      for (int j = 0; j < num_aerosol_populations; ++j) {
        aero_tend(j, pack_info::pack_idx(i))[pack_info::vec_idx(i)] = 0;
        for (int k = 0; k < num_levels; ++k) {
          aero_tend(j, pack_info::pack_idx(i))[pack_info::vec_idx(i)] +=
              int_aerosols(0, pack_info::pack_idx(i))[pack_info::vec_idx(i)] *
              temp(pack_info::pack_idx(k))[pack_info::vec_idx(k)];
        }
      }
      for (int j = 0; j < num_populations; ++j) {
        first_aersol(j, pack_info::pack_idx(i))[pack_info::vec_idx(i)] = 0;
        for (int k = 0; k < num_levels; ++k) {
          first_aersol(j, pack_info::pack_idx(i))[pack_info::vec_idx(i)] +=
              i * j * temp(pack_info::pack_idx(k))[pack_info::vec_idx(k)];
        }
        second_aersol(j, pack_info::pack_idx(i))[pack_info::vec_idx(i)] = j * i;
      }
    }
  };

 private:
  const int num_levels;
  const Diagnostics::Token aersol_0;
  const Diagnostics::Token aersol_1;
  const Diagnostics::Token generic_0;
};

TEST_CASE("process_tests", "aerosol_process") {
  const int num_levels = 72;
  const int num_gases = 1;
  const int num_modes = 1;
  const Real t = 2.3;
  const Real dt = 0.15;
  int num_vert_packs = num_levels / HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels) {
    num_vert_packs++;
  }

  SpeciesColumnView dev_gases;
  {
    // example of filling a device view from a std::vector
    std::vector<std::vector<Real>> host_gases(num_gases,
                                              std::vector<Real>(num_levels));
    for (int i = 0; i < num_levels; ++i) {
      for (int j = 0; j < num_gases; ++j) {
        host_gases[j][i] = i + j;
      }
    }
    dev_gases = vectors_to_row_packed_2dview(host_gases, "gases");
  }

  SpeciesColumnView dev_int_aerosols("interstitial aerosols", 1,
                                     num_vert_packs);
  {
    // example of filling a device view from a host view
    auto host_int_aerosols = Kokkos::create_mirror_view(dev_int_aerosols);
    for (int i = 0; i < num_levels; ++i) {
      host_int_aerosols(0, pack_info::pack_idx(i))[pack_info::vec_idx(i)] = i;
    }
    Kokkos::deep_copy(dev_int_aerosols, host_int_aerosols);
  }

  SpeciesColumnView dev_cld_aerosols("cloudborne aerosols", 1, num_vert_packs);
  ModeColumnView dev_int_num_mix_ratios("interstitial number mix ratios",
                                        num_modes, num_vert_packs);
  ModeColumnView dev_cld_num_mix_ratios("cloud borne number mix ratios",
                                        num_modes, num_vert_packs);

  auto host_cld_aerosols = Kokkos::create_mirror_view(dev_cld_aerosols);
  auto host_int_num_mix_ratios =
      Kokkos::create_mirror_view(dev_int_num_mix_ratios);
  auto host_cld_num_mix_ratios =
      Kokkos::create_mirror_view(dev_cld_num_mix_ratios);

  for (int i = 0; i < num_levels; ++i) {
    host_cld_aerosols(0, pack_info::pack_idx(i))[pack_info::vec_idx(i)] = i;
    host_int_num_mix_ratios(0, pack_info::pack_idx(i))[pack_info::vec_idx(i)] =
        i;
    host_cld_num_mix_ratios(0, pack_info::pack_idx(i))[pack_info::vec_idx(i)] =
        i;
  }
  Kokkos::deep_copy(dev_cld_aerosols, host_cld_aerosols);
  Kokkos::deep_copy(dev_int_num_mix_ratios, host_int_num_mix_ratios);
  Kokkos::deep_copy(dev_cld_num_mix_ratios, host_cld_num_mix_ratios);

  Prognostics progs(num_modes, {1}, num_gases, num_levels, dev_int_aerosols,
                    dev_cld_aerosols, dev_int_num_mix_ratios,
                    dev_cld_num_mix_ratios, dev_gases);

  Kokkos::View<PackType *> temp("temperature", num_vert_packs);
  Kokkos::View<PackType *> press("pressure", num_vert_packs);
  Kokkos::View<PackType *> qv("vapor mixing ratio", num_vert_packs);
  Kokkos::View<PackType *> pdel("hydrostatic_dp", num_vert_packs);
  int num_iface_packs = (num_levels + 1) / HAERO_PACK_SIZE;
  if (num_iface_packs * HAERO_PACK_SIZE < (num_levels + 1)) {
    num_iface_packs++;
  }
  Kokkos::View<PackType *> ht("height", num_iface_packs);
  {
    auto host_temp = Kokkos::create_mirror_view(temp);
    for (int i = 0; i < num_levels; ++i) {
      host_temp(pack_info::pack_idx(i))[pack_info::vec_idx(i)] = i;
    }
    Kokkos::deep_copy(temp, host_temp);
  }
  Real pblh = 100.0;
  Atmosphere atmos(num_levels, temp, press, qv, ht, pdel, pblh);

  std::vector<int> num_aero_species(num_modes);
  std::vector<Mode> modes = create_mam4_modes();
  std::map<std::string, std::vector<std::string>> mode_species =
      create_mam4_mode_species();
  for (int m = 0; m < num_modes; ++m) {
    num_aero_species[m] = mode_species[modes[m].name()].size();
  }

  HostDiagnostics diagnostics_register(num_modes, num_aero_species, num_gases,
                                       num_levels);
  auto aersol_0 = diagnostics_register.create_aerosol_var("First Aerosol");
  auto aersol_1 = diagnostics_register.create_aerosol_var("Second Aerosol");
  auto generic_0 = diagnostics_register.create_var("Generic Aerosol");

  Tendencies tends(progs);
  {
    const int num_populations = progs.num_aerosol_populations();
    SpeciesColumnView aero_tend = tends.interstitial_aerosols;
    auto host_aero_tend = Kokkos::create_mirror_view(aero_tend);
    for (int i = 0; i < num_levels; ++i) {
      for (int j = 0; j < num_populations; ++j) {
        host_aero_tend(j, pack_info::pack_idx(i))[pack_info::vec_idx(i)] =
            .1 * i + .1 * j;
      }
    }
    Kokkos::deep_copy(aero_tend, host_aero_tend);
  }

  // Create and initialize our process.
  AerosolProcessType type = CloudBorneWetRemovalProcess;
  const std::string name = "CloudProcess";
  MyAerosolProcess pp(type, name, num_levels, aersol_0, aersol_1, generic_0);
  std::vector<AerosolSpecies> aero_species = create_mam4_aerosol_species();
  std::vector<GasSpecies> gas_species = create_mam4_gas_species();
  ModalAerosolConfig aero_config(modes, aero_species, mode_species,
                                 gas_species);
  pp.init(aero_config);

  // Move the process to the device and run it.
  auto device_pp = pp.copy_to_device();
  const Diagnostics &diagnostics = diagnostics_register.GetDiagnostics();
  const auto &team_policy = haero::TeamPolicy(1u, Kokkos::AUTO);
  Kokkos::parallel_for(
      team_policy, KOKKOS_LAMBDA(const TeamType &team) {
        // Const cast because everything in lambda is const. Need to
        // google how to fix.
        Tendencies *tendency = const_cast<Tendencies *>(&tends);
        device_pp->run(team, t, dt, progs, atmos, diagnostics, *tendency);
      });
  AerosolProcess::delete_on_device(device_pp);

  {
    using fp_helper = FloatingPoint<float>;
    const SpeciesColumnView int_aerosols = progs.interstitial_aerosols;
    const ColumnView temp = atmos.temperature;
    const SpeciesColumnView first_aersol = diagnostics.aerosol_var(aersol_0);
    const SpeciesColumnView second_aersol = diagnostics.aerosol_var(aersol_1);
    const ColumnView generic_var = diagnostics.var(generic_0);
    SpeciesColumnView aero_tend = tends.interstitial_aerosols;

    auto host_first_aersol = Kokkos::create_mirror_view(first_aersol);
    auto host_second_aersol = Kokkos::create_mirror_view(second_aersol);
    auto host_generic_var = Kokkos::create_mirror_view(generic_var);
    auto host_aero_tend = Kokkos::create_mirror_view(aero_tend);
    Kokkos::deep_copy(host_first_aersol, first_aersol);
    Kokkos::deep_copy(host_second_aersol, second_aersol);
    Kokkos::deep_copy(host_generic_var, generic_var);
    Kokkos::deep_copy(host_aero_tend, aero_tend);
    const int num_populations = first_aersol.extent(0);
    const int num_aerosol_populations = aero_tend.extent(0);
    for (int i = 0; i < num_levels; ++i) {
      for (int k = 0; k < num_levels; ++k) {
        const Real val = num_levels * k;
        const Real tst =
            host_generic_var(pack_info::pack_idx(k))[pack_info::vec_idx(k)];
        REQUIRE(fp_helper::equiv(tst, val));
      }
      for (int j = 0; j < num_aerosol_populations; ++j) {
        for (int i = 0; i < num_levels; ++i) {
          const Real val = 2556 * i;
          const Real tst =
              host_aero_tend(j, pack_info::pack_idx(i))[pack_info::vec_idx(i)];
          REQUIRE(fp_helper::equiv(tst, val));
        }
      }
      for (int j = 0; j < num_populations; ++j) {
        {
          const Real val = 2556 * (i * j);
          const Real tst = host_first_aersol(
              j, pack_info::pack_idx(i))[pack_info::vec_idx(i)];
          REQUIRE(fp_helper::equiv(tst, val));
        }
        {
          const Real val = j * i;
          const Real tst = host_second_aersol(
              j, pack_info::pack_idx(i))[pack_info::vec_idx(i)];
          REQUIRE(fp_helper::equiv(tst, val));
        }
      }
    }
  }
}
