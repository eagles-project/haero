#include "haero/aerosol_process.hpp"

#include "catch2/catch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <cmath>
#include <iostream>

using namespace haero;

class MyAerosolProcess final : public DeviceAerosolProcess<MyAerosolProcess> {
 public:
  MyAerosolProcess(const std::string &name, const int num_lev,
                   const Diagnostics::Token aer_0,
                   const Diagnostics::Token aer_1,
                   const Diagnostics::Token gen_0)
      : DeviceAerosolProcess<MyAerosolProcess>(name),
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
    const int num_vert_packs = temp.extent(0);
    const int nk = HAERO_PACK_SIZE * num_vert_packs;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nk), [=](int k) {
      generic_var(pack_info::pack_idx(k))[pack_info::vec_idx(k)] = 0;
      for (int p = 0; p < num_aerosol_populations; ++p) {
        aero_tend(p, pack_info::pack_idx(k))[pack_info::vec_idx(k)] = 0;
        first_aersol(p, pack_info::pack_idx(k))[pack_info::vec_idx(k)] = 0;
      }
    });
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nk), [=](int k) {
      for (int n = 0; n < num_levels; ++n) {
        generic_var(pack_info::pack_idx(k))[pack_info::vec_idx(k)] +=
            temp(pack_info::pack_idx(k))[pack_info::vec_idx(k)];
      }
    });

    for (int n = 0; n < num_levels; ++n) {
      for (int p = 0; p < num_aerosol_populations; ++p) {
        Real reduced = 0;
        Kokkos::Sum<Real> reducer_real(reduced);
        Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team, nk),
            [=](int k, Real &v) {
              v += int_aerosols(0,
                                pack_info::pack_idx(n))[pack_info::vec_idx(n)] *
                   temp(pack_info::pack_idx(k))[pack_info::vec_idx(k)];
            },
            reducer_real);
        aero_tend(p, pack_info::pack_idx(n))[pack_info::vec_idx(n)] = reduced;
      }
    };
    for (int n = 0; n < num_levels; ++n) {
      for (int p = 0; p < num_populations; ++p) {
        Real reduced = 0;
        Kokkos::Sum<Real> reducer_real(reduced);
        Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team, nk),
            [=](int k, Real &v) {
              v += n * p * temp(pack_info::pack_idx(k))[pack_info::vec_idx(k)];
            },
            reducer_real);
        first_aersol(p, pack_info::pack_idx(n))[pack_info::vec_idx(n)] =
            reduced;
        second_aersol(p, pack_info::pack_idx(n))[pack_info::vec_idx(n)] = p * n;
      }
    }
  };

 private:
  const int num_levels;
  const Diagnostics::Token aersol_0;
  const Diagnostics::Token aersol_1;
  const Diagnostics::Token generic_0;
};

TEST_CASE("aerosol_process_multi_column_run", "") {
  // Define an arbitrary number of columns for this rank to
  // process.  This would be out of the total number of
  // columns in the atmosphere from the whole-earth model.
  const int num_atm_columns = 1245;
  using View1D = kokkos_device_type::view_1d<PackType>;
  using View2D = kokkos_device_type::view_2d<PackType>;
  using MultiColumnSpeciesView = kokkos_device_type::view_3d<PackType>;
  using MultiColumnModeView = DeviceType::view_3d<PackType>;

  // Create a default aerosol configuration for testing
  // purposes only.
  auto aero_config = ModalAerosolConfig::create_mam4_config();
  int num_levels = 72;  // number of levels
  int num_vert_packs = PackInfo::num_packs(num_levels);
  int num_iface_packs = PackInfo::num_packs(num_levels+1);
  int num_gases = aero_config.num_gases();
  int num_modes = aero_config.num_modes();

  // Set up some prognostics aerosol data views
  const int num_aero_populations = aero_config.num_aerosol_populations;

  // Define a pseudo-random number generator [0,1) based on a few primes.
  static constexpr unsigned p0{987659};
  static constexpr unsigned p1{12373};
  long unsigned seed{54319};
  auto random = [&]() {
    seed = (seed * p1) % p0;
    return Real(seed) / p0;
  };

  // interstitial aerosols mmr [kg/kg(of air)]
  MultiColumnSpeciesView mult_col_int_aerosols(
      "multi column interstitial aerosols", num_atm_columns,
      num_aero_populations, num_vert_packs);
  // cloud borne aerosols mmr [kg/kg(of air)]
  MultiColumnSpeciesView mult_col_cld_aerosols(
      "multi column cloudborne aerosols", num_atm_columns, num_aero_populations,
      num_vert_packs);
  {
    // Set initial conditions
    // aerosols mass mixing ratios
    auto h_int_aerosols = Kokkos::create_mirror_view(mult_col_int_aerosols);
    for (int column = 0; column < num_atm_columns; ++column) {
      for (std::size_t i = 0; i < num_levels; ++i) {
        h_int_aerosols(column, 0,
                       pack_info::pack_idx(i))[pack_info::vec_idx(i)] = i;
      }
    }
    Kokkos::deep_copy(mult_col_int_aerosols, h_int_aerosols);
  }
  MultiColumnSpeciesView mult_col_gases("multi column gases", num_atm_columns,
                                        num_gases, num_vert_packs);

  MultiColumnModeView mult_col_int_num_nmrs(
      "multi column interstitial number nmrs", num_atm_columns,
      num_modes,        // interstitial aerosols number
      num_vert_packs);  // mixing ratios [#/kg(of air)]
  MultiColumnModeView mult_col_cld_num_nmrs(
      "multi column cloud borne number nmrs", num_atm_columns,
      num_modes,        // cloud borne aerosols number
      num_vert_packs);  // mixing ratios [#/kg(of air)]
  {
    // aerosols number mixing ratios
    auto h_int_nmrs = Kokkos::create_mirror_view(mult_col_int_num_nmrs);
    auto h_cld_nmrs = Kokkos::create_mirror_view(mult_col_cld_num_nmrs);
    for (int column = 0; column < num_atm_columns; ++column) {
      for (std::size_t imode = 0; imode < num_modes; ++imode) {
        for (std::size_t k = 0; k < num_vert_packs; ++k) {
          h_int_nmrs(column, imode, k) = 1e8 + random();
          h_cld_nmrs(column, imode, k) = 1e8 + random();
        }
      }
    }
    Kokkos::deep_copy(mult_col_int_num_nmrs, h_int_nmrs);
    Kokkos::deep_copy(mult_col_cld_num_nmrs, h_cld_nmrs);
  }
  // Set up atmospheric data and populate it.
  // Each view is vertical_levels x number_of_columns
  // These 2D Views are then sub-viewed for each
  // column to set up the classes.
  View2D mult_col_temp("temperature", num_atm_columns, num_vert_packs);  //[K]
  {
    auto host_temp = Kokkos::create_mirror_view(mult_col_temp);
    for (int j = 0; j < num_atm_columns; ++j)
      for (int i = 0; i < num_levels; ++i) {
        host_temp(j, pack_info::pack_idx(i))[pack_info::vec_idx(i)] = i;
      }
    Kokkos::deep_copy(mult_col_temp, host_temp);
  }
  View2D mult_col_press("pressure", num_atm_columns, num_vert_packs);  //[Pa]
  View2D mult_col_rel_hum("relative humidity", num_atm_columns, num_vert_packs);
  View2D mult_col_qv("vapor mixing ratio", num_atm_columns, num_vert_packs);
  View2D mult_col_pdel("hydrostatic_dp", num_atm_columns,
                       num_vert_packs);  //[Pa]
  View2D mult_col_ht("height", num_atm_columns, num_iface_packs);

  // Views for the Hearo classes.  These are views of class instances
  // not pointers to instances.  This means that the classes must
  // have default constructors and copy operators to define them.
  // There is one class instance for every column of the atmosphere.
  // These are device views so will have to be deep copied in order
  // to define the class instances.
  std::vector<Diagnostics> host_mult_col_diagnostics(num_atm_columns);
  std::vector<Atmosphere> host_mult_col_atmosphere(num_atm_columns);
  std::vector<Prognostics> host_mult_col_prognostics(num_atm_columns);
  std::vector<Tendencies> host_mult_col_tendencies(num_atm_columns);
  kokkos_device_type::view_1d<Diagnostics> mult_col_diagnostics(
      "multi column diags", num_atm_columns);
  kokkos_device_type::view_1d<Atmosphere> mult_col_atmosphere(
      "multi column Atmosphere", num_atm_columns);
  kokkos_device_type::view_1d<Prognostics> mult_col_prognostics(
      "multi column Prognostics", num_atm_columns);
  kokkos_device_type::view_1d<Tendencies> mult_col_tendencies(
      "multi column Tendencies", num_atm_columns);

  Diagnostics::Token aersol_0 = Diagnostics::VAR_NOT_FOUND;
  Diagnostics::Token aersol_1 = Diagnostics::VAR_NOT_FOUND;
  Diagnostics::Token generic_0 = Diagnostics::VAR_NOT_FOUND;
  // Loop over all of the columns to define the class data.
  for (int column = 0; column < num_atm_columns; ++column) {
    // Kokkos views can be sub-viewed in the constructor using
    // the special ALL flag and selecting a single column.
    //
    // All of these Kokkos views are Device Views!  Do not
    // try to access them on Host without a create_mirrow_view.
    // So, even though host_atmos is an instance of an
    // Atmosphere that is defined on the host, the Views
    // used to define it are Device views.
    SpeciesColumnView int_aerosols(mult_col_int_aerosols, column, Kokkos::ALL,
                                   Kokkos::ALL);
    SpeciesColumnView cld_aerosols(mult_col_cld_aerosols, column, Kokkos::ALL,
                                   Kokkos::ALL);
    SpeciesColumnView gases(mult_col_gases, column, Kokkos::ALL, Kokkos::ALL);
    ModeColumnView int_num_nmrs(mult_col_int_num_nmrs, column, Kokkos::ALL,
                                Kokkos::ALL);
    ModeColumnView cld_num_nmrs(mult_col_cld_num_nmrs, column, Kokkos::ALL,
                                Kokkos::ALL);

    View1D temp(mult_col_temp, column, Kokkos::ALL);    //[K]
    View1D press(mult_col_press, column, Kokkos::ALL);  //[Pa]
    View1D rel_hum(mult_col_rel_hum, column, Kokkos::ALL);
    View1D qv(mult_col_qv, column, Kokkos::ALL);
    View1D pdel(mult_col_pdel, column, Kokkos::ALL);  //[Pa]
    View1D ht(mult_col_ht, column, Kokkos::ALL);
    Real pblh{100.0};  // planetary BL height [m]

    // Create the Hearo class instances from the Kokkos
    // Views.  The diagnostics and tendencies classes
    // define and hold their own Kokkos Views.  In an acutal
    // application these two classes would be checked on
    // the host to make sure the correct Views are defined.
    Prognostics progs(aero_config, num_levels, int_aerosols, cld_aerosols,
        int_num_nmrs, cld_num_nmrs, gases);
    HostDiagnostics diags(aero_config, num_levels);

    const auto t0 = diags.create_aerosol_var("First Aerosol");
    const auto t1 = diags.create_aerosol_var("Second Aerosol");
    const auto g0 = diags.create_var("Generic Aerosol");
    EKAT_REQUIRE(Diagnostics::VAR_NOT_FOUND == aersol_0 || t0 == aersol_0);
    EKAT_REQUIRE(Diagnostics::VAR_NOT_FOUND == aersol_1 || t1 == aersol_1);
    EKAT_REQUIRE(Diagnostics::VAR_NOT_FOUND == generic_0 || g0 == generic_0);
    aersol_0 = t0;
    aersol_1 = t1;
    generic_0 = g0;

    const Tendencies tends(progs);
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
    const Atmosphere atm(num_levels, temp, press, qv, ht, pdel, pblh);
    host_mult_col_atmosphere[column] = atm;
    host_mult_col_tendencies[column] = tends;
    host_mult_col_prognostics[column] = progs;
    host_mult_col_diagnostics[column] = diags;
    // Create the device copies of diagnostics, atmosphere, prognostics and
    // tendencies
    // These are created on host and copied to device through a
    // lambda copy.  This works for POD and Kokkos::Views but should
    // not be used to copy things like std::vector.  An example is the
    // HostDiagnostics class which has std:: member data so has to be cast as
    // a Diagnostics instance to copy to device.
    const Diagnostics &dev_diags = diags;
    Kokkos::parallel_for(
        "Put into device storage", 1, KOKKOS_LAMBDA(const int) {
          mult_col_atmosphere(column) = atm;
          mult_col_tendencies(column) = tends;
          mult_col_prognostics(column) = progs;
          mult_col_diagnostics(column) = dev_diags;
        });
  }

  SECTION("my_aerosol_process_multi_column_run") {
    // Create and initialize our process.
    const std::string name = "CloudProcess";
    auto process = new MyAerosolProcess(name, num_levels, aersol_0,
                                        aersol_1, generic_0);
    process->init(aero_config);
    auto d_process = process->copy_to_device();
    typedef ekat::ExeSpaceUtils<>::TeamPolicy::member_type TeamHandleType;
    const auto &teamPolicy = ekat::ExeSpaceUtils<>::get_default_team_policy(
        num_atm_columns, num_vert_packs);
    Real t = 0.0, dt = 30.0;
    // Parallel launch on GPU or over threads.
    Kokkos::parallel_for(
        teamPolicy, KOKKOS_LAMBDA(const TeamHandleType &team) {
          const int column = team.league_rank();

          // Since the Kokkos
          Diagnostics &diags = mult_col_diagnostics(column);
          Atmosphere &atmos = mult_col_atmosphere(column);
          Prognostics &progs = mult_col_prognostics(column);
          Tendencies &tends = mult_col_tendencies(column);
          d_process->run(team, t, dt, progs, atmos, diags, tends);
        });
    AerosolProcess::delete_on_device(d_process);
  }

  for (int column = 0; column < num_atm_columns; ++column) {
    const Atmosphere &atm = host_mult_col_atmosphere[column];
    // To view the data in a column, copy the View back to host
    auto host_temp = Kokkos::create_mirror_view(atm.temperature);
    Kokkos::deep_copy(host_temp, atm.temperature);
    for (int lev = 0; lev < num_vert_packs; ++lev) {
      for (int s = 0; s < HAERO_PACK_SIZE; ++s) {
        const Real temp = host_temp(lev)[s];
        REQUIRE(temp == lev * HAERO_PACK_SIZE + s);
      }
    }

    using fp_helper = FloatingPoint<float>;
    const Diagnostics &diags = host_mult_col_diagnostics[column];
    const Tendencies &tends = host_mult_col_tendencies[column];
    const SpeciesColumnView first_aersol = diags.aerosol_var(aersol_0);
    const SpeciesColumnView second_aersol = diags.aerosol_var(aersol_1);
    const ColumnView generic_var = diags.var(generic_0);
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
    for (int k = 0; k < num_levels; ++k) {
      const Real val = num_levels * k;
      const int pack_idx = pack_info::pack_idx(k);
      const int vec_idx = pack_info::vec_idx(k);
      const Real tst = host_generic_var(pack_idx)[vec_idx];
      std::cout << __LINE__ << " " << tst << " " << val << std::endl;
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
    for (int i = 0; i < num_levels; ++i) {
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
  // All the class instances have been stored in std::vector and should be
  // deleted when the views go out of scope.
}
