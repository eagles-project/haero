#include "chemUtil.hpp"

namespace chemUtil {

using namespace from_tchem;

chemFiles::chemFiles(std::string chemDir){
  prefixPath = "tests/" + chemDir;
  chemFile = prefixPath + "chem.inp";
  thermFile = prefixPath + "therm.dat";
  outputFile = prefixPath + "omega.dat";
}

// chemSolver::chemSolver(SimulationInput sim_inp)
// {
//   k1 = sim_inp.reactions.products[]
//   parse_problem_specific(SimulationInput sim_inp);
// }

chemSolver::chemSolver(std::string chemDir, bool detail, int inBatch,
                       bool iverbose,
                       Real k1, Real k2,
                       Real initX, Real initX2)
                       :
                       cfiles(chemDir),
                       verbose(iverbose),
                       nBatch(inBatch),
                       policy(TChem::exec_space(), nBatch, Kokkos::AUTO()) {
  // optionally print some configuration info
    // (note: this doesn't appear to do anything when compiled in serial)
  // TChem::exec_space::print_configuration(std::cout, detail);
  // TChem::host_exec_space::print_configuration(std::cout, detail);

  // construct kmd and use the view for testing
  kmd = TChem::KineticModelData(cfiles.chemFile, cfiles.thermFile);
  std::ofstream output(cfiles.thermFile);
  // create a const version
  kmcd = kmd.createConstData<TChem::exec_space>();


  // 2d view for reaction rates (calculated in team_invoke_detail())
  reactRate = real_type_2d_view("ReactionRates", nBatch, kmcd.nSpec);

  // input: state vectors: temperature, pressure and mass fraction
  state = real_type_2d_view("StateVector", nBatch, kmcd.nSpec);
  // output: omega, reaction rates
  omega = real_type_2d_view("NetProductionRate", nBatch, kmcd.nSpec);

  // create mirror views on host to set the initial values of
  // reactRate then deep copy to device
  auto reactRate_host = Kokkos::create_mirror_view(reactRate);

  // assign the constructor arguments to the corresponding mirror view
  for (int i = 0; i < nBatch; ++i)
  {
    reactRate_host(i, 0) = k1;
    reactRate_host(i, 1) = k2;
  }

  // deep copy to device
  Kokkos::deep_copy(reactRate, reactRate_host);

  /// create a mirror view to store input from a file
  auto state_host = Kokkos::create_mirror_view(state);

  // set the initial condition for the state variables (X and X2)
  // this is done on host and then cloned to the mirror view on device
  {
    auto state_host_at_0 = Kokkos::subview(state_host, 0, Kokkos::ALL());
    state_host_at_0(0) = initX;
    state_host_at_0(1) = initX2;

    // FIXME(mjs): I don't know what this does just keeping it around from what I
              // picked out of TChem, itself
    // TChem::Test::cloneView(state_host);
  }

  // reset timer and time the deep copy from host to device
  timer.reset();
  Kokkos::deep_copy(state, state_host);
  t_deepcopy = timer.seconds();

  // this bit is a little over my head
  const ordinal_type level = 1;
  const ordinal_type per_team_extent = getWorkSpaceSize(kmcd);
  const ordinal_type per_team_scratch =
    TChem::Scratch<real_type_1d_view>::shmem_size(per_team_extent);
  policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));

} // end chemSolver constructor

void chemSolver::print_summary(const chemFiles& cfiles){
  printf("---------------------------------------------------\n");
  printf("Testing Arguments: \n batch size %d\n chemfile %s\n thermfile %s\n outputfile %s\n verbose %s\n",
         nBatch,
         cfiles.chemFile.c_str(),
         cfiles.thermFile.c_str(),
         // inputFile.c_str(),
         cfiles.outputFile.c_str(),
         verbose ? "true" : "false");
  printf("---------------------------------------------------\n");
  printf("Time deep copy      %e [sec] %e [sec/sample]\n",
         t_deepcopy,
         t_deepcopy / Real(nBatch));
  printf("Time reaction rates %e [sec] %e [sec/sample]\n",
         t_device_batch,
         t_device_batch / Real(nBatch));
} // end chemSolver::print_summary

real_type_2d_view chemSolver::get_results(){
  // reset timer
  timer.reset();
  // run the model
  from_tchem::SourceTermToyProblem::runDeviceBatch(policy, reactRate, state, omega, kmcd);
  Kokkos::fence(); /// timing purpose
  t_device_batch = timer.seconds();
  /// create a mirror view of omega (output) to export a file
  if (verbose) {
    print_summary(cfiles);
    auto omega_host = Kokkos::create_mirror_view(omega);
    Kokkos::deep_copy(omega_host, omega);

    /// all values are same (print only the first one)
    {
      auto omega_host_at_0 = Kokkos::subview(omega_host, 0, Kokkos::ALL());
      TChem::Test::writeReactionRates(
        cfiles.outputFile, kmcd.nSpec, omega_host_at_0);
    }
  }
  return omega;
} // end chemSolver::get_results

} // end chemUtil namespace

/*
***NOTE: everything in this namespace is copied directly from TChem***
*/
namespace from_tchem {

  using Real = haero::Real;

  template<typename KineticModelConstDataType>
  KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize(
    const KineticModelConstDataType& kmcd)
  {
    return 6 * kmcd.nReac;
  }

    template<typename MemberType,
             typename RealType1DViewType,
             typename OrdinalType1DViewType,
             typename KineticModelConstDataType>
    KOKKOS_INLINE_FUNCTION static void team_invoke_detail(
      const MemberType& member,
      /// input
      const RealType1DViewType& reactRate,
      const RealType1DViewType& concX,
      /// output
      const RealType1DViewType& omega, /// (kmcd.nSpec)
      const RealType1DViewType& kfor,
      const RealType1DViewType& krev,
      const RealType1DViewType& ropFor,
      const RealType1DViewType& ropRev,
      const OrdinalType1DViewType& iter,
      /// const input from kinetic model
      const KineticModelConstDataType& kmcd)
    {

      kfor(0) = reactRate(0);
      kfor(1) = reactRate(1);
      krev(0) = 0;
      krev(1) = 0;

      member.team_barrier();


      /// 2. compute rate-of-progress
      TChem::Impl::RateOfProgress::team_invoke(member,
                                  kfor,
                                  krev,
                                  concX, /// input
                                  ropFor,
                                  ropRev, /// output
                                  iter,   /// workspace for iterators
                                  kmcd);

      member.team_barrier();

      /// 3. assemble reaction rates
      auto rop = ropFor;
      Kokkos::parallel_for(
        Kokkos::TeamVectorRange(member, kmcd.nReac), [&](const ordinal_type& i) {
          rop(i) -= ropRev(i);
          const Real rop_at_i = rop(i);
          for (ordinal_type j = 0; j < kmcd.reacNreac(i); ++j) {
            const ordinal_type kspec = kmcd.reacSidx(i, j);
            // omega(kspec) += kmcd.reacNuki(i,j)*rop_at_i;
            const Real val = kmcd.reacNuki(i, j) * rop_at_i;
            Kokkos::atomic_fetch_add(&omega(kspec), val);
          }
          const ordinal_type joff = kmcd.reacSidx.extent(1) / 2;
          for (ordinal_type j = 0; j < kmcd.reacNprod(i); ++j) {
            const ordinal_type kspec = kmcd.reacSidx(i, j + joff);
            // omega(kspec) += kmcd.reacNuki(i,j+joff)*rop_at_i;
            const Real val = kmcd.reacNuki(i, j + joff) * rop_at_i;
            Kokkos::atomic_fetch_add(&omega(kspec), val);
          }
        });


      member.team_barrier();

  #if defined(TCHEM_ENABLE_SERIAL_TEST_OUTPUT) && !defined(__CUDA_ARCH__)
      if (member.league_rank() == 0) {
        FILE* fs = fopen("SourceTermToyProblem.team_invoke.test.out", "a+");
        fprintf(fs, ":: SourceTermToyProblem::team_invoke\n");
        fprintf(fs, ":::: input\n");
        fprintf(fs,
                "     nSpec %3d, nReac %3d\n",
                kmcd.nSpec,
                kmcd.nReac);
        //
        fprintf(fs, ":: concX\n");
        for (int i = 0; i < int(concX.extent(0)); ++i)
          fprintf(fs,
                  "     i %3d, kfor %e\n",
                  i,
                  concX(i));
        fprintf(fs, ":: kfor\n");
        for (int i = 0; i < int(kfor.extent(0)); ++i)
          fprintf(fs,
                  "     i %3d, kfor %e\n",
                  i,
                  kfor(i));
        fprintf(fs, ":: krev\n");
        for (int i = 0; i < int(krev.extent(0)); ++i)
          fprintf(fs,
                  "     i %3d, krev %e\n",
                  i,
                  krev(i));
        //
        fprintf(fs, ":: ropFor\n");
        for (int i = 0; i < int(ropFor.extent(0)); ++i)
          fprintf(fs,
                  "     i %3d, kfor %e\n",
                  i,
                  ropFor(i));
        fprintf(fs, ":: ropRev\n");
        for (int i = 0; i < int(ropRev.extent(0)); ++i)
          fprintf(fs,
                  "     i %3d, krev %e\n",
                  i,
                  ropRev(i));
        fprintf(fs, "\n");
        fprintf(fs, ":::: output\n");
        for (int i = 0; i < int(omega.extent(0)); ++i)
          fprintf(fs, "     i %3d, omega %e\n", i, omega(i));
        fprintf(fs, "\n");
      }
  #endif

      // if (member.league_rank() == 0) {
      //   FILE *fs = fopen("SourceTermToyProblem.team_invoke.test.out", "a+");
      //   for (int i=0;i<int(Crnd.extent(0));++i)
      //     fprintf(fs, " %d %e\n", i , Real(1e-3)*Crnd(i)*kfor(i));
      // }
    }


    template<typename MemberType,
             typename WorkViewType,
             typename RealType1DViewType,
             typename KineticModelConstDataType>
    KOKKOS_FORCEINLINE_FUNCTION static void team_invoke(
      const MemberType& member,
      /// input
      const RealType1DViewType& reactRate,
      const RealType1DViewType& X, /// (kmcd.nSpec)
      /// output
      const RealType1DViewType& omega, /// (kmcd.nSpec)
      /// workspace
      const WorkViewType& work,
      /// const input from kinetic model
      const KineticModelConstDataType& kmcd)
    {

      ///
      auto w = (Real*)work.data();
      auto kfor = RealType1DViewType(w, kmcd.nReac);
      w += kmcd.nReac;
      auto krev = RealType1DViewType(w, kmcd.nReac);
      w += kmcd.nReac;
      auto ropFor = RealType1DViewType(w, kmcd.nReac);
      w += kmcd.nReac;
      auto ropRev = RealType1DViewType(w, kmcd.nReac);
      w += kmcd.nReac;
      auto iter = Kokkos::View<ordinal_type*,
                               Kokkos::LayoutRight,
                               typename WorkViewType::memory_space>(
        (ordinal_type*)w, kmcd.nReac * 2);
      w += kmcd.nReac * 2;

      team_invoke_detail(member,
                         reactRate,
                         X,
                         omega,
                         kfor,
                         krev,
                         ropFor,
                         ropRev,
                         iter,
                         kmcd);
    }

  template<typename PolicyType,
          // typename RealType1DViewType,
          typename RealType2DViewType,
          typename KineticModelConstType>
  //
  void SourceTermToyProblem_TemplateRun( /// input
    const std::string& profile_name,
    /// team size setting
    const PolicyType& policy,
    const RealType2DViewType& reactRate,
    const RealType2DViewType& state,
    const RealType2DViewType& SourceTermToyProblem,
    const KineticModelConstType& kmcd
  )
  {
    Kokkos::Profiling::pushRegion(profile_name);
    using policy_type = PolicyType;

    const ordinal_type level = 1;
    const ordinal_type per_team_extent = getWorkSpaceSize(kmcd);

    Kokkos::parallel_for(
      profile_name,
      policy,
      KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
        const ordinal_type i = member.league_rank();
        const real_type_1d_view state_at_i =
          Kokkos::subview(state, i, Kokkos::ALL());
        const real_type_1d_view reactRate_at_i =
          Kokkos::subview(reactRate, i, Kokkos::ALL());
        const real_type_1d_view SourceTermToyProblem_at_i =
          Kokkos::subview(SourceTermToyProblem, i, Kokkos::ALL());

        Scratch<real_type_1d_view> work(member.team_scratch(level),
                                        per_team_extent);

        team_invoke(member, reactRate_at_i, state_at_i,
                    SourceTermToyProblem_at_i, work, kmcd);

      });
    Kokkos::Profiling::popRegion();
  }

  void SourceTermToyProblem::runDeviceBatch( /// input
    typename UseThisTeamPolicy<exec_space>::type& policy,
    const real_type_2d_view& reactRate,
    const real_type_2d_view& state,
    /// output
    const real_type_2d_view& SourceTermToyProblem,
    /// const data from kinetic model
    const KineticModelConstDataDevice& kmcd)
  {

    SourceTermToyProblem_TemplateRun( /// input
      "SourceTermToyProblem::runDeviceBatch",
      /// team size setting
      policy,
      reactRate,
      state,
      SourceTermToyProblem,
      kmcd);

  }

} // end from_tchem namespace
