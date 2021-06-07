#include "chemDriver.hpp"

#include <cstdarg>
#include <cstdio>

#include "haero/physical_constants.hpp"

namespace haero {
namespace chemDriver {

using namespace from_tchem;

// anonymous namespace to hold this YamlException class
namespace {
/// This exception class stores information about errors encountered in reading
/// data from a YAML file.
class YamlException : public std::exception {
 public:
  /// Constructs an exception containing the given descriptive message.
  YamlException(const std::string& message) : _message(message) {}

  /// Constructs an exception containing the given formatting string and
  /// C-style variadic arguments (a la printf).
  YamlException(const char* fmt, ...);

  const char* what() const throw() { return _message.c_str(); }

 private:
  std::string _message;
};

YamlException::YamlException(const char* fmt, ...) {
  char ss[256];
  va_list args;
  va_start(args, fmt);
  vsnprintf(ss, 255, fmt, args);
  va_end(args);
  _message.assign(ss);
}
}  // namespace

ChemicalSpecies::ChemicalSpecies(std::string mname, Real minitial_value,
                                 std::string munits)
    : name(mname), initial_value(minitial_value), units(munits) {}

EnvironmentalConditions::EnvironmentalConditions(Real T0, std::string T_units,
                                                 Real P0, std::string P_units)
    : initial_temp(T0),
      units_temp(T_units),
      initial_pressure(P0),
      units_pressure(P_units) {}
Reaction::Reaction(std::string mtype_str,
                   std::map<std::string, Real> mreactants,
                   std::map<std::string, Real> mproducts,
                   std::map<std::string, Real> mrate_coefficients)
    : type_str(mtype_str), reactants(mreactants), products(mproducts) {
  // convert reaction type string to lowercase
  transform(type_str.begin(), type_str.end(), type_str.begin(), ::tolower);
  // use the type string to set the enumerated type for the reaction and assign
  // the rate coefficients, either based on input or to default values
  if (type_str.compare("arrhenius") == 0) {
    type = arrhenius;
    rate_coefficients["A"] =
        (mrate_coefficients.find("A") != mrate_coefficients.end())
            ? mrate_coefficients["A"]
            : 1.0;
    rate_coefficients["Ea"] =
        (mrate_coefficients.find("Ea") != mrate_coefficients.end())
            ? mrate_coefficients["Ea"]
            : 0.0;
    rate_coefficients["B"] =
        (mrate_coefficients.find("B") != mrate_coefficients.end())
            ? mrate_coefficients["B"]
            : 0.0;
    rate_coefficients["D"] =
        (mrate_coefficients.find("D") != mrate_coefficients.end())
            ? mrate_coefficients["D"]
            : 300.0;
    rate_coefficients["E"] =
        (mrate_coefficients.find("E") != mrate_coefficients.end())
            ? mrate_coefficients["E"]
            : 0.0;
  } else if (type_str.compare("troe") == 0) {
    type = troe;
    rate_coefficients["k0_A"] =
        (mrate_coefficients.find("k0_A") != mrate_coefficients.end())
            ? mrate_coefficients["k0_A"]
            : 1.0;
    rate_coefficients["k0_B"] =
        (mrate_coefficients.find("k0_B") != mrate_coefficients.end())
            ? mrate_coefficients["k0_B"]
            : 0.0;
    rate_coefficients["k0_C"] =
        (mrate_coefficients.find("k0_C") != mrate_coefficients.end())
            ? mrate_coefficients["k0_C"]
            : 0.0;
    rate_coefficients["kinf_A"] =
        (mrate_coefficients.find("kinf_A") != mrate_coefficients.end())
            ? mrate_coefficients["kinf_A"]
            : 1.0;
    rate_coefficients["kinf_B"] =
        (mrate_coefficients.find("kinf_B") != mrate_coefficients.end())
            ? mrate_coefficients["kinf_B"]
            : 0.0;
    rate_coefficients["kinf_C"] =
        (mrate_coefficients.find("kinf_C") != mrate_coefficients.end())
            ? mrate_coefficients["kinf_C"]
            : 0.0;
    rate_coefficients["Fc"] =
        (mrate_coefficients.find("Fc") != mrate_coefficients.end())
            ? mrate_coefficients["Fc"]
            : 0.6;
    rate_coefficients["N"] =
        (mrate_coefficients.find("N") != mrate_coefficients.end())
            ? mrate_coefficients["N"]
            : 1.0;
  } else {
    std::cout << "ERROR: reaction type currently unsupported."
              << "\n";
  }
}

ChemSolver::ChemSolver(SimulationInput& sim_inp)
    : reactions(sim_inp.reactions) {
  parse_tchem_inputs(sim_inp);

  policy = policy_type(TChem::exec_space(), nBatch, Kokkos::AUTO());

  temperature = sim_inp.env_conditions.initial_temp;
  pressure = sim_inp.env_conditions.initial_pressure;
  units_temp = sim_inp.env_conditions.units_temp;
  units_pressure = sim_inp.env_conditions.units_pressure;

  // set the initial state, based on input
  int nSpecies = sim_inp.species.size();
  state = real_type_2d_view("StateVector", nBatch, nSpecies);
  auto state_host = Kokkos::create_mirror_view(state);
  for (int i = 0; i < nBatch; ++i) {
    for (int j = 0; j < nSpecies; ++j) {
      state_host(i, j) = sim_inp.species[j].initial_value;
    }
  }
  Kokkos::deep_copy(state, state_host);

  // optionally print some configuration info
  // (note: this doesn't appear to do anything when compiled in serial)
  TChem::exec_space::print_configuration(std::cout, detail);
  TChem::host_exec_space::print_configuration(std::cout, detail);

  // construct kmd and use the view for testing
  kmd = TChem::KineticModelData(cfiles.chemFile, cfiles.thermFile);
  std::ofstream output(cfiles.thermFile);
  // create a const version
  kmcd = kmd.createConstData<TChem::exec_space>();

  // initialize omega (tendency output)
  omega = real_type_2d_view("NetProductionRate", nBatch, kmcd.nSpec);

  // FIXME: this bit is a little over my head--should probably consider the
  // implications of what's happening here
  const ordinal_type level = 1;
  const ordinal_type per_team_extent = getWorkSpaceSize(kmcd);
  const ordinal_type per_team_scratch =
      TChem::Scratch<real_type_1d_view>::shmem_size(per_team_extent);
  policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));
}

void ChemSolver::print_summary(const ChemFiles& cfiles) {
  printf("---------------------------------------------------\n");
  printf(
      "Testing Arguments: \n batch size %d\n chemfile %s\n thermfile %s\n "
      "outputfile %s\n verbose %s\n",
      nBatch, cfiles.chemFile.c_str(), cfiles.thermFile.c_str(),
      // inputFile.c_str(),
      cfiles.outputFile.c_str(), verbose ? "true" : "false");
  printf("---------------------------------------------------\n");
  printf("Time reaction rates %e [sec] %e [sec/sample]\n", t_device_batch,
         t_device_batch / Real(nBatch));
}  // end ChemSolver::print_summary

real_type_2d_view ChemSolver::get_results() {
  // reset timer
  timer.reset();

  // calculate the reaction rates
  set_reaction_rates();

  // run the model
  SourceTermToyProblem::runDeviceBatch(policy, kfor, krev, state, omega, kmcd);
  Kokkos::fence();  /// timing purpose
  t_device_batch = timer.seconds();
  /// create a mirror view of omega (output) to export a file
  if (verbose) {
    print_summary(cfiles);
    auto omega_host = Kokkos::create_mirror_view(omega);
    Kokkos::deep_copy(omega_host, omega);

    /// all values are same (print only the first one)
    {
      auto omega_host_at_0 = Kokkos::subview(omega_host, 0, Kokkos::ALL());
      TChem::Test::writeReactionRates(cfiles.outputFile, kmcd.nSpec,
                                      omega_host_at_0);
    }
  }
  return omega;
}  // end ChemSolver::get_results

void ChemSolver::parse_tchem_inputs(SimulationInput& sim_inp) {
  // Try to load the input from the yaml file
  try {
    auto root = YAML::LoadFile(sim_inp.input_file);
    if (root["tchem_inputs"] and root["tchem_inputs"].IsMap()) {
      auto node = root["tchem_inputs"];
      if (not node["detail"]) {
        throw YamlException(
            "problem specific entry has no detail boolean (detail).");
      } else if (not node["nbatch"]) {
        throw YamlException(
            "problem specific entry does not specify number "
            "of batches (nbatch).");
      } else if (not node["verbose"]) {
        throw YamlException(
            "problem specific entry has no verbose boolean (verbose).");
      } else {
        detail = node["detail"].as<bool>();
        nBatch = node["nbatch"].as<int>();
        verbose = node["verbose"].as<bool>();
      }
    } else {
      throw YamlException("No tchem_inputs section was found!");
    }
  } catch (YAML::BadFile& e) {
    throw YamlException(e.what());
  } catch (YAML::ParserException& e) {
    throw YamlException(e.what());
  }
}

void ChemSolver::set_reaction_rates() {
  int nRxn = reactions.size();
  std::vector<Real> rxn_rate(nRxn);

  for (int i = 0; i < nRxn; ++i) {
    switch (reactions[i].type) {
      case arrhenius: {
        Real A = reactions[i].rate_coefficients["A"];
        Real Ea = reactions[i].rate_coefficients["Ea"];
        Real B = reactions[i].rate_coefficients["B"];
        Real D = reactions[i].rate_coefficients["D"];
        Real E = reactions[i].rate_coefficients["E"];
        Real kb = constants::boltzmann;
        rxn_rate[i] = A * exp(-Ea / (kb * temperature)) *
                      pow((temperature / D), B) * (1.0 + E * pressure);
        break;
      }
      case troe: {
        rxn_rate[i] = 2.5;
        break;
      }
    }
  }

  // NOTE: so far, we aren't considering reversible reactions,
  // so we set those rates to zero--this could change
  kfor = real_type_2d_view("ForwardRate", nBatch, reactions.size());
  krev = real_type_2d_view("ReverseRate", nBatch, reactions.size());

  /// create a mirror view to store input
  auto kfor_host = Kokkos::create_mirror_view(kfor);
  auto krev_host = Kokkos::create_mirror_view(krev);

  // assign the constructor arguments to the corresponding mirror view
  for (int i = 0; i < nBatch; ++i) {
    for (int j = 0; j < nRxn; ++j) {
      kfor_host(i, j) = rxn_rate[j];
      krev_host(i, j) = 0;
    }
  }

  // deep copy to device
  Kokkos::deep_copy(kfor, kfor_host);
  Kokkos::deep_copy(krev, krev_host);
}

ChemSolver::~ChemSolver() {
  // delete the temporary chem files
  // NOTE: for now, removing the output file, too
  remove(cfiles.chemFile.data());
  remove(cfiles.thermFile.data());
  remove(cfiles.outputFile.data());
}

}  // namespace chemDriver
}  // namespace haero

/*
***NOTE: everything in this namespace is copied directly from TChem***
*/
namespace from_tchem {

using Real = haero::Real;

template <typename KineticModelConstDataType>
KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize(
    const KineticModelConstDataType& kmcd) {
  return 6 * kmcd.nReac;
}

template <typename MemberType, typename RealType1DViewType,
          typename OrdinalType1DViewType, typename KineticModelConstDataType>
KOKKOS_INLINE_FUNCTION static void team_invoke_detail(
    const MemberType& member,
    /// input
    const RealType1DViewType& concX,
    /// output
    const RealType1DViewType& omega,  /// (kmcd.nSpec)
    const RealType1DViewType& kfor, const RealType1DViewType& krev,
    const RealType1DViewType& ropFor, const RealType1DViewType& ropRev,
    const OrdinalType1DViewType& iter,
    /// const input from kinetic model
    const KineticModelConstDataType& kmcd) {
  member.team_barrier();

  /// 2. compute rate-of-progress
  TChem::Impl::RateOfProgress::team_invoke(member, kfor, krev,
                                           concX,  /// input
                                           ropFor,
                                           ropRev,  /// output
                                           iter,    /// workspace for iterators
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
    fprintf(fs, "     nSpec %3d, nReac %3d\n", kmcd.nSpec, kmcd.nReac);
    //
    fprintf(fs, ":: concX\n");
    for (int i = 0; i < int(concX.extent(0)); ++i)
      fprintf(fs, "     i %3d, kfor %e\n", i, concX(i));
    fprintf(fs, ":: kfor\n");
    for (int i = 0; i < int(kfor.extent(0)); ++i)
      fprintf(fs, "     i %3d, kfor %e\n", i, kfor(i));
    fprintf(fs, ":: krev\n");
    for (int i = 0; i < int(krev.extent(0)); ++i)
      fprintf(fs, "     i %3d, krev %e\n", i, krev(i));
    //
    fprintf(fs, ":: ropFor\n");
    for (int i = 0; i < int(ropFor.extent(0)); ++i)
      fprintf(fs, "     i %3d, kfor %e\n", i, ropFor(i));
    fprintf(fs, ":: ropRev\n");
    for (int i = 0; i < int(ropRev.extent(0)); ++i)
      fprintf(fs, "     i %3d, krev %e\n", i, ropRev(i));
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

template <typename MemberType, typename WorkViewType,
          typename RealType1DViewType, typename KineticModelConstDataType>
KOKKOS_FORCEINLINE_FUNCTION static void team_invoke(
    const MemberType& member,
    /// input
    const RealType1DViewType& kfor, const RealType1DViewType& krev,
    const RealType1DViewType& X,  /// (kmcd.nSpec)
    /// output
    const RealType1DViewType& omega,  /// (kmcd.nSpec)
    /// workspace
    const WorkViewType& work,
    /// const input from kinetic model
    const KineticModelConstDataType& kmcd) {
  /// FIXME: not exactly clear on what's happening here with work/iter
  auto w = (Real*)work.data();
  // auto mkfor = RealType1DViewType(w, kmcd.nReac);
  // w += kmcd.nReac;
  // auto mkrev = RealType1DViewType(w, kmcd.nReac);
  // w += kmcd.nReac;
  auto ropFor = RealType1DViewType(w, kmcd.nReac);
  w += kmcd.nReac;
  auto ropRev = RealType1DViewType(w, kmcd.nReac);
  w += kmcd.nReac;
  auto iter = Kokkos::View<ordinal_type*, Kokkos::LayoutRight,
                           typename WorkViewType::memory_space>(
      (ordinal_type*)w, kmcd.nReac * 2);
  w += kmcd.nReac * 2;

  team_invoke_detail(member, X, omega, kfor, krev, ropFor, ropRev, iter, kmcd);
}

template <typename PolicyType,
          // typename RealType1DViewType,
          typename RealType2DViewType, typename KineticModelConstType>
//
void SourceTermToyProblem_TemplateRun(  /// input
    const std::string& profile_name,
    /// team size setting
    const PolicyType& policy, const RealType2DViewType& kfor,
    const RealType2DViewType& krev, const RealType2DViewType& state,
    const RealType2DViewType& SourceTermToyProblem,
    const KineticModelConstType& kmcd) {
  Kokkos::Profiling::pushRegion(profile_name);

  using policy_type = PolicyType;

  const ordinal_type level = 1;
  const ordinal_type per_team_extent = getWorkSpaceSize(kmcd);

  Kokkos::parallel_for(
      profile_name, policy,
      KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
        const ordinal_type i = member.league_rank();
        const real_type_1d_view state_at_i =
            Kokkos::subview(state, i, Kokkos::ALL());
        const real_type_1d_view kfor_at_i =
            Kokkos::subview(kfor, i, Kokkos::ALL());
        const real_type_1d_view krev_at_i =
            Kokkos::subview(krev, i, Kokkos::ALL());
        const real_type_1d_view SourceTermToyProblem_at_i =
            Kokkos::subview(SourceTermToyProblem, i, Kokkos::ALL());

        Scratch<real_type_1d_view> work(member.team_scratch(level),
                                        per_team_extent);

        team_invoke(member, kfor_at_i, krev_at_i, state_at_i,
                    SourceTermToyProblem_at_i, work, kmcd);
      });
  Kokkos::Profiling::popRegion();
}

void SourceTermToyProblem::runDeviceBatch(  /// input
    typename UseThisTeamPolicy<exec_space>::type& policy,
    const real_type_2d_view& kfor, const real_type_2d_view& krev,
    const real_type_2d_view& state,
    /// output
    const real_type_2d_view& SourceTermToyProblem,
    /// const data from kinetic model
    const KineticModelConstDataDevice& kmcd) {
  SourceTermToyProblem_TemplateRun(  /// input
      "SourceTermToyProblem::runDeviceBatch",
      /// team size setting
      policy, kfor, krev, state, SourceTermToyProblem, kmcd);
}

}  // namespace from_tchem
