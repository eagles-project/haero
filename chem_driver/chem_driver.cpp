#include "chem_driver.hpp"

#include <cstdarg>
#include <cstdio>

#include "haero/physical_constants.hpp"

namespace haero {
namespace chem_driver {

namespace {
// default reaction rate coefficients
// Note: these are the default values used by Musica/MusicBox
const std::map<std::string, Real> default_rxn_coeffs = {
    {"A", 1.0},      {"Ea", 0.0},     {"B", 0.0},      {"D", 300.0},
    {"E", 0.0},      {"k0_A", 1.0},   {"k0_B", 0.0},   {"k0_C", 0.0},
    {"kinf_A", 1.0}, {"kinf_B", 0.0}, {"kinf_C", 0.0}, {"Fc", 0.6},
    {"N", 1.0}};
}  // namespace

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
}  // end anonymous namespace

/// ChemicalSpecies constructor
ChemicalSpecies::ChemicalSpecies(const std::string& mname, Real minitial_value,
                                 const std::string& munits)
    : name(mname), initial_value(minitial_value), units(munits) {}

/// EnvironmentalConditions constructor
EnvironmentalConditions::EnvironmentalConditions(Real T0,
                                                 const std::string& T_units,
                                                 Real P0,
                                                 const std::string& P_units)
    : initial_temp(T0),
      units_temp(T_units),
      initial_pressure(P0),
      units_pressure(P_units) {}

/// Reaction constructor: enumerates the reaction type and sets default rate
/// coefficients if values are not provided in the yaml input
// FIXME: currently only have arrhenius and troe, but we'll cross that bridge
// when we come to it
Reaction::Reaction(const std::string& mtype_str,
                   const std::map<std::string, Real>& mreactants,
                   const std::map<std::string, Real>& mproducts,
                   const std::map<std::string, Real>& mrate_coefficients)
    : type_str(mtype_str), reactants(mreactants), products(mproducts) {
  // convert reaction type string to lowercase
  transform(type_str.begin(), type_str.end(), type_str.begin(), ::tolower);
  // use the type string to set the enumerated type for the reaction and assign
  // the rate coefficients, either based on input or to default values
  if (type_str.compare("arrhenius") == 0) {
    // assign the corresponding enum
    type = arrhenius;
    // get the rate coefficient from the SimulationInput, otherwise use the
    // default value
    get_or_default(mrate_coefficients, "A");
    get_or_default(mrate_coefficients, "Ea");
    get_or_default(mrate_coefficients, "B");
    get_or_default(mrate_coefficients, "D");
    get_or_default(mrate_coefficients, "E");
  } else if (type_str.compare("troe") == 0) {
    type = troe;
    get_or_default(mrate_coefficients, "k0_A");
    get_or_default(mrate_coefficients, "k0_B");
    get_or_default(mrate_coefficients, "k0_C");
    get_or_default(mrate_coefficients, "kinf_A");
    get_or_default(mrate_coefficients, "kinf_B");
    get_or_default(mrate_coefficients, "kinf_C");
    get_or_default(mrate_coefficients, "Fc");
    get_or_default(mrate_coefficients, "N");
  } else {
    std::cout << "ERROR: reaction type currently unsupported."
              << "\n";
  }
}

/// method for either getting the coefficient from SimulationInput
/// or using defaults
void Reaction::get_or_default(
    const std::map<std::string, Real>& mrate_coefficients,
    const std::string& name) {
  // determine whether a rate coefficient was provided by the input yaml
  // if not, assign that coefficient the default value
  // Note: find() returns vec.end() if it is not found
  rate_coefficients[name] =
      (mrate_coefficients.find(name) != mrate_coefficients.end())
          ? mrate_coefficients.at(name)
          : default_rxn_coeffs.at(name);
}

/// ChemSolver constructor: initializes all the required views on device and
/// sets some kokkos-related parameters
ChemSolver::ChemSolver(SimulationInput& sim_inp)
    : reactions(sim_inp.reactions) {
  // parse the tchem section of the input yaml
  parse_tchem_inputs(sim_inp);

  // set the kokkos parallel policy
  // FIXME: determine whether this lines up with the overall haero goals, as
  // what's here was copied over from the TChem implementation
  policy = policy_type(TChem::exec_space(), nbatch, Kokkos::AUTO());

  // set the temp and pressure that we got from the simulation input
  temperature = sim_inp.env_conditions.initial_temp;
  pressure = sim_inp.env_conditions.initial_pressure;
  units_temp = sim_inp.env_conditions.units_temp;
  units_pressure = sim_inp.env_conditions.units_pressure;

  // set the initial state, based on input
  int nSpecies = sim_inp.species.size();
  state = real_type_2d_view("StateVector", nbatch, nSpecies);
  auto state_host = Kokkos::create_mirror_view(state);
  for (int i = 0; i < nbatch; ++i) {
    for (int j = 0; j < nSpecies; ++j) {
      state_host(i, j) = sim_inp.species[j].initial_value;
    }
  }
  Kokkos::deep_copy(state, state_host);

  // optionally print some configuration info
  // (note: this doesn't appear to do anything when compiled for serial)
  if (verbose) {
    TChem::exec_space::print_configuration(std::cout, verbose);
    TChem::host_exec_space::print_configuration(std::cout, verbose);
  }

  // construct the KMD object
  // create kmd via yaml input?
  kmd = TChem::KineticModelData(cfiles.chemFile, cfiles.thermFile);
  std::ofstream output(cfiles.thermFile);
  // create a const version, as required by TChem/kokkos
  kmcd = kmd.createConstData<TChem::exec_space>();

  // initialize omega (tendency output)
  omega = real_type_2d_view("NetProductionRate", nbatch, kmcd.nSpec);

  // FIXME: this bit is a little over my head--should probably consider the
  // implications of what's happening here
  const ordinal_type level = 1;
  const ordinal_type per_team_extent = getWorkSpaceSize(kmcd);
  const ordinal_type per_team_scratch =
      TChem::Scratch<real_type_1d_view>::shmem_size(per_team_extent);
  policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));
}

/// print summary information about ChemSolver, as it relates to TChem
void ChemSolver::print_summary(const ChemFiles& cfiles) {
  printf("---------------------------------------------------\n");
  printf(
      "Testing Arguments: \n batch size %d\n chemfile %s\n thermfile %s\n "
      "outputfile %s\n verbose %s\n",
      nbatch, cfiles.chemFile.c_str(), cfiles.thermFile.c_str(),
      // inputFile.c_str(),
      cfiles.outputFile.c_str(), verbose ? "true" : "false");
  printf("---------------------------------------------------\n");
  printf("Time reaction rates %e [sec] %e [sec/sample]\n", t_device_batch,
         t_device_batch / Real(nbatch));
}  // end ChemSolver::print_summary

/// get
real_type_2d_view ChemSolver::get_tendencies() {
  // reset timer
  timer.reset();

  // calculate the reaction rates (based on current temp/pressure and the given
  // rate coefficients)
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

    /// print the first (of nbatch) values
    {
      auto omega_host_at_0 = Kokkos::subview(omega_host, 0, Kokkos::ALL());
      TChem::Test::writeReactionRates(cfiles.outputFile, kmcd.nSpec,
                                      omega_host_at_0);
    }
  }
  return omega;
}  // end ChemSolver::get_tendencies

/// get the TChem-specific inputs from the input yaml
void ChemSolver::parse_tchem_inputs(SimulationInput& sim_inp) {
  auto root = YAML::LoadFile(sim_inp.input_file);
  if (root["tchem"] and root["tchem"].IsMap()) {
    auto node = root["tchem"];
    if (not node["nbatch"]) {
      throw YamlException(
          "problem specific entry does not specify number "
          "of batches (nbatch).");
    } else if (not node["verbose"]) {
      throw YamlException(
          "problem specific entry has no verbose boolean (verbose).");
    } else {
      nbatch = node["nbatch"].as<int>();
      verbose = node["verbose"].as<bool>();
    }
  } else {
    printf(
        "No tchem section was found--using default values: verbose = false"
        ", nbatch = 1.\n");
    // FIXME: depending on how we ultimately parallelize, the nbatch default
    // may have to be set more cleverly
    nbatch = 1;
    verbose = false;
  }
}

/// set the reaction rates, given the coefficient from input and the current
/// temperature and pressure
// FIXME: currently the ordering here corresponds to the hard-coded input files
// that are generated in toy_problem.cpp (and consequently the input yaml)
// --will probably have to be careful when we update the TChem input file
// handling
// FIXME: may want to consider doing this in parallel, closer to the tendency
// calculations, depending on how we intend to parallelize over batches in haero
// FIXME: make this work for Troe reactions, when the time comes
void ChemSolver::set_reaction_rates() {
  int nRxn = reactions.size();
  std::vector<Real> rxn_rate(nRxn);

  // make the rate calculation for each reaction in the reactions vector
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
      // FIXME: WIP
      case troe: {
        rxn_rate[i] = 2.5;
        break;
      }
    }
  }

  // NOTE: so far, we aren't considering reversible reactions,
  // so we set those rates to zero--this could change
  kfor = real_type_2d_view("ForwardRate", nbatch, reactions.size());
  krev = real_type_2d_view("ReverseRate", nbatch, reactions.size());

  /// create a mirror view to store input
  auto kfor_host = Kokkos::create_mirror_view(kfor);
  auto krev_host = Kokkos::create_mirror_view(krev);

  // assign the calculated rates to the corresponding mirror view
  for (int i = 0; i < nbatch; ++i) {
    for (int j = 0; j < nRxn; ++j) {
      kfor_host(i, j) = rxn_rate[j];
      krev_host(i, j) = 0;
    }
  }

  // deep copy to device
  Kokkos::deep_copy(kfor, kfor_host);
  Kokkos::deep_copy(krev, krev_host);
}

/// ChemSolver destructor that removes the temporary TChem input files from disk
ChemSolver::~ChemSolver() {
  // delete the temporary chem files
  // NOTE: for now, removing the output file, too
  remove(cfiles.chemFile.data());
  remove(cfiles.thermFile.data());
  remove(cfiles.outputFile.data());
}

}  // end namespace chem_driver
}  // end namespace haero

/*
NOTE: everything in this namespace is copied directly from TChem
*/
namespace from_tchem {

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
          const Real val = kmcd.reacNuki(i, j) * rop_at_i;
          Kokkos::atomic_fetch_add(&omega(kspec), val);
        }
        const ordinal_type joff = kmcd.reacSidx.extent(1) / 2;
        for (ordinal_type j = 0; j < kmcd.reacNprod(i); ++j) {
          const ordinal_type kspec = kmcd.reacSidx(i, j + joff);
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
  // mjs: keep an eye on this to make sure it's working properly
  auto w = (Real*)work.data();
  auto mkfor = RealType1DViewType(w, kmcd.nReac);
  w += kmcd.nReac;
  auto mkrev = RealType1DViewType(w, kmcd.nReac);
  w += kmcd.nReac;
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

template <typename PolicyType, typename RealType2DViewType,
          typename KineticModelConstType>
void SourceTermToyProblem_TemplateRun(
    /// input
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

void SourceTermToyProblem::runDeviceBatch(
    /// input
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

}  // end namespace from_tchem
