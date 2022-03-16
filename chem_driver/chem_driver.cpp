#include "chem_driver.hpp"

namespace haero {
namespace chem_driver {

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
};  // end anonymous namespace

// the fields written to file are: iteration, t, dt, Density, Pressure,
// Temperature, "concentration(s)" (in the relevant, specified units)
void static writeState(const ordinal_type iter, const Real_1d_view_host _t,
                       const Real_1d_view_host _dt,
                       const Real_2d_view_host _state_at_i, FILE* fout_) {
  for (size_t sp = 0; sp < _state_at_i.extent(0); sp++) {
    fprintf(fout_, "%d \t %15.10e \t  %15.10e \t ", iter, _t(sp), _dt(sp));
    //
    for (ordinal_type k = 0, kend = _state_at_i.extent(1); k < kend; ++k)
      fprintf(fout_, "%15.10e \t", _state_at_i(sp, k));

    fprintf(fout_, "\n");
  }
};

void static printState(const time_advance_type _tadv, const double _t,
                       const Real_1d_view_host _state_at_i) {
  /// iter, t, dt, rho, pres, temp, Ys ....
  printf("%e %e %e %e %e", _t, _t - _tadv._tbeg, _state_at_i(0), _state_at_i(1),
         _state_at_i(2));
  for (ordinal_type k = 3, kend = _state_at_i.extent(0); k < kend; ++k)
    printf(" %e", _state_at_i(k));
  printf("\n");
};

/// ChemSolver constructor: initializes all the required views on device and
/// sets some kokkos-related parameters
ChemSolver::ChemSolver(std::string input_file) {
  // parse the tchem section of the input yaml
  parse_tchem_inputs_(input_file);
  // read the parameters from file
  solver_params_.set_params(input_file, verbose_);

  TChem::exec_space::print_configuration(std::cout, verbose_);
  TChem::host_exec_space::print_configuration(std::cout, verbose_);
  using device_type = typename Tines::UseThisDevice<exec_space>::type;

  /// construct kmd and use the view for testing
  kmd_ = TChem::KineticModelData(input_file);
  kmcd_ = TChem::createNCAR_KineticModelConstData<device_type>(kmd_);

  const ordinal_type stateVecDim = TChem::Impl::getStateVectorSize(kmcd_.nSpec);

  if (verbose_) {
    printf("Number of Species %d \n", kmcd_.nSpec);
    printf("Number of Reactions %d \n", kmcd_.nReac);
  }
  const auto speciesNamesHost = Kokkos::create_mirror_view(kmcd_.speciesNames);
  Kokkos::deep_copy(speciesNamesHost, kmcd_.speciesNames);

  fout_ = fopen(solver_params_.outputfile.c_str(), "w");
  // read scenario condition from yaml file
  TChem::AtmChemistry::setScenarioConditions(input_file, speciesNamesHost,
                                             kmcd_.nSpec, state_host_, nbatch_);
  state_ = Real_2d_view("StateVector Devices", nbatch_, stateVecDim);

  // Kokkos::Impl::Timer timer;

  // timer.reset();
  Kokkos::deep_copy(state_, state_host_);
  // const double t_deepcopy = timer.seconds();

  const auto exec_space_instance = TChem::exec_space();

  /// team policy
  policy_ = policy_type(exec_space_instance, nbatch_, Kokkos::AUTO());

  if (team_size_ > 0 && vector_size_ > 0) {
    policy_ = policy_type(exec_space_instance, nbatch_, team_size_, vector_size_);
  }

  const ordinal_type level = 1;
  ordinal_type per_team_extent(0);

  per_team_extent = TChem::AtmosphericChemistry::getWorkSpaceSize(kmcd_);

  const ordinal_type per_team_scratch =
      TChem::Scratch<Real_1d_view>::shmem_size(per_team_extent);
  policy_.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));
}

// time integrator that gets tbeg and tend from yaml input (stored in
// solver_params_)
void ChemSolver::time_integrate() {
  Real_1d_view t("time", nbatch_);
  Kokkos::deep_copy(t, solver_params_.tbeg);
  Real_1d_view dt("delta time", nbatch_);
  Kokkos::deep_copy(dt, solver_params_.dtmax);

  Real_1d_view_host t_host;
  Real_1d_view_host dt_host;

  t_host = Real_1d_view_host("time host", nbatch_);
  dt_host = Real_1d_view_host("dt host", nbatch_);

  ordinal_type number_of_equations(0);

  using device_type = typename Tines::UseThisDevice<exec_space>::type;
  using problem_type =
      TChem::Impl::AtmosphericChemistry_Problem<double, device_type>;
  number_of_equations = problem_type::getNumberOfTimeODEs(kmcd_);

  Real_2d_view tol_time("tol time", number_of_equations, 2);
  Real_1d_view tol_newton("tol newton", 2);

  Real_2d_view fac("fac", nbatch_, number_of_equations);

  /// copy tolerances to device
  {
    auto tol_time_host = Kokkos::create_mirror_view(tol_time);
    auto tol_newton_host = Kokkos::create_mirror_view(tol_newton);

    for (ordinal_type i = 0, iend = tol_time.extent(0); i < iend; ++i) {
      tol_time_host(i, 0) = solver_params_.atol_time;
      tol_time_host(i, 1) = solver_params_.tol_time;
    }
    tol_newton_host(0) = solver_params_.atol_newton;
    tol_newton_host(1) = solver_params_.rtol_newton;

    Kokkos::deep_copy(tol_time, tol_time_host);
    Kokkos::deep_copy(tol_newton, tol_newton_host);
  }

  time_advance_type tadv_default;
  tadv_default._tbeg = solver_params_.tbeg;
  tadv_default._tend = solver_params_.tend;
  tadv_default._dt = solver_params_.dtmin;
  tadv_default._dtmin = solver_params_.dtmin;
  tadv_default._dtmax = solver_params_.dtmax;
  tadv_default._max_num_newton_iterations = solver_params_.max_newton_iterations;
  tadv_default._num_time_iterations_per_interval =
      solver_params_.num_time_iterations_per_interval;
  tadv_default._jacobian_interval = solver_params_.jacobian_interval;

  time_advance_type_1d_view tadv("tadv", nbatch_);
  Kokkos::deep_copy(tadv, tadv_default);

  /// host views to print QOI
  const auto tadv_at_i = Kokkos::subview(tadv, 0);
  const auto t_at_i = Kokkos::subview(t, 0);
  const auto state_at_i = Kokkos::subview(state_, 0, Kokkos::ALL());

  auto tadv_at_i_host = Kokkos::create_mirror_view(tadv_at_i);
  auto t_at_i_host = Kokkos::create_mirror_view(t_at_i);
  auto state_at_i_host = Kokkos::create_mirror_view(state_at_i);

  ordinal_type iter = 0;
  /// print of store QOI for the first sample
  if (print_qoi_) {
    Kokkos::deep_copy(tadv_at_i_host, tadv_at_i);
    Kokkos::deep_copy(t_at_i_host, t_at_i);
    Kokkos::deep_copy(state_at_i_host, state_at_i);
    printState(tadv_at_i_host(), t_at_i_host(), state_at_i_host);
  }

  Kokkos::deep_copy(dt_host, dt);
  Kokkos::deep_copy(t_host, t);

  fprintf(fout_, "%s \t %s \t %s \t ", "iter", "t", "dt");
  fprintf(fout_, "%s \t %s \t %s \t", "Density[kg/m3]", "Pressure[Pascal]",
          "Temperature[K]");

  const auto speciesNamesHost = Kokkos::create_mirror_view(kmcd_.speciesNames);
  Kokkos::deep_copy(speciesNamesHost, kmcd_.speciesNames);

  for (ordinal_type k = 0; k < kmcd_.nSpec; k++) {
    fprintf(fout_, "%s \t", &speciesNamesHost(k, 0));
  }
  fprintf(fout_, "\n");

  writeState(-1, t_host, dt_host, state_host_, fout_);

  double tsum(0);
  for (; iter < solver_params_.max_time_iterations &&
         tsum <= solver_params_.tend * 0.9999;
       ++iter) {
    TChem::AtmosphericChemistry::runDeviceBatch(
        policy_, tol_newton, tol_time, fac, tadv, state_, t, dt, state_, kmcd_);

    if (print_qoi_) {
      /// could use cuda streams
      Kokkos::deep_copy(tadv_at_i_host, tadv_at_i);
      Kokkos::deep_copy(t_at_i_host, t_at_i);
      Kokkos::deep_copy(state_at_i_host, state_at_i);
      printState(tadv_at_i_host(), t_at_i_host(), state_at_i_host);
    }

    Kokkos::deep_copy(dt_host, dt);
    Kokkos::deep_copy(t_host, t);
    Kokkos::deep_copy(state_host_, state_);

    writeState(iter, t_host, dt_host, state_host_, fout_);

    /// carry over time and dt computed in this step
    double tsum(0);
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<TChem::exec_space>(0, nbatch_),
        KOKKOS_LAMBDA(const ordinal_type& i, double& update) {
          tadv(i)._tbeg = t(i);
          tadv(i)._dt = dt(i);
          update += t(i);
        },
        tsum);
    Kokkos::fence();
    tsum /= nbatch_;
  }                 // end for
  Kokkos::fence();  /// timing purpose

  if (print_qoi_) {
    Kokkos::parallel_for(
        Kokkos::RangePolicy<TChem::exec_space>(0, nbatch_),
        KOKKOS_LAMBDA(const ordinal_type& i) {
          printf("Devices:: Solution sample No %d\n", i);
          auto state_at_i = Kokkos::subview(state_, i, Kokkos::ALL());
          for (ordinal_type k = 0, kend = state_at_i.extent(0); k < kend; ++k)
            printf(" %e", state_at_i(k));
          printf("\n");
        });
  }
}

// time integrator that takes tbeg and tend as arguments--useful for debugging
void ChemSolver::time_integrate(const double& tbeg, const double& tend) {
  Real_1d_view t("time", nbatch_);
  Kokkos::deep_copy(t, tbeg);
  Real_1d_view dt("delta time", nbatch_);
  Kokkos::deep_copy(dt, solver_params_.dtmax);

  Real_1d_view_host t_host;
  Real_1d_view_host dt_host;

  t_host = Real_1d_view_host("time host", nbatch_);
  dt_host = Real_1d_view_host("dt host", nbatch_);

  ordinal_type number_of_equations(0);

  using device_type = typename Tines::UseThisDevice<exec_space>::type;
  using problem_type =
      TChem::Impl::AtmosphericChemistry_Problem<double, device_type>;
  number_of_equations = problem_type::getNumberOfTimeODEs(kmcd_);

  Real_2d_view tol_time("tol time", number_of_equations, 2);
  Real_1d_view tol_newton("tol newton", 2);

  Real_2d_view fac("fac", nbatch_, number_of_equations);

  /// copy tolerances to device
  {
    auto tol_time_host = Kokkos::create_mirror_view(tol_time);
    auto tol_newton_host = Kokkos::create_mirror_view(tol_newton);

    for (ordinal_type i = 0, iend = tol_time.extent(0); i < iend; ++i) {
      tol_time_host(i, 0) = solver_params_.atol_time;
      tol_time_host(i, 1) = solver_params_.tol_time;
    }
    tol_newton_host(0) = solver_params_.atol_newton;
    tol_newton_host(1) = solver_params_.rtol_newton;

    Kokkos::deep_copy(tol_time, tol_time_host);
    Kokkos::deep_copy(tol_newton, tol_newton_host);
  }

  time_advance_type tadv_default;
  tadv_default._tbeg = tbeg;
  tadv_default._tend = tend;
  tadv_default._dt = solver_params_.dtmin;
  tadv_default._dtmin = solver_params_.dtmin;
  tadv_default._dtmax = solver_params_.dtmax;
  tadv_default._max_num_newton_iterations = solver_params_.max_newton_iterations;
  tadv_default._num_time_iterations_per_interval =
      solver_params_.num_time_iterations_per_interval;
  tadv_default._jacobian_interval = solver_params_.jacobian_interval;

  time_advance_type_1d_view tadv("tadv", nbatch_);
  Kokkos::deep_copy(tadv, tadv_default);

  /// host views to print QOI
  const auto tadv_at_i = Kokkos::subview(tadv, 0);
  const auto t_at_i = Kokkos::subview(t, 0);
  const auto state_at_i = Kokkos::subview(state_, 0, Kokkos::ALL());

  auto tadv_at_i_host = Kokkos::create_mirror_view(tadv_at_i);
  auto t_at_i_host = Kokkos::create_mirror_view(t_at_i);
  auto state_at_i_host = Kokkos::create_mirror_view(state_at_i);

  /// print of store QOI for the first sample
  if (print_qoi_) {
    Kokkos::deep_copy(tadv_at_i_host, tadv_at_i);
    Kokkos::deep_copy(t_at_i_host, t_at_i);
    Kokkos::deep_copy(state_at_i_host, state_at_i);
    printState(tadv_at_i_host(), t_at_i_host(), state_at_i_host);
  }

  Kokkos::deep_copy(dt_host, dt);
  Kokkos::deep_copy(t_host, t);

  fprintf(fout_, "%s \t %s \t %s \t ", "iter", "t", "dt");
  fprintf(fout_, "%s \t %s \t %s \t", "Density[kg/m3]", "Pressure[Pascal]",
          "Temperature[K]");

  const auto speciesNamesHost = Kokkos::create_mirror_view(kmcd_.speciesNames);
  Kokkos::deep_copy(speciesNamesHost, kmcd_.speciesNames);

  for (ordinal_type k = 0; k < kmcd_.nSpec; k++) {
    fprintf(fout_, "%s \t", &speciesNamesHost(k, 0));
  }
  fprintf(fout_, "\n");

  writeState(-1, t_host, dt_host, state_host_, fout_);

  double tsum(0);
  ordinal_type iter = 0;
  // note that this stops according to the tend passed to the function
  for (; iter < solver_params_.max_time_iterations && tsum <= tend * 0.9999;
       ++iter) {
    TChem::AtmosphericChemistry::runDeviceBatch(
        policy_, tol_newton, tol_time, fac, tadv, state_, t, dt, state_, kmcd_);

    if (print_qoi_) {
      /// could use cuda streams
      Kokkos::deep_copy(tadv_at_i_host, tadv_at_i);
      Kokkos::deep_copy(t_at_i_host, t_at_i);
      Kokkos::deep_copy(state_at_i_host, state_at_i);
      printState(tadv_at_i_host(), t_at_i_host(), state_at_i_host);
    }

    Kokkos::deep_copy(dt_host, dt);
    Kokkos::deep_copy(t_host, t);
    Kokkos::deep_copy(state_host_, state_);

    writeState(iter, t_host, dt_host, state_host_, fout_);

    /// carry over time and dt computed in this step
    double tsum(0);
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<TChem::exec_space>(0, nbatch_),
        KOKKOS_LAMBDA(const ordinal_type& i, double& update) {
          tadv(i)._tbeg = t(i);
          tadv(i)._dt = dt(i);
          update += t(i);
        },
        tsum);
    Kokkos::fence();
    tsum /= nbatch_;
  }                 // end for
  Kokkos::fence();  /// timing purpose

  if (print_qoi_) {
    Kokkos::parallel_for(
        Kokkos::RangePolicy<TChem::exec_space>(0, nbatch_),
        KOKKOS_LAMBDA(const ordinal_type& i) {
          printf("Devices:: Solution sample No %d\n", i);
          auto state_at_i = Kokkos::subview(state_, i, Kokkos::ALL());
          for (ordinal_type k = 0, kend = state_at_i.extent(0); k < kend; ++k)
            printf(" %e", state_at_i(k));
          printf("\n");
        });
  }
}

/// get the TChem-specific inputs from the input yaml
void ChemSolver::parse_tchem_inputs_(const std::string& input_file) {
  auto root = YAML::LoadFile(input_file);
  if (root["tchem"] and root["tchem"].IsMap()) {
    auto node = root["tchem"];
    if (not node["nbatch"]) {
      throw YamlException(
          "problem specific entry does not specify number "
          "of batches (nbatch).");
    } else if (not node["verbose"]) {
      throw YamlException(
          "problem specific entry has no verbose boolean (verbose).");
    } else if (not node["team_size"]) {
      throw YamlException(
          "problem specific entry has no team_size entry (team_size).");
    } else if (not node["vector_size"]) {
      throw YamlException(
          "problem specific entry has no vector_size entry (vector_size).");
    } else if (not node["print_qoi"]) {
      throw YamlException(
          "problem specific entry has no print_qoi boolean (verbose).");
    } else {
      nbatch_ = node["nbatch"].as<int>();
      verbose_ = node["verbose"].as<bool>();
      team_size_ = node["team_size"].as<int>();
      vector_size_ = node["vector_size"].as<int>();
      print_qoi_ = node["print_qoi"].as<bool>();
    }
  } else {
    printf(
        "No tchem section was found--using default values: verbose = false"
        ", nbatch_ = 1.\n");
    // FIXME: depending on how we ultimately parallelize, the nbatch_ default
    // may have to be set more cleverly
    nbatch_ = 1;
    verbose_ = false;
    team_size_ = -1;
    vector_size_ = -1;
  }
}

/// function to read the chemistry input yaml file and construct a
/// SolverParams struct from what is found there
void SolverParams::set_params(const std::string& filename,
                              const bool& verbose_) {
  // load the input from the yaml file
  auto root = YAML::LoadFile(filename);
  if (root["solver_parameters"] and root["solver_parameters"].IsMap()) {
    auto node = root["solver_parameters"];
    if (not node["dtmin"]) {
      throw YamlException(
          "solver_parameters contains no dtmin "
          "entry (dtmin).");
    }
    if (not node["dtmax"]) {
      throw YamlException(
          "solver_parameters contains no dtmax "
          "entry (dtmax).");
    }
    if (not node["tbeg"]) {
      throw YamlException(
          "solver_parameters contains no tbeg "
          "entry (tbeg).");
    }
    if (not node["tend"]) {
      throw YamlException(
          "solver_parameters contains no tend "
          "entry (tend).");
    }
    if (not node["num_time_iterations_per_interval"]) {
      throw YamlException(
          "solver_parameters contains no num_time_iterations_per_interval "
          "entry (num_time_iterations_per_interval).");
    }
    if (not node["max_time_iterations"]) {
      throw YamlException(
          "solver_parameters contains no max_time_iterations "
          "entry (max_time_iterations).");
    }
    if (not node["max_newton_iterations"]) {
      throw YamlException(
          "solver_parameters contains no max_newton_iterations "
          "entry (max_newton_iterations).");
    }
    if (not node["atol_newton"]) {
      throw YamlException(
          "solver_parameters contains no atol_newton "
          "entry (atol_newton).");
    }
    if (not node["rtol_newton"]) {
      throw YamlException(
          "solver_parameters contains no rtol_newton "
          "entry (rtol_newton).");
    }
    if (not node["atol_time"]) {
      throw YamlException(
          "solver_parameters contains no atol_time "
          "entry (atol_time).");
    }
    if (not node["tol_time"]) {
      throw YamlException(
          "solver_parameters contains no tol_time "
          "entry (tol_time).");
    }
    if (not node["jacobian_interval"]) {
      throw YamlException(
          "solver_parameters contains no jacobian_interval "
          "entry (jacobian_interval).");
    }
    if (not node["outputfile"]) {
      throw YamlException(
          "solver_parameters contains no outputfile "
          "entry (outputfile).");
    } else {
      dtmin = node["dtmin"].as<double>();
      dtmax = node["dtmax"].as<double>();
      tbeg = node["tbeg"].as<double>();
      tend = node["tend"].as<double>();
      num_time_iterations_per_interval =
          node["num_time_iterations_per_interval"].as<int>();
      max_time_iterations = node["max_time_iterations"].as<int>();
      max_newton_iterations = node["max_newton_iterations"].as<int>();
      atol_newton = node["atol_newton"].as<double>();
      rtol_newton = node["rtol_newton"].as<double>();
      atol_time = node["atol_time"].as<double>();
      tol_time = node["tol_time"].as<double>();
      jacobian_interval = node["jacobian_interval"].as<int>();
      outputfile = node["outputfile"].as<std::string>();
    }
  } else {
    dtmin = 1.0e-8;
    dtmax = 1.0e-1;
    tend = 0.0;
    tend = 1.0;
    num_time_iterations_per_interval = 10;
    max_time_iterations = 1.0e3;
    max_newton_iterations = 100;
    atol_newton = 1.0e-10;
    rtol_newton = 1.0e-6;
    atol_time = 1.0e-12;
    tol_time = 1.0e-4;
    jacobian_interval = 1;
    outputfile = "chem.dat";
    if (verbose_) {
      printf("No solver_parameters section was found--using defaults\n");
    }
  }
}

/// ChemSolver destructor
ChemSolver::~ChemSolver() { fclose(fout_); }

}  // end namespace chem_driver
}  // end namespace haero
