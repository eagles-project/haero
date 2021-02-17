#include "chemUtil.hpp"

namespace chemUtil {

chemFiles::chemFiles(std::string chemDir){
  prefixPath = "data/" + chemDir;
  chemFile = prefixPath + "chem.inp";
  thermFile = prefixPath + "therm.dat";
  outputFile = prefixPath + "omega.dat";
}

chemSolver::chemSolver(std::string chemDir, bool detail, int inBatch,
                       bool iverbose, real_type itheta, real_type ilambda,
                       real_type initX, real_type initX2)
                       :
                       cfiles(chemDir),
                       verbose(iverbose),
                       nBatch(inBatch),
                       policy(TChem::exec_space(), nBatch, Kokkos::AUTO()),
                       theta("latitude", nBatch),
                       lambda("longitude", nBatch){

  // NOTE: these appear to do nothing, but keeping them around, just in case
  TChem::exec_space::print_configuration(std::cout, detail);
  TChem::host_exec_space::print_configuration(std::cout, detail);

  // construct kmd and use the view for testing
  std::ofstream output(cfiles.thermFile);
  TChem::KineticModelData kmd(cfiles.chemFile, cfiles.thermFile);
  // create a const version
  kmcd = kmd.createConstData<TChem::exec_space>();

  // input: state vectors: temperature, pressure and mass fraction
  state = real_type_2d_view("StateVector", nBatch, kmcd.nSpec);
  // output: omega, reaction rates
  omega = real_type_2d_view("NetProductionRate", nBatch, kmcd.nSpec);

  // assign the constructor arguments to the corresponding views
  theta(0) = itheta;
  lambda(0) = ilambda;

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
    TChem::Test::cloneView(state_host);
  }

  // reset timer and time the deep copy from host to device
  timer.reset();
  Kokkos::deep_copy(state, state_host);
  t_deepcopy = timer.seconds();

  // this bit is a little over my head
  const ordinal_type level = 1;
  const ordinal_type per_team_extent =
    TChem::SourceTermToyProblem::getWorkSpaceSize(kmcd);
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
         t_deepcopy / real_type(nBatch));
  printf("Time reaction rates %e [sec] %e [sec/sample]\n",
         t_device_batch,
         t_device_batch / real_type(nBatch));
} // end chemSolver::print_summary

real_type_2d_view chemSolver::get_results(){
  // reset timer
  timer.reset();
  // run the model
  TChem::SourceTermToyProblem::runDeviceBatch(policy, theta, lambda, state, omega, kmcd);
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

}
