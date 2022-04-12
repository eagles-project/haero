// This is a box model used to validate Haero against MAM's box model. It
// behaves identically to that box model.

#include "box_input.hpp"
#include "box_output.hpp"

#include <haero/aerosol_process.hpp>
#include <haero/atmosphere.hpp>
#include <haero/diagnostics.hpp>
#include <haero/haero.hpp>
#include <haero/modal_aerosol_config.hpp>
#include <haero/prognostics.hpp>

#include <ekat/ekat_session.hpp>

#include <cstdio>
#include <cstdlib>
#include <limits>

using namespace std;
using namespace haero;

namespace {

// Print version banner.
void print_banner() {
  string version = haero::version();
  string revision = haero::revision();
  bool not_for_science = haero::has_uncommitted_changes();
  printf("haero v%s (git revision %s)\n", version.c_str(), revision.c_str());
  if (not_for_science) {
    fprintf(stderr, "WARNING: haero was built with uncommitted changeѕ!\n");
    fprintf(stderr, "WARNING: DO NOT USE FOR SCIENCE!\n");
  }
}

// Print driver usage information and exit.
void usage(const char* exe) {
  fprintf(stderr, "%s: usage:\n", exe);
  fprintf(stderr, "%s <input.nml>\n", exe);
  exit(1);
}

// Create the processes selected within the given input.
std::vector<AerosolProcess*> create_processes(const ModalAerosolConfig& config,
                                              const BoxInput& input) {
  std::vector<AerosolProcess*> processes;
  return processes;
}

// Create prognostic variables from input.
Prognostics create_prognostics(const ModalAerosolConfig& config,
                               const BoxInput& input) {
  Prognostics progs(config, 1);
  return progs;
}

// Create atmospheric variables from input.
Atmosphere create_atmosphere(const BoxInput& input) {
  static const Real pblh = 1100.0; // from the box model
  Atmosphere atm(1, pblh);
  return atm;
}

// Create diagnostic variables from input.
HostDiagnostics create_diagnostics(const ModalAerosolConfig& config,
                                   const BoxInput& input) {
  HostDiagnostics diags(config, 1);
  return diags;
}

// Advance the selected processes using the same sequential splitting approach
// used in MAM4.
void advance(const ModalAerosolConfig& config, Real t, Real dt,
             const std::vector<AerosolProcess*>& processes,
             Prognostics& prognostics,
             Atmosphere& atmosphere,
             HostDiagnostics& diagnostics,
             std::map<std::string, Tendencies>& process_tendencies) {
  // Call the processes in order, computing their tendencies and integrating
  // prognostic variables in place. We can do things this way because we happen
  // to know that the MAM4 processes do sequential updates and then back out the
  // tendencies using a finite difference.
  auto team_policy = haero::TeamPolicy(1u, Kokkos::AUTO);
  for (AerosolProcess *process: processes) {
    // Compute the tendencies for the prognostic variables.
    const std::string& name = process->name();
    Tendencies& tendencies = process_tendencies[name];
    Kokkos::parallel_for(
      team_policy, KOKKOS_LAMBDA(const TeamType& team) {
        process->run(team, t, dt, prognostics, atmosphere, diagnostics,
                     tendencies);
      });

    // Integrate the prognostic variables.
    for (int m = 0; m < config.num_aerosol_modes(); ++m) {
      prognostics.interstitial_num_mix_ratios(m, 0) +=
        dt * tendencies.interstitial_num_mix_ratios(m, 0);
      prognostics.cloud_num_mix_ratios(m, 0) +=
        dt * tendencies.cloud_num_mix_ratios(m, 0);
    }
    for (int p = 0; p < config.num_aerosol_populations; ++p) {
      prognostics.interstitial_aerosols(p, 0) +=
        dt * tendencies.interstitial_aerosols(p, 0);
      prognostics.cloud_aerosols(p, 0) +=
        dt * tendencies.cloud_aerosols(p, 0);
    }
    for (int g = 0; g < config.num_gases(); ++g) {
      prognostics.gases(g, 0) += dt * tendencies.gases(g, 0);
    }
  }
}

}  // anonymous namespace

int main(int argc, const char** argv) {
  ekat::initialize_ekat_session(argc, const_cast<char**>(argv), false);
  print_banner();

  if (argc < 2) usage(argv[0]);

  // Read the input file and extract input.
  std::string input_file(argv[1]);
  auto input = Box::read_namelist(input_file);

  // Initialize a set of processes with the given input.
  auto config = ModalAerosolConfig::create_mam4_config();
  std::vector<AerosolProcess*> processes = create_processes(config, input);

  // Initialize aerosol prognostics and atmospheric conditions with input.
  Prognostics prognostics = create_prognostics(config, input);
  Atmosphere atmosphere = create_atmosphere(input);
  HostDiagnostics diagnostics = create_diagnostics(config, input);

  // Create process-specific tendencies.
  std::map<std::string, Tendencies> process_tendencies;
  for (AerosolProcess* process: processes) {
    process_tendencies.emplace(process->name(), Tendencies(prognostics));
  }

  // Open the output file and define quantities.
  Box::BoxOutput output(config);

  Real t = 0.0;
  for (int n = 0; n < input.mam_nstep; ++n) {
    Real dt = static_cast<Real>(input.mam_dt);
    advance(config, t, dt, processes, prognostics, atmosphere, diagnostics,
            process_tendencies);
    t += dt;
    if ((n % input.mam_output_intvl) == 0) {
      output.append(prognostics, diagnostics, process_tendencies);
    }
  }

  // Write output.
  output.write("haero_output.nc");

  // Clean up.
  for (size_t i = 0; i < processes.size(); ++i) {
    delete processes[i];
  }
  ekat::finalize_ekat_session();
  return 0;
}
