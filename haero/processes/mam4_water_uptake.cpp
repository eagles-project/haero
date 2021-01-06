#include "haero/model.hpp"
#include "haero/processes/mam4_water_uptake.hpp"

using haero::Real;

extern "C" {

// Fortran subroutines.
extern void* wateruptake_init(void* model);
extern void* wateruptake_init_with_bisection(void* model);
extern void wateruptake_update(void* process,
                               void* model,
                               Real t,
                               void* prognostics,
                               void* diagnostics);
extern void wateruptake_destroy(void* process);

} // extern "C"

namespace haero {

Mam4WaterUptake::Mam4WaterUptake(bool use_bisection):
  FDiagnosticProcess(haero::WaterUptakeProcess, "MAM4 water uptake (Fortran)",
                     use_bisection ? wateruptake_init_with_bisection
                                   : wateruptake_init,
                     wateruptake_update, wateruptake_destroy) {

  // List the diagnostic variables needed/updated by this process.
  auto vars = std::vector<std::string>({"total_aero_water"});
  auto modal_vars = std::vector<std::string>({"aero_water", "mean_diameter",
    "mean_wet_diameter", "wet_density"});
  set_diag_vars(vars, modal_vars);
}

} // end haero namespace

