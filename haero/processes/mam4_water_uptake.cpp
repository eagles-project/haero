#include "haero/model.hpp"
#include "haero/processes/mam4_water_uptake.hpp"

using haero::Real;

extern "C" {

// Fortran subroutines.
extern void wateruptake_init(FModel* model,
                             bool use_bisection);
extern void wateruptake_update(FModel* model,
                               Real t,
                               FPrognostics* prognostics,
                               FDiagnostics* diagnostics);

} // extern "C"

namespace haero {

Mam4WaterUptake::Mam4WaterUptake(bool use_bisection):
  DiagnosticProcess(haero::WaterUptakeProcess, "MAM4 water uptake (Fortran)"),
  use_bisection_(use_bisection) {

  // List the diagnostic variables needed/updated by this process.
  auto vars = std::vector<std::string>({"total_aero_water"});
  auto modal_vars = std::vector<std::string>({"aero_water", "mean_diameter",
    "mean_wet_diameter", "wet_density"});
  set_diag_vars(vars, modal_vars);
}

Mam4WaterUptake::~Mam4WaterUptake() {
}

void Mam4WaterUptake::init(const Model& model) {
  wateruptake_init(model.to_fortran(), use_bisection_);
}

void Mam4WaterUptake::update(const Model& model, Real t,
                             const Prognostics& prognostics,
                             Diagnostics& diagnostics) const {
  wateruptake_update(model.to_fortran(), t, prognostics.to_fortran(),
                     diagnostics.to_fortran());
}

} // end haero namespace

