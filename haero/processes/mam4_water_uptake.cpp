#include "haero/model.hpp"
#include "haero/processes/mam4_water_uptake.hpp"

using haero::Real;

extern "C" {

// Fortran subroutine that computes water uptake.
extern void modal_aero_wateruptake_sub(int* num_cols,
                                       int* top_lev,
                                       int* num_levels,
                                       int* num_modes,
                                       bool* use_bisection,
                                       Real* rhcrystal,
                                       Real* rhdeliques,
                                       Real* dryrad,
                                       Real* naer,
                                       Real* hygro,
                                       Real* rh,
                                       Real* dryvol,
                                       Real* drymass,
                                       Real* specdens_1,
                                       Real* dgncur_a,
                                       Real* dgncur_awet,
                                       Real* qaerwat,
                                       Real* wetdens);

} // extern "C"

namespace {

/*
void compute_relative_humidity() {
  for (int i = 0; i < num_cols; ++i) {
    qsat_water(t(:ncol,k), pmid(:ncol,k), es(:ncol), qs(:ncol))
    for (int k = 0; k < num_levels; ++k) {
      if (qs(i) > h2ommr(i,k)) {
        rh(i,k) = h2ommr(i,k)/qs(i)
      } else {
        rh(i,k) = 0.98_wp
      }
      rh(i,k) = std::max(rh(i,k), 0.0_wp)
      rh(i,k) = std::min(rh(i,k), 0.98_wp)
      if(pergro_mods) {
        cldn_thresh = 0.9998_wp
      } else {
        cldn_thresh = 1.0_wp !original code
      }
      if (cldn(i,k) .lt. cldn_thresh) { // BSINGH
        rh(i,k) = (rh(i,k) - cldn(i,k)) / (1.0_wp - cldn(i,k))  ! clear portion
      }
      rh(i,k) = std::max(rh(i,k), 0.0_wp)
    }
  }
}
*/

}

namespace haero {

Mam4WaterUptake::Mam4WaterUptake(bool use_bisection):
  DiagnosticProcess(haero::WaterUptakeProcess, "MAM4 water uptake (Fortran)"),
  use_bisection_(use_bisection), rh(nullptr), maer(nullptr), hygro(nullptr),
  naer(nullptr), dryvol(nullptr), drymass(nullptr), dryrad(nullptr),
  rhcrystal(nullptr), rhdeliques(nullptr), specdens_1(nullptr) {
}

Mam4WaterUptake::~Mam4WaterUptake() {
  if (rh != nullptr) {
    delete [] rh;
    delete [] maer;
    delete [] hygro;
    delete [] naer;
    delete [] dryvol;
    delete [] drymass;
    delete [] dryrad;
    delete [] rhcrystal;
    delete [] rhdeliques;
    delete [] specdens_1;
  }
}

void Mam4WaterUptake::init(const Model& model) {
  int num_cols = model.num_columns();
  int num_levels = model.num_levels();
  int num_modes = model.modes().size();

  // Allocate memory for working arrays if we haven't already. All of these
  // working arrays have identical names to their counterparts in MAM4.
  if (rh == nullptr) {
    rh = new Real[num_cols * num_levels];
    maer = new Real[num_cols * num_levels * num_modes];
    hygro = new Real[num_cols * num_levels * num_modes];
    naer = new Real[num_cols * num_levels * num_modes];
    dryvol = new Real[num_cols * num_levels * num_modes];
    drymass = new Real[num_cols * num_levels * num_modes];
    dryrad = new Real[num_cols * num_levels * num_modes];
    rhcrystal = new Real[num_modes];
    rhdeliques = new Real[num_modes];
    specdens_1 = new Real[num_modes];
  }
}

void Mam4WaterUptake::update(const Model& model, Real t,
                             const Prognostics& prognostics,
                             Diagnostics& diagnostics) const {
  int num_cols = model.num_columns();
  int num_levels = model.num_levels();
  int num_modes = model.modes().size();

  // This implementation follows the logic in modal_aero_wateruptake_dr
  // (see modal_aero_wateruptake.F90 in our box model repo), which handles
  // the coupling between MAM and CAM.
  //
  // We follow the path in which list_idx = 0, where data are
  // extracted from the physics buffer. Here we extract our own data from
  // our model and state.

  // Extract diagnostic variables from the state.
  auto qaerwat = diagnostics.modal_var("aero_water");
  auto dgncur_a = diagnostics.modal_var("mean_diameter");
  auto dgncur_awet = diagnostics.modal_var("mean_wet_diameter");
  auto wetdens = diagnostics.modal_var("wet_density");
  auto aerosol_water = diagnostics.var("total_aero_water");

  // Compute the relative humidity.

  // Call the Fortran water uptake subroutine.
  int top_lev = 1; // top of atmosphere is the first entry in the column
  bool bisect = use_bisection_;
  modal_aero_wateruptake_sub(&num_cols, &top_lev, &num_levels, &num_modes,
                             &bisect, &rhcrystal[0], &rhdeliques[0], dryrad,
                             naer, hygro, rh, dryvol, drymass, specdens_1,
                             &(*dgncur_a.data())[0],
                             &(*dgncur_awet.data())[0],
                             &(*qaerwat.data())[0],
                             &(*wetdens.data())[0]);
}

}

