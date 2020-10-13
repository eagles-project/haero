#include "haero/aero_model.hpp"
#include "haero/processes/mam4_water_uptake.hpp"

using haero::Real;

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

namespace {

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

namespace haero {

Mam4WaterUptake::Mam4WaterUptake(bool use_bisection):
  DiagnosticAeroProcess(WaterUptakeProcess, "MAM4 water uptake (Fortran)"),
  use_bisection_(use_bisection), rh(nullptr), maer(nullptr), hygro(nullptr),
  naer(nullptr), dryvol(nullptr), drymass(nullptr), dryrad(nullptr),
  rhcrystal(nullptr), rhdeliques(nullptr), ѕpecdens_1(nullptr) {
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

void Mam4WaterUptake::update(const AeroModel& model, Real t, AeroState& state) const {
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

  // Extract diagnostic variables from the state.
  auto qaerwat = state.diagnostic("aero_water");
  auto aerosol_water = state.diagnostic("total_aero_water");
  auto dgncur_awet = state.diagnostic("mean_wet_diameter");
  auto wetdens = state.diagnostic("wet_density");

  // Compute the relative humidity.
  compute_relative_humidity(

  // Call the Fortran water uptake subroutine.
  int top_lev = 1; // top of atmosphere is the first entry in the column
  modal_aero_wateruptake_sub(&num_cols, &top_lev, &num_levels, &num_modes,
                             &use_bisection_, &rhcrystal[0], &rhdeliques[0],
                             dryrad, naer, hygro, rh, dryvol, drymass,
                             specdens_1, dgncur_a, dgncur_awet, qaerwat, wetdens);

  // Store the updated quantities in the state.
  state.put("aero_water", qaerwat);
  state.put("total_aero_water", aerosol_water);
  state.put("mean_wet_diameter", dgncur_awet);
  state.put("wet_density", wetdens);
}

}

