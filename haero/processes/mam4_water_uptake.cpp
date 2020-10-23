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

namespace haero {

Mam4WaterUptake::Mam4WaterUptake(bool use_bisection):
  DiagnosticProcess(haero::WaterUptakeProcess, "MAM4 water uptake (Fortran)"),
  use_bisection_(use_bisection), rh(nullptr), maer(nullptr), hygro(nullptr),
  naer(nullptr), dryvol(nullptr), drymass(nullptr), dryrad(nullptr),
  rhcrystal(nullptr), rhdeliques(nullptr), specdens_1(nullptr) {

  // List the diagnostic variables needed/updated by this process.
  auto vars = std::vector<std::string>({"total_aero_water"});
  auto modal_vars = std::vector<std::string>({"aero_water", "mean_diameter",
    "mean_wet_diameter", "wet_density"});
  set_diag_vars(vars, modal_vars);
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
  int num_columns = model.num_columns();
  int num_levels = model.num_levels();
  const std::vector<Mode>& modes = model.modes();
  int num_modes = modes.size();

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

/*
  // Compute the volume-mean hygroscopicity and the dry radius, volume,
  // and mass.
  for (int m = 0; m < num_modes; ++m) {
    dryvolmr(:,:) = 0._wp

    // Get mode properties and aerosol prognostics.
    Real sigmag = modes[m].mean_std_dev;
    Real rhcrystal = modes[m].rh_crystal;
    Real rhdeliques = modes[m].rh_deliques;

    // Get species information for this mode.
    const std::vector<Species>& species = model.aerosol_species_for_mode(m);
    int num_species = species.size();

    // Interstitial aerosol mixing fractions.
    const auto& int_aero_mfrac = prognostics.interstitial_aerosols(m);

    // Loop over species
    for (int s = 0; s < num_species; ++s) {

      // Get species properties.
      Real specdens = species.mass_density;
      Real spechygro = species.hygroscopicity;

      // get species interstitial mixing ratio ('a')

      // get species mass mixing ratio
      call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, raer)

      if (l == 0) {
        // save off these values to be used as defaults
        specdens_1(m)  = specdens
        spechygro_1    = spechygro
      }

      for (int i = 0; i < num_columns; ++i) {
        for (int k = 0; k < num_levels; ++k) {
          Real raer_ik = raer(i,k);
          maer(i,k,m) += raeri_ik;
          Real = duma/specdens
          dryvolmr(i,k) += dryvolmr(i,k) + 
          hygro(i,k,m)  = hygro(i,k,m) + dumb*spechygro
        } // levels
      } // columns
    } // species

    Real alnsg = log(sigmag);

    for (int i = 0; i < num_columns; ++i) {
      for (int k = 0; k < num_levels; ++k) {

        if (dryvolmr(i,k) > 1.0e-30_wp) {
          // volume-mean hygroscopicity
          hygro(i,k,m) = hygro(i,k,m)/dryvolmr(i,k);
        } else {
          hygro(i,k,m) = spechygro_1;
        }

        // dry aerosol properties
        v2ncur_a = 1._wp / ( (pi/6._wp)*(dgncur_a(i,k,m)**3._wp)*exp(4.5_wp*alnsg**2._wp) )

        // naer = aerosol number (#/kg)
        naer(i,k,m) = dryvolmr(i,k)*v2ncur_a

        // compute mean (1 particle) dry volume and mass for each mode
        // old coding is replaced because the new (1/v2ncur_a) is equal to
        // the mean particle volume
        // also moletomass forces maer >= 1.0e-30, so (maer/dryvolmr)
        // should never cause problems (but check for maer < 1.0e-31 anyway)
        if (maer(i,k,m) > 1.0e-31) {
          drydens = maer(i,k,m)/dryvolmr(i,k)
        } else {
          drydens = 1.0_wp
        }
        dryvol(i,k,m)   = 1.0_wp/v2ncur_a
        drymass(i,k,m)  = drydens*dryvol(i,k,m)
        dryrad(i,k,m)   = (dryvol(i,k,m)/pi43)**third
      }
    }
  }

  // Compute the clear-sky relative humidity.
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
*/

  // Call the Fortran water uptake subroutine.
  int top_lev = 1; // top of atmosphere is the first entry in the column
  bool bisect = use_bisection_;
  modal_aero_wateruptake_sub(&num_columns, &top_lev, &num_levels, &num_modes,
                             &bisect, &rhcrystal[0], &rhdeliques[0], dryrad,
                             naer, hygro, rh, dryvol, drymass, specdens_1,
                             &(*dgncur_a.data())[0],
                             &(*dgncur_awet.data())[0],
                             &(*qaerwat.data())[0],
                             &(*wetdens.data())[0]);

  // Accumulate aerosol water over modes.
  for (int m = 0; m < num_modes; ++m) {
    for (int i = 0; i < num_columns; ++i) {
      for (int k = 0; k < num_levels; ++k) {
        aerosol_water(i, k) += qaerwat(m, i, k);
      }
    }
  }
}

}

