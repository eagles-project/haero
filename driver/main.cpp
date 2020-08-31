// This is Haero's standalone C++ driver. It reads a YAML file containing input
// parameters for a simulation (or an ensemble of simulations) and dispatches
// the runs to a Fortran driver subroutine.

#include <cstdio>
#include <cstdlib>
#include <limits>
#include "haero/yaml_file.hpp"
#include "driver/driver.hpp"

using namespace std;
using namespace haero;

namespace {

// Global simulation data. Needed for providing input to the Fortran portion
// of the driver. We can scuttle this when we've moved the bulk of the driver
// to C++.
Sim_input_data* _sim_data = nullptr;

// Are we a member of an ensemble? If so, this defines a unique index for the
// active ensemble member.
int _member_index = -1;

// Are we conducting a convergence study using multiple time step sizes? If
// so, this identifies a unique index for setting the time step.
int _dt_index = -1;

// Print driver usage information and exit.
void usage(const char* exe)
{
  fprintf(stderr, "%s: usage:\n", exe);
  fprintf(stderr, "%s <input.yml>\n", exe);
  exit(1);
}

} // anonymous namespace

extern "C" {

// Here's a function that delivers input to the current member of an ensemble
// in the Fortran driver.
void read_cxx_yaml(int* mam_dt, int* mam_nstep,
                   int* do_gaschem, int* do_cloudchem, int* do_gasaerexch,
                   int* do_rename, int* do_newnuc, int* do_coag,
                   int* do_calcsize, int* num_unit, int* frac_unit,
                   int* gas_unit, Real* temp, Real* press, Real* RH_CLEA,
                   Real* hgt, Real* cld_frac, Real* numc1, Real* numc2,
                   Real* numc3, Real* numc4, Real* mfso41, Real* mfpom1,
                   Real* mfsoa1, Real* mfbc1, Real* mfdst1, Real* mfncl1,
                   Real* mfso42, Real* mfsoa2, Real* mfncl2, Real* mfdst3,
                   Real* mfncl3, Real* mfso43, Real* mfbc3, Real* mfpom3,
                   Real* mfsoa3, Real* mfpom4, Real* mfbc4, Real* qso2,
                   Real* qh2so4, Real* qsoag, Real* num_factor,
                   Real* gas_factor, Real* h2so4_chem_prod_rate,
                   int* mam_outfrq, char output_file[255])
{
  // Read in parameters.
  *mam_dt = _sim_data->dts[_dt_index];
  *mam_nstep = int(_sim_data->duration / (*mam_dt));
  *do_gaschem = _sim_data->do_gaschem;
  *do_cloudchem = _sim_data->do_cloudchem;
  *do_gasaerexch = _sim_data->do_gasaerexch;
  *do_rename = _sim_data->do_rename;
  *do_newnuc = _sim_data->do_newnuc;
  *do_coag = _sim_data->do_coag;
  *do_calcsize = _sim_data->do_calcsize;
  *num_unit = _sim_data->num_unit;
  *frac_unit = _sim_data->frac_unit;
  *gas_unit = _sim_data->gas_unit;
  *temp = _sim_data->temp;
  *press = _sim_data->press;
  *RH_CLEA = _sim_data->RH_CLEA;
  *hgt = _sim_data->hgt;
  *cld_frac = _sim_data->cld_frac;

  // Pick the right set of initial conditions based on our ensemble index.
  // NOTE: We use a non-const reference for ics because of the idiotic behavior
  // of the STL's map bracket operator, which creates an entry if it's not
  // found.
  map<string, Real>& ics = _sim_data->initial_conditions[_member_index];

  // Read in initial conditions (values default to 0 if not found).
  *numc1 = ics["numc1"];
  *numc2 = ics["numc2"];
  *numc3 = ics["numc3"];
  *numc4 = ics["numc4"];
  *mfso41 = ics["mfso41"];
  *mfpom1 = ics["mfpom1"];
  *mfsoa1 = ics["mfsoa1"];
  *mfbc1 = ics["mfbc1"];
  *mfdst1 = ics["mfdst1"];
  *mfncl1 = ics["mfncl1"];
  *mfso42 = ics["mfso42"];
  *mfsoa2 = ics["mfsoa2"];
  *mfncl2 = ics["mfncl2"];
  *mfdst3 = ics["mfdst3"];
  *mfncl3 = ics["mfncl3"];
  *mfso43 = ics["mfso43"];
  *mfbc3 = ics["mfbc3"];
  *mfpom3 = ics["mfpom3"];
  *mfsoa3 = ics["mfsoa3"];
  *mfpom4 = ics["mfpom4"];
  *mfbc4 = ics["mfbc4"];
  *qso2 = ics["qso2"];
  *qh2so4 = ics["qh2so4"];
  *qsoag = ics["qsoag"];
  *num_factor = ics["num_factor"];
  *gas_factor = ics["gas_factor"];
  *h2so4_chem_prod_rate = ics["h2so4_chem_prod_rate"];

  // Apply perturbations if needed.
  if (_sim_data->perturb_factor != 0.0)
  {
    constexpr Real eps = numeric_limits<Real>::epsilon();
    Real pert = _sim_data->perturb_factor * eps;
    *numc1 += pert;
    *numc2 += pert;
    *numc3 += pert;
    *numc4 += pert;
    *mfso41 += pert;
    *mfpom1 += pert;
    *mfsoa1 += pert;
    *mfbc1 += pert;
    *mfdst1 += pert;
    *mfncl1 += pert;
    *mfso42 += pert;
    *mfsoa2 += pert;
    *mfncl2 += pert;
    *mfdst3 += pert;
    *mfncl3 += pert;
    *mfso43 += pert;
    *mfbc3 += pert;
    *mfpom3 += pert;
    *mfsoa3 += pert;
    *mfpom4 += pert;
    *mfbc4 += pert;
    *qso2 += pert;
    *qh2so4 += pert;
    *qsoag += pert;
    *num_factor += pert;
    *gas_factor += pert;
    *h2so4_chem_prod_rate += pert;
  }

  // Determine the output file name. If we're a member of an ensemble, the
  // output file is labeled with the ensemble index. We also add a label for
  // output files that are associated with convergence studies.
  string tag;
  if (_sim_data->dts.size() > 1)
    tag += string("_dt=") + to_string(*mam_dt);
  if (_sim_data->initial_conditions.size() > 1)
    tag += string("_") + to_string(_member_index+1);
  snprintf(output_file, 254, "%s/%s%s.nc", _sim_data->output_dir.c_str(),
           _sim_data->output_prefix.c_str(), tag.c_str());

  // Get the output frequency. Perversely, we set it to 1 by default.
  *mam_outfrq = _sim_data->output_freq;
  if (*mam_outfrq < 0)
    *mam_outfrq = 1;
}

}

int main(int argc, const char* argv[])
{
  if (argc < 2)
    usage(argv[0]);

  // Read the input file.
  int ensemble_size = 0;
  int num_timesteps = 0;
  try
  {
    // Read our Yaml_file.
    Yaml_file input(argv[1]);

    // Populate our global ensemble data.
    _sim_data = new Sim_input_data;

    _sim_data->initial_conditions = input.read_initial_conditions();
    ensemble_size = _sim_data->initial_conditions.size();
    _sim_data->perturb_factor = input.read_perturb_factor();

    input.read_control_params(_sim_data->do_gaschem,
                              _sim_data->do_cloudchem,
                              _sim_data->do_gasaerexch,
                              _sim_data->do_rename,
                              _sim_data->do_newnuc,
                              _sim_data->do_coag,
                              _sim_data->do_calcsize,
                              _sim_data->num_unit,
                              _sim_data->frac_unit,
                              _sim_data->gas_unit);

    _sim_data->dts = input.read_timesteps();
    num_timesteps = _sim_data->dts.size();
    _sim_data->duration = input.read_duration();

    input.read_met_input(_sim_data->temp, _sim_data->press,
                         _sim_data->RH_CLEA, _sim_data->hgt,
                         _sim_data->cld_frac);

    input.read_output_params(_sim_data->output_dir,
                             _sim_data->output_prefix,
                             _sim_data->output_freq);
  }
  catch (Yaml_exception& e)
  {
    fprintf(stderr, "Input error: %s\n", e.what());
    exit(1);
  }
//  catch (std::exception& e)
//  {
//    fprintf(stderr, "Error: %s\n", e.what());
//    exit(1);
//  }

  // Call the Fortran driver with each set of initial conditions. Here, we
  // loop over two indices:
  // 1. _dt_index is a global timestepping index that determines which timestep
  //    is used in a simulation or simulation ensemble.
  // 2. _member_index is a global ensemble member index that determines which
  // initial data is parceled to the Fortran driver via read_cxx_yaml(), above.
  // This is crummy, but it's a stepping stone toward an ensemble-aware C++
  // driver.
  for (_dt_index = 0; _dt_index < num_timesteps; ++_dt_index)
    for (_member_index = 0; _member_index < ensemble_size; ++_member_index)
      haero_driver(*_sim_data);

  // Clean up.
  delete _sim_data;

  return 0;
}
