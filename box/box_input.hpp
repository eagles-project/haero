#ifndef HAERO_BOX_INPUT_HPP
#define HAERO_BOX_INPUT_HPP

#include <haero/haero_config.hpp>

#include <string>
#include <iostream>

extern "C" {

// This struct holds input data read from a namelist. The namelist is
// identical to that used by the legacy MAM box model.
typedef struct BoxInput {
  // time_input
  int         mam_dt, mam_nstep, mam_output_intvl;

  // cntl_input
  int         mdo_gaschem, mdo_cloudchem, mdo_gasaerexch, mdo_test_siamgs19,
              mdo_rename, mdo_newnuc, mdo_coag, mdo_calcsize,
              newnuc_method_user_choice, pbl_nuc_wang2008_user_choice;
  bool        cluster_growth_arg_write;
  char        cluster_growth_arg_write_fname[128];
  int         num_unit, frac_unit, gas_unit;
  haero::Real ic_perturb_factor;

  // met_input
  haero::Real temp, press, RH_CLEA, hgt, cld_frac;

  // chem_input
  haero::Real numc1, numc2, numc3, numc4,
              mfso41, mfpom1, mfsoa1, mfbc1, mfdst1, mfncl1,
              mfso42, mfsoa2, mfncl2,
              mfdst3, mfncl3, mfso43, mfbc3, mfpom3, mfsoa3,
              mfpom4, mfbc4,
              qso2, qh2so4, qsoag, num_factor, gas_factor,
              h2so4_chem_prod_rate;
} BoxInput;

} // extern "C"

namespace Box {

// Reads input from a Fortran namelist with the given filename.
BoxInput read_namelist(const std::string& filename);

// Use this to output a BoxInput struct to check namelist data.
std::ostream& operator<<(std::ostream& s, const BoxInput& input);

} // namespace Box

#endif

