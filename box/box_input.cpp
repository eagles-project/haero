#include "box_input.hpp"

#include <cstring>

extern "C" {

// This function reaches into the Fortran universe and extracts a BoxInput's
// worth of input parameters.
void read_nml(char filename[255], BoxInput *input);

} // extern "C"

namespace Box {

BoxInput read_namelist(const std::string& filename) {
  BoxInput input;
  char nml_file[255];
  std::strncpy(nml_file, filename.c_str(), 255);
  read_nml(nml_file, &input);
  return input;
}

std::ostream& operator<<(std::ostream& s, const BoxInput& input) {
  s << "&time_input\n"
    << "mam_dt           = " << input.mam_dt << ",\n"
    << "mam_nstep        = " << input.mam_nstep << ",\n"
    << "mam_output_intvl = " << input.mam_output_intvl << ",\n"
    << "/\n"
    << "&cntl_input\n"
    << "mdo_calcsize                 = " << input.mdo_calcsize << ",\n"
    << "mdo_gaschem                  = " << input.mdo_gaschem << ",\n"
    << "mdo_cloudchem                = " << input.mdo_cloudchem << ",\n"
    << "mdo_gasaerexch               = " << input.mdo_gasaerexch << ",\n"
    << "mdo_rename                   = " << input.mdo_rename << ",\n"
    << "mdo_newnuc                   = " << input.mdo_newnuc << ",\n"
    << "newnuc_method_user_choice    = " << input.newnuc_method_user_choice << ",\n"
    << "pbl_nuc_wang2008_user_choice = " << input.pbl_nuc_wang2008_user_choice << ",\n"
    << "mdo_coag                     = " << input.mdo_coag << ",\n"
    << "mdo_test_siamgs19            = " << input.mdo_test_siamgs19 << ",\n"
    << "num_unit                     = " << input.num_unit << ",\n"
    << "frac_unit                    = " << input.frac_unit << ",\n"
    << "gas_unit                     = " << input.gas_unit << ",\n"
    << "ic_perturb_factor            = " << input.ic_perturb_factor << ",\n"
    << "/\n"
    << "&met_input\n"
    << "temp     = " << input.temp << ",\n"
    << "press    = " << input.press << ",\n"
    << "RH_CLEA  = " << input.RH_CLEA << ",\n"
    << "hgt      = " << input.hgt << ",\n"
    << "cld_frac = " << input.cld_frac << ",\n"
    << "/\n"
    << "&chem_input\n"
    << "!\n"
    << "h2so4_chem_prod_rate = " << input.h2so4_chem_prod_rate << ",\n"
    << "!\n"
    << "! initial number concentrations\n"
    << "!\n"
    << "numc1 = " << input.numc1 << ",\n"
    << "numc2 = " << input.numc2 << ",\n"
    << "numc3 = " << input.numc3 << ",\n"
    << "numc4 = " << input.numc4 << ",\n"
    << "!\n"
    << "! mfABCx: mass fraction of species ABC in mode x.\n"
    << "!\n"
    << "! The mass fraction of mom is calculated by\n"
    << "! 1 - sum(mfABCx). If sum(mfABCx) > 1, an error\n"
    << "! is issued by the test driver. number of species\n"
    << "! ABC in each mode x comes from the MAM4 with mom.\n"
    << "!\n"
    << "mfso41     = " << input.mfso41 << ",\n"
    << "mfpom1     = " << input.mfpom1 << ",\n"
    << "mfsoa1     = " << input.mfsoa1 << ",\n"
    << "mfbc1      = " << input.mfbc1 << ",\n"
    << "mfdst1     = " << input.mfdst1 << ",\n"
    << "mfncl1     = " << input.mfncl1 << ",\n"
    << "mfso42     = " << input.mfso42 << ",\n"
    << "mfsoa2     = " << input.mfsoa2 << ",\n"
    << "mfncl2     = " << input.mfncl2 << ",\n"
    << "mfdst3     = " << input.mfdst3 << ",\n"
    << "mfncl3     = " << input.mfncl3 << ",\n"
    << "mfso43     = " << input.mfso43 << ",\n"
    << "mfbc3      = " << input.mfbc3 << ",\n"
    << "mfpom3     = " << input.mfpom3 << ",\n"
    << "mfsoa3     = " << input.mfsoa3 << ",\n"
    << "mfpom4     = " << input.mfpom4 << ",\n"
    << "mfbc4      = " << input.mfbc4 << ",\n"
    << "qso2       = " << input.qso2 << ",\n"
    << "qh2so4     = " << input.qh2so4 << ",\n"
    << "qsoag      = " << input.qsoag << ",\n"
    << "num_factor = " << input.num_factor << ",\n"
    << "gas_factor = " << input.gas_factor << ",\n"
    << "/\n";
  return s;
}

} // namespace Box

