#include "box_output.hpp"

#include <time.h>

#if HAERO_DOUBLE_PRECISION
#define NC_REAL NC_DOUBLE
#define nc_put_var_real nc_put_var_double
#else
#define NC_REAL NC_FLOAT
#define nc_put_var_real nc_put_var_float
#endif

namespace Box {

BoxOutput::BoxOutput(const haero::ModalAerosolConfig& config,
                     int output_interval):
  config_(config),
  output_interval_(output_interval),
  n_steps(0) {
}

void BoxOutput::append(const haero::Prognostics& prognostics,
                       const haero::HostDiagnostics& diagnostics,
                       const std::vector<haero::Tendencies>& tendencies) {
  size_t num_modes = static_cast<size_t>(config_.num_modes());
  size_t old_size = num_aer_.size();
  size_t new_size = old_size + num_modes;

  // Resize the containers.
  num_aer_.resize(new_size);
  so4_aer_.resize(new_size);
  soa_aer_.resize(new_size);
  h2so4_.resize(new_size);
  soag_.resize(new_size);
  dgn_a_.resize(new_size);
  dgn_awet_.resize(new_size);

  qtend_cond_aging_so4_.resize(new_size);
  qtend_rename_so4_.resize(new_size);
  qtend_newnuc_so4_.resize(new_size);
  qtend_coag_so4_.resize(new_size);
  qtend_cond_aging_soa_.resize(new_size);
  qtend_rename_soa_.resize(new_size);
  qtend_newnuc_soa_.resize(new_size);
  qtend_coag_soa_.resize(new_size);
  qtend_cond_aging_h2so4_.resize(new_size);
  qtend_rename_h2so4_.resize(new_size);
  qtend_newnuc_h2so4_.resize(new_size);
  qtend_coag_h2so4_.resize(new_size);
  qtend_cond_aging_soag_.resize(new_size);
  qtend_rename_soag_.resize(new_size);
  qtend_newnuc_soag_.resize(new_size);
  qtend_coag_soag_.resize(new_size);

  // Extract the new data.

  ++nstep_;
}

void BoxOutput::write(const std::string& filename) const {
  int err;

  int ncid;
  err = nc_create(filename, NC_CLOBBER, &ncid); // returns NC_NOERR, hopefully

  // Define dimensions.
  int nstep_dimid, mode_dimid;
  err = nc_def_dim(ncid, "nsteps", nstep_, &nstep_dimid);
  err = nc_def_dim(ncid, "mode", config_.num_modes(), &mode_dimid);

  // The dimids array is used to pass the IDs of the dimensions of
  // the variables.
  int dimids[2] = {nstep_dimid, mode_dimid};

  // Define variables.
  int varid[23];
  err = nc_def_var(ncid, "num_aer", NC_REAL, dimids, &varid[0]);
  err = nc_def_var(ncid, "so4_aer", NC_REAL, dimids, &varid[1]);
  err = nc_def_var(ncid, "soa_aer", NC_REAL, dimids, &varid[2]);
  err = nc_def_var(ncid, "h2so4_gas", NC_REAL, nstep_dimid, &varid[3]);
  err = nc_def_var(ncid, "soag_gas", NC_REAL, nstep_dimid, &varid[4]);
  err = nc_def_var(ncid, "dgn_a", NC_REAL, dimids, &varid[5]);
  err = nc_def_var(ncid, "dgn_awet", NC_REAL, dimids, &varid[6]);
  err = nc_def_var(ncid, "qtend_cond_aging_so4", NC_REAL, dimids, &varid[7]);
  err = nc_def_var(ncid, "qtend_rename_so4", NC_REAL, dimids, &varid[8]);
  err = nc_def_var(ncid, "qtend_newnuc_so4", NC_REAL, dimids, &varid[9]);
  err = nc_def_var(ncid, "qtend_coag_so4", NC_REAL, dimids, &varid[10]);
  err = nc_def_var(ncid, "qtend_cond_aging_soa", NC_REAL, dimids, &varid[11]);
  err = nc_def_var(ncid, "qtend_rename_soa", NC_REAL, dimids, &varid[12]);
  err = nc_def_var(ncid, "qtend_newnuc_soa", NC_REAL, dimids, &varid[13]);
  err = nc_def_var(ncid, "qtend_coag_soa", NC_REAL, dimids, &varid[14]);
  err = nc_def_var(ncid, "qtend_cond_aging_h2so4", NC_REAL, nstep_dimid, &varid[15]);
  err = nc_def_var(ncid, "qtend_rename_h2so4", NC_REAL, nstep_dimid, &varid[16]);
  err = nc_def_var(ncid, "qtend_newnuc_h2so4", NC_REAL, nstep_dimid, &varid[17]);
  err = nc_def_var(ncid, "qtend_coag_h2so4", NC_REAL, nstep_dimid, &varid[18]);
  err = nc_def_var(ncid, "qtend_cond_aging_soag", NC_REAL, nstep_dimid, &varid[19]);
  err = nc_def_var(ncid, "qtend_rename_soag", NC_REAL, nstep_dimid, &varid[20]);
  err = nc_def_var(ncid, "qtend_newnuc_soag", NC_REAL, nstep_dimid, &varid[21]);
  err = nc_def_var(ncid, "qtend_coag_soag", NC_REAL, nstep_dimid, &varid[22]);

  // Assign attributes to fields.
  err = nc_put_att(ncid, varid[0], "units", "#/kg-air");
  err = nc_put_att(ncid, varid[1], "units", "kg-aer/kg-air");
  err = nc_put_att(ncid, varid[2], "units", "kg-aer/kg-air");
  err = nc_put_att(ncid, varid[3], "units", "kg-gas/kg-air");
  err = nc_put_att(ncid, varid[4], "units", "kg-gas/kg-air");
  err = nc_put_att(ncid, varid[5], "units", "meter");
  err = nc_put_att(ncid, varid[6], "units", "meter");
  err = nc_put_att(ncid, varid[7], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[7], "descr", "condensation-aging tendency");
  err = nc_put_att(ncid, varid[8], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[8], "descr", "rename tendency");
  err = nc_put_att(ncid, varid[9], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[9], "descr", "nucleation tendency");
  err = nc_put_att(ncid, varid[10], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[10], "descr", "coagulation tendency");
  err = nc_put_att(ncid, varid[11], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[11], "descr", "condensation-aging tendency");
  err = nc_put_att(ncid, varid[12], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[12], "descr", "rename tendency");
  err = nc_put_att(ncid, varid[13], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[13], "descr", "nucleation tendency");
  err = nc_put_att(ncid, varid[14], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[14], "descr", "coagulation tendency");
  err = nc_put_att(ncid, varid[15], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[15], "descr", "condensation-aging tendency");
  err = nc_put_att(ncid, varid[16], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[16], "descr", "rename tendency");
  err = nc_put_att(ncid, varid[17], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[17], "descr", "nucleation tendency");
  err = nc_put_att(ncid, varid[18], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[18], "descr", "coagulation tendency");
  err = nc_put_att(ncid, varid[19], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[19], "descr", "condensation-aging tendency");
  err = nc_put_att(ncid, varid[20], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[20], "descr", "rename tendency");
  err = nc_put_att(ncid, varid[21], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[21], "descr", "nucleation tendency");
  err = nc_put_att(ncid, varid[22], "units", "mol mol-1 s-1");
  err = nc_put_att(ncid, varid[22], "descr", "coagulation tendency");

  // Add global attributes.
  nc_put_att(ncid, NC_GLOBAL, & "Created_by", "HAERO");
  time_t current_time = time(NULL);
  const char *date = ctime(&current_time);
  err = nc_put_att(ncid, NC_GLOBAL, "Created_date", date);

  // End define mode. This tells netCDF we are done defining metadata.
  err = nc_enddef(ncid);

  // Write variable data.
  err = nc_put_var(ncid, varid[0], &num_aer_[0]);
  err = nc_put_var(ncid, varid[1], &so4_aer_[0]);
  err = nc_put_var(ncid, varid[2], &soa_aer_[0]);
  err = nc_put_var(ncid, varid[3], &h2so4_[0]);
  err = nc_put_var(ncid, varid[4], &soag_[0]);
  err = nc_put_var(ncid, varid[5], &dgn_a_[0]);
  err = nc_put_var(ncid, varid[6], &dgn_awet_[0]);

  err = nc_put_var(ncid, varid[7],  &qtend_cond_aging_so4_[0]);
  err = nc_put_var(ncid, varid[8],  &qtend_rename_so4_[0]);
  err = nc_put_var(ncid, varid[9],  &qtend_newnuc_so4_[0]);
  err = nc_put_var(ncid, varid[10], &qtend_coag_so4_[0]);
  err = nc_put_var(ncid, varid[11], &qtend_cond_aging_soa_[0]);
  err = nc_put_var(ncid, varid[12], &qtend_rename_soa_[0]);
  err = nc_put_var(ncid, varid[13], &qtend_newnuc_soa_[0]);
  err = nc_put_var(ncid, varid[14], &qtend_coag_soa_[0]);
  err = nc_put_var(ncid, varid[15], &qtend_cond_aging_h2so4_[0]);
  err = nc_put_var(ncid, varid[16], &qtend_rename_h2so4_[0]);
  err = nc_put_var(ncid, varid[17], &qtend_newnuc_h2so4_[0]);
  err = nc_put_var(ncid, varid[18], &qtend_coag_h2so4_[0]);
  err = nc_put_var(ncid, varid[19], &qtend_cond_aging_soag_[0]);
  err = nc_put_var(ncid, varid[20], &qtend_rename_soag_[0]);
  err = nc_put_var(ncid, varid[21], &qtend_newnuc_soag_[0]);
  err = nc_put_var(ncid, varid[22], &qtend_coag_soag_[0]);

  // Close the file. This frees up any internal netCDF resources
  // associated with the file, and flushes any buffers.
  err = nc_close(ncid);
}

} // namespace Box


