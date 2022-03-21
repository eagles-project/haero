#ifndef HAERO_BOX_OUTPUT_HPP
#define HAERO_BOX_OUTPUT_HPP

#include <haero/haero.hpp>
#include <haero/diagnostics.hpp>
#include <haero/modal_aerosol_config.hpp>
#include <haero/prognostics.hpp>
#include <haero/tendencies.hpp>

#include <map>
#include <string>

namespace Box {

using Real = haero::Real;

class BoxOutput {
  public:

  // Constructs a BoxOutput for the given aerosol configuration.
  explicit BoxOutput(const haero::ModalAerosolConfig& config);

  // Appends the given prognostic and diagnostic data to that maintained by
  // this BoxOutput. The tendencies are indexed by their process names, since
  // contributions from different processes are tracked separately.
  void append(const haero::Prognostics& prognostics,
              const haero::HostDiagnostics& diagnostics,
              const std::map<std::string, haero::Tendencies>& tendencies);

  // Writes NetCDF data to the file with the given name.
  void write(const std::string& filename) const;

  private:

  const haero::ModalAerosolConfig& config_; // aerosol configuration

  std::vector<int> iso4_, isoa_; // indices for SO4 and SOA within each mode
  int ih2so4_, isoag_;           // indices for H2SO4 and SOAG gases

  int nstep_; // number of steps written

  // Aerosol data (identical in form to that of the MAM box model)
  std::vector<Real> num_aer_, so4_aer_, soa_aer_, h2so4_, soag_,
                    dgn_a_, dgn_awet_;

  // Tendency data (also identical to that of the MAM box model)
  std::vector<Real> qtend_cond_aging_so4_, qtend_rename_so4_,
                    qtend_newnuc_so4_, qtend_coag_so4_;
  std::vector<Real> qtend_cond_aging_soa_, qtend_rename_soa_,
                    qtend_newnuc_soa_, qtend_coag_soa_;
  std::vector<Real> qtend_cond_aging_h2so4_, qtend_rename_h2so4_,
                    qtend_newnuc_h2so4_, qtend_coag_h2so4_;
  std::vector<Real> qtend_cond_aging_soag_, qtend_rename_soag_,
                    qtend_newnuc_soag_, qtend_coag_soag_;
};

} // namespace Box

#endif

