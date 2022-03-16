#ifndef HAERO_BOX_OUTPUT_HPP
#define HAERO_BOX_OUTPUT_HPP

#include <haero/haero.hpp>
#include <haero/diagnostics.hpp>
#include <haero/prognostics.hpp>
#include <haero/modal_aerosol_config.hpp>

#include <string>

namespace Box {

class BoxOutput {
  public:

  // Constructs a BoxOutput for the given aerosol configuration.
  explicit BoxOutput(const haero::ModalAerosolConfig& config);

  // Appends the given prognostic and diagnostic data to that maintained by
  // this BoxOutput.
  void append(const haero::Prognostics& prognostics,
              const haero::HostDiagnostics& diagnostics);

  // Writes NetCDF data to the file with the given name.
  void write(const std::string& filename) const;

  private:

  ModalAerosolConfig& config_;

};

} // namespace Box

#endif

