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

  BoxOutput(const haero::ModalAerosolConfig& config);

  void append(const haero::Prognostics& prognostics,
              const haero::HostDiagnostics& diagnostics);

  void write(const std::string& filename) const;
};

} // namespace Box

#endif

