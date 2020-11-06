#include "haero/process.hpp"
#include "haero/parameterizations.hpp"

namespace haero {

PrognosticProcess* select_prognostic_process(ProcessType type,
                                             const Parameterizations& selections)
{
  PrognosticProcess* process = nullptr;
  if (type == ActivationProcess) {
    if (selections.activation == Parameterizations::NoActivation) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == CloudBorneWetRemovalProcess) {
    if (selections.cloudborne_wet_removal == Parameterizations::NoCloudBorneWetRemoval) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == CoagulationProcess) {
    if (selections.coagulation == Parameterizations::NoCoagulation) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == CondensationProcess) {
    if (selections.condensation == Parameterizations::NoCondensation) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == DryDepositionProcess) {
    if (selections.cloudborne_wet_removal == Parameterizations::NoCloudBorneWetRemoval) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == EmissionsProcess) {
    if (selections.emissions == Parameterizations::NoEmissions) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == InterstitialWetRemovalProcess) {
    if (selections.interstitial_wet_removal == Parameterizations::NoInterstitialWetRemoval) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == NucleationProcess) {
    if (selections.nucleation == Parameterizations::NoNucleation) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == ResuspensionProcess) {
    if (selections.resuspension == Parameterizations::NoResuspension) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == WaterUptakeProcess) {
    EKAT_REQUIRE_MSG(false, "Water uptake is a diagnostic aerosol process!");
  }
  EKAT_REQUIRE_MSG(process != nullptr, "No prognostic process was selected!");
  return process;
}

DiagnosticProcess* select_diagnostic_process(ProcessType type,
                                             const Parameterizations& selections)
{
  DiagnosticProcess* process = nullptr;
  if (type == ActivationProcess) {
    EKAT_REQUIRE_MSG(process != nullptr, "Activation is a prognostic aerosol process!");
  }
  else if (type == CloudBorneWetRemovalProcess) {
    EKAT_REQUIRE_MSG(process != nullptr, "Cloud-borne wet removal is a prognostic aerosol process!");
  }
  else if (type == CoagulationProcess) {
    EKAT_REQUIRE_MSG(process != nullptr, "Coagulation is a prognostic aerosol process!");
  }
  else if (type == CondensationProcess) {
    EKAT_REQUIRE_MSG(process != nullptr, "Condensation is a prognostic aerosol process!");
  }
  else if (type == DryDepositionProcess) {
    EKAT_REQUIRE_MSG(process != nullptr, "Dry deposition is a prognostic aerosol process!");
  }
  else if (type == EmissionsProcess) {
    EKAT_REQUIRE_MSG(process != nullptr, "Emissions are a prognostic aerosol process!");
  }
  else if (type == InterstitialWetRemovalProcess) {
    EKAT_REQUIRE_MSG(process != nullptr, "Interstitial wet removal is a prognostic aerosol process!");
  }
  else if (type == NucleationProcess) {
    EKAT_REQUIRE_MSG(process != nullptr, "Nucleation is a prognostic aerosol process!");
  }
  else if (type == ResuspensionProcess) {
    EKAT_REQUIRE_MSG(process != nullptr, "Resuspension is a prognostic aerosol process!");
  }
  else if (type == WaterUptakeProcess) {
    if (selections.water_uptake == Parameterizations::NoWaterUptake) {
      process = new NullDiagnosticProcess(type);
    }
  }
  EKAT_REQUIRE_MSG(process != nullptr, "No diagnostic process was selected!");
  return process;
}

}
