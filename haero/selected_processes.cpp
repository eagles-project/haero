#include "haero/process.hpp"
#include "haero/selected_processes.hpp"
#include "haero/available_processes.hpp"

namespace haero {

SelectedProcesses::SelectedProcesses():
  activation(NoActivation),
  cloudborne_wet_removal(NoCloudBorneWetRemoval),
  coagulation(NoCoagulation),
  condensation(NoCondensation),
  dry_deposition(NoDryDeposition),
  emissions(NoEmissions),
  interstitial_wet_removal(NoInterstitialWetRemoval),
  nucleation(NoNucleation),
  resuspension(NoResuspension),
  water_uptake(NoWaterUptake) {
}

PrognosticProcess* select_prognostic_process(ProcessType type,
                                             const SelectedProcesses& selections) {
  PrognosticProcess* process = nullptr;
  if (type == ActivationProcess) {
    if (selections.activation == SelectedProcesses::NoActivation) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == CloudBorneWetRemovalProcess) {
    if (selections.cloudborne_wet_removal == SelectedProcesses::NoCloudBorneWetRemoval) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == CoagulationProcess) {
    if (selections.coagulation == SelectedProcesses::NoCoagulation) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == CondensationProcess) {
    if (selections.condensation == SelectedProcesses::NoCondensation) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == DryDepositionProcess) {
    if (selections.cloudborne_wet_removal == SelectedProcesses::NoCloudBorneWetRemoval) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == EmissionsProcess) {
    if (selections.emissions == SelectedProcesses::NoEmissions) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == InterstitialWetRemovalProcess) {
    if (selections.interstitial_wet_removal == SelectedProcesses::NoInterstitialWetRemoval) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == NucleationProcess) {
    if (selections.nucleation == SelectedProcesses::MAMNucleation) {
      process = new MAMNucleationProcess();
#if HAERO_FORTRAN
    } else if (selections.nucleation == SelectedProcesses::MAMFNucleation) {
      process = new MAMNucleationFProcess();
#endif
    } else if (selections.nucleation == SelectedProcesses::NoNucleation) {
      process = new NullPrognosticProcess(type);
    }
  }
  else if (type == ResuspensionProcess) {
    if (selections.resuspension == SelectedProcesses::NoResuspension) {
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
                                             const SelectedProcesses& selections) {
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
    if (selections.water_uptake == SelectedProcesses::NoWaterUptake) {
      process = new NullDiagnosticProcess(type);
    }
  }
  EKAT_REQUIRE_MSG(process != nullptr, "No diagnostic process was selected!");
  return process;
}

}

