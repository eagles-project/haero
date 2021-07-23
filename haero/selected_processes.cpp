#include "haero/selected_processes.hpp"

#include "haero/available_processes.hpp"

namespace haero {

SelectedProcesses::SelectedProcesses()
    : activation(NoActivation),
      cloudborne_wet_removal(NoCloudBorneWetRemoval),
      coagulation(NoCoagulation),
      condensation(NoCondensation),
      dry_deposition(NoDryDeposition),
      emissions(NoEmissions),
      interstitial_wet_removal(NoInterstitialWetRemoval),
      nucleation(NoNucleation),
      calcsize(NoCalcsize),
      rename(NoRename),
      resuspension(NoResuspension) {}

AerosolProcess* select_aerosol_process(AerosolProcessType type,
                                       const SelectedProcesses& selections) {
  AerosolProcess* process = nullptr;
  if (type == ActivationProcess) {
    if (selections.activation == SelectedProcesses::NoActivation) {
      process = new NullAerosolProcess(type);
    }
  } else if (type == CloudBorneWetRemovalProcess) {
    if (selections.cloudborne_wet_removal ==
        SelectedProcesses::NoCloudBorneWetRemoval) {
      process = new NullAerosolProcess(type);
    }
  } else if (type == CoagulationProcess) {
    if (selections.coagulation == SelectedProcesses::NoCoagulation) {
      process = new NullAerosolProcess(type);
    }
  } else if (type == CondensationProcess) {
    if (selections.condensation == SelectedProcesses::NoCondensation) {
      process = new NullAerosolProcess(type);
    }
  } else if (type == DryDepositionProcess) {
    if (selections.cloudborne_wet_removal ==
        SelectedProcesses::NoCloudBorneWetRemoval) {
      process = new NullAerosolProcess(type);
    }
  } else if (type == EmissionsProcess) {
    if (selections.emissions == SelectedProcesses::NoEmissions) {
      process = new NullAerosolProcess(type);
    }
  } else if (type == InterstitialWetRemovalProcess) {
    if (selections.interstitial_wet_removal ==
        SelectedProcesses::NoInterstitialWetRemoval) {
      process = new NullAerosolProcess(type);
    }
  } else if (type == NucleationProcess) {
    if (selections.nucleation == SelectedProcesses::MAMNucleation) {
      process = new MAMNucleationProcess();
#if HAERO_FORTRAN
    } else if (selections.nucleation == SelectedProcesses::MAMFNucleation) {
      process = new MAMNucleationFProcess();
#endif
    } else if (selections.nucleation == SelectedProcesses::NoNucleation) {
      process = new NullAerosolProcess(type);
    }
  } else if (type == CalcsizeProcess) {
    if (selections.calcsize == SelectedProcesses::MAMCalcsize) {
      process = new MAMCalcsizeProcess();
#if HAERO_FORTRAN
    } else if (selections.calcsize == SelectedProcesses::MAMFCalcsize) {
      process = new MAMCalcsizeFProcess();
#endif
    } else if (selections.calcsize == SelectedProcesses::NoCalcsize) {
      process = new NullAerosolProcess(type);
    }
  } else if (type == RenameProcess) {
    if (selections.rename == SelectedProcesses::MAMRename) {
      process = new MAMRenameProcess();
#if HAERO_FORTRAN
    } else if (selections.rename == SelectedProcesses::MAMFRename) {
      process = new MAMRenameFProcess();
#endif
    } else if (selections.rename == SelectedProcesses::NoRename) {
      process = new NullAerosolProcess(type);
    }
  } else if (type == GasAerosolExchangeProcess) {
    if (selections.gasaerosolexchange == SelectedProcesses::MAMGasAerosolExchange) {
      process = new MAMGasAerosolExchangeProcess();
#if HAERO_FORTRAN
    } else if (selections.gasaerosolexchange == SelectedProcesses::MAMFGasAerosolExchange) {
      process = new MAMGasAerosolExchangeFProcess();
#endif
    } else if (selections.gasaerosolexchange == SelectedProcesses::NoGasAerosolExchange) {
      process = new NullAerosolProcess(type);
    }
  } else if (type == ResuspensionProcess) {
    if (selections.resuspension == SelectedProcesses::NoResuspension) {
      process = new NullAerosolProcess(type);
    }
  }
  EKAT_REQUIRE_MSG(process != nullptr, "No aerosol process was selected!");
  return process;
}

}  // namespace haero
