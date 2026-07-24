#pragma once

#include "../../config/CaseConfig.h"
#include "../../config/MmgConfig.h"
#include "ILowFreqManeuverSolver.h"

#include <memory>
#include <string>

/// Builds the coupled slow-time manoeuvring solver from `Manoeuvring.lowFreqManeuverModel`.
std::unique_ptr<ICoupledLowFreqManeuverSolver> makeLowFreqManeuverSolver(
    const ShipConfig& ship,
    const MmgConfig& mmg);

/// Human-readable tag for logs (stable identifier).
std::string lowFreqManeuverModelTag(const MmgConfig& mmg);

/// Subfolder under the case directory for coupled CSV output, e.g. "" -> `coupled_turning`,
/// `"_sjtu"` -> `coupled_turning_sjtu` (paired with `lowFreqManeuverModel`).
std::string lowFreqManeuverOutputFolderSuffix(const MmgConfig& mmg);
