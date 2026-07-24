#pragma once
#include "../common/CoupledTypes.h"

class ICoupledLowFreqManeuverSolver
{
public:
    virtual ~ICoupledLowFreqManeuverSolver() = default;

    virtual CoupledSlowState3DOF step(
        const CoupledSlowState3DOF& s,
        const CoupledExternalLoads3DOF& ext,
        double deltaCmd,
        double dt) = 0;

    virtual double referenceSpeedFromFn(double Fn) const = 0;
};