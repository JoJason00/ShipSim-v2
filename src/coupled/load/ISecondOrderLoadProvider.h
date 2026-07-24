#pragma once

#include "../common/CoupledTypes.h"
#include "../environment/WaveEnvironment.h"

class ISecondOrderLoadProvider
{
public:
    virtual ~ISecondOrderLoadProvider() = default;

    virtual CoupledExternalLoads3DOF sample(
        const CoupledSlowState3DOF& s,
        const CoupledEncounterState& enc) const = 0;
};
