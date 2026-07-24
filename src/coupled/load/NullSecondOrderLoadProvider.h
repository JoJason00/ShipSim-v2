#pragma once

#include "ISecondOrderLoadProvider.h"

class NullSecondOrderLoadProvider : public ISecondOrderLoadProvider
{
public:
    CoupledExternalLoads3DOF sample(
        const CoupledSlowState3DOF&,
        const CoupledEncounterState&) const override
    {
        return CoupledExternalLoads3DOF{};
    }
};
