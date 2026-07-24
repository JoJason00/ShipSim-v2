#pragma once
#include "../common/CoupledTypes.h"
#include <string>

struct CoupledManeuverCommand
{
    double deltaCmd = 0.0;
    bool started = false;
    bool finished = false;
};

class ICoupledManeuverScenario
{
public:
    virtual ~ICoupledManeuverScenario() = default;

    virtual void reset(const CoupledSlowState3DOF& s0) = 0;

    virtual CoupledManeuverCommand evaluate(
        const CoupledSlowState3DOF& s,
        double t) = 0;

    virtual std::string tag() const = 0;
};