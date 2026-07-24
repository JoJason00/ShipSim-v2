#include "CoupledZigzagScenario.h"
#include <cmath>

namespace
{
    constexpr double PI_local = 3.14159265358979323846;
}

double CoupledZigzagScenario::deg2rad(double deg)
{
    return deg * PI_local / 180.0;
}

double CoupledZigzagScenario::wrapToPi(double x)
{
    while (x <= -PI_local) x += 2.0 * PI_local;
    while (x > PI_local) x -= 2.0 * PI_local;
    return x;
}

CoupledZigzagScenario::CoupledZigzagScenario(const ZigzagCfg& cfg, double approachTime)
    : cfg_(cfg), approachTime_(approachTime)
{
    zigAngle_ =
        (std::abs(cfg_.rudder_angle) > 1.0e-12)
        ? std::abs(cfg_.rudder_angle)
        : deg2rad(10.0);

    switchTarget_ =
        (cfg_.ZigzagCaseCicle > 0)
        ? static_cast<int>(cfg_.ZigzagCaseCicle)
        : 4;
}

void CoupledZigzagScenario::reset(const CoupledSlowState3DOF& s0)
{
    started_ = false;
    finished_ = false;
    psi0_ = s0.psi;
    currentCmd_ = +zigAngle_;
    switchCount_ = 0;
}

CoupledManeuverCommand CoupledZigzagScenario::evaluate(
    const CoupledSlowState3DOF& s,
    double t)
{
    CoupledManeuverCommand cmd{};

    if (finished_)
    {
        cmd.deltaCmd = 0.0;
        cmd.started = started_;
        cmd.finished = true;
        return cmd;
    }

    if (!started_)
    {
        if (t >= approachTime_)
        {
            started_ = true;
            psi0_ = s.psi;
            currentCmd_ = +zigAngle_;
            switchCount_ = 0;
        }
    }
    else
    {
        const double psiRel = wrapToPi(s.psi - psi0_);

        if (currentCmd_ > 0.0 && psiRel >= +zigAngle_)
        {
            currentCmd_ = -zigAngle_;
            ++switchCount_;
        }
        else if (currentCmd_ < 0.0 && psiRel <= -zigAngle_)
        {
            currentCmd_ = +zigAngle_;
            ++switchCount_;
        }

        if (switchCount_ >= switchTarget_)
            finished_ = true;
    }

    cmd.deltaCmd = started_ ? currentCmd_ : 0.0;
    cmd.started = started_;
    cmd.finished = finished_;
    return cmd;
}
