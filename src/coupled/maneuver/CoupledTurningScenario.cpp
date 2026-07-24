#include "CoupledTurningScenario.h"
#include <cmath>

namespace
{
    constexpr double PI_local = 3.14159265358979323846;
}

double CoupledTurningScenario::deg2rad(double deg)
{
    return deg * PI_local / 180.0;
}

double CoupledTurningScenario::wrapToPi(double x)
{
    while (x <= -PI_local) x += 2.0 * PI_local;
    while (x > PI_local) x -= 2.0 * PI_local;
    return x;
}

CoupledTurningScenario::CoupledTurningScenario(const TuringCfg& cfg, bool port, double approachTime, double tMaxAbs)
    : cfg_(cfg), port_(port), approachTime_(approachTime), tMaxAbs_(tMaxAbs)
{
    const double deltaMag =
        (std::abs(cfg_.rudder_angle) > 1.0e-12)
        ? std::abs(cfg_.rudder_angle)
        : deg2rad(35.0);

    // δ > 0 in this MMG implementation drives a starboard (right) turn:
    // NR = -(xR + aH·xH)·FN·cos(δ) with xR, xH < 0 ⇒ NR > 0 ⇒ r > 0 ⇒ ψ̇ > 0
    // (z-down, right-hand rule). So port (left) turn must command -|δ|.
    deltaTarget_ = port_ ? -deltaMag : deltaMag;

    const double turningCircle =
        (cfg_.TurningCaseCicle > 0.0)
        ? cfg_.TurningCaseCicle
        : 540.0 / 360.0;

    turnAngle_ = turningCircle * 2.0 * PI_local;
    tag_ = port_ ? "port" : "starboard";
}

void CoupledTurningScenario::reset(const CoupledSlowState3DOF& s0)
{
    started_ = false;
    finished_ = false;
    psiPrev_ = s0.psi;
    psiAcc_ = 0.0;
}

CoupledManeuverCommand CoupledTurningScenario::evaluate(
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
            psiPrev_ = s.psi;
            psiAcc_ = 0.0;
        }
    }
    else
    {
        const double dpsi = wrapToPi(s.psi - psiPrev_);
        psiAcc_ += dpsi;
        psiPrev_ = s.psi;

        // 圈数判据（仅在未启用时间上限时生效，由 stop_by 选定）
        if (tMaxAbs_ <= 0.0 && std::abs(psiAcc_) >= turnAngle_)
            finished_ = true;
    }

    // 时间判据：按无量纲时间换算得到的有量纲上限（stop_by="t_star" 时 tMaxAbs_>0）
    if (tMaxAbs_ > 0.0 && t >= tMaxAbs_)
        finished_ = true;

    cmd.deltaCmd = started_ ? deltaTarget_ : 0.0;
    cmd.started = started_;
    cmd.finished = finished_;
    return cmd;
}