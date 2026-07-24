#pragma once
#include "ICoupledManeuverScenario.h"
#include "../../config/MmgConfig.h"

class CoupledTurningScenario : public ICoupledManeuverScenario
{
public:
    // tMaxAbs：有量纲时间上限 [s]，>0 时启用"按时间停"判据（由 t_star_max 换算而来）；
    //          ≤0 时仅按圈数停（原行为）。两种判据并存，谁先满足谁停。
    CoupledTurningScenario(const TuringCfg& cfg, bool port,
                           double approachTime = 0.0,
                           double tMaxAbs = 0.0);

    void reset(const CoupledSlowState3DOF& s0) override;

    CoupledManeuverCommand evaluate(
        const CoupledSlowState3DOF& s,
        double t) override;

    std::string tag() const override { return tag_; }

private:
    static double deg2rad(double deg);
    static double wrapToPi(double x);

private:
    TuringCfg cfg_{};
    bool port_ = true;
    double approachTime_ = 0.0;
    double tMaxAbs_ = 0.0;

    bool started_ = false;
    bool finished_ = false;

    double deltaTarget_ = 0.0;
    double turnAngle_ = 0.0;
    double psiPrev_ = 0.0;
    double psiAcc_ = 0.0;

    std::string tag_;
};