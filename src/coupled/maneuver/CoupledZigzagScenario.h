#pragma once
#include "ICoupledManeuverScenario.h"
#include "../../config/MmgConfig.h"

class CoupledZigzagScenario : public ICoupledManeuverScenario
{
public:
    CoupledZigzagScenario(const ZigzagCfg& cfg, double approachTime = 0.0);

    void reset(const CoupledSlowState3DOF& s0) override;

    CoupledManeuverCommand evaluate(
        const CoupledSlowState3DOF& s,
        double t) override;

    std::string tag() const override { return tag_; }

private:
    static double deg2rad(double deg);
    static double wrapToPi(double x);

private:
    ZigzagCfg cfg_{};
    double approachTime_ = 0.0;

    bool started_ = false;
    bool finished_ = false;

    double psi0_ = 0.0;
    double currentCmd_ = 0.0;
    double zigAngle_ = 0.0;
    int switchCount_ = 0;
    int switchTarget_ = 4;

    std::string tag_ = "zigzag";
};
