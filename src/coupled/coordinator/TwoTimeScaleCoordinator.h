#pragma once

#include "../../config/CaseConfig.h"
#include "../common/CoupledTypes.h"
#include "../io/CoupledWriter.h"
#include "../environment/WaveEnvironment.h"
#include "../seakeeping/IWindowSeakeepingSolver.h"
#include "../maneuver/ILowFreqManeuverSolver.h"
#include "../maneuver/ICoupledManeuverScenario.h"
#include "../load/ISecondOrderLoadProvider.h"

class TwoTimeScaleCoordinator
{
public:
    TwoTimeScaleCoordinator(
        const CaseConfig& cfg,
        const std::string& casePath,
        IWindowSeakeepingSolver& fastSolver,
        ICoupledLowFreqManeuverSolver& slowSolver,
        ISecondOrderLoadProvider& secondOrderProvider,
        CoupledWaveEnvironment& waveEnv);

    void run();

private:
    CoupledWindowRequest buildWindowRequest(
        const CoupledSlowState3DOF& s0,
        const CoupledSlowState3DOF& s1pred) const;

    void runScenario(
        ICoupledManeuverScenario& scenario,
        double Fn,
        int iWave);

private:
    CaseConfig cfg_;
    std::string casePath_;

    IWindowSeakeepingSolver& fastSolver_;
    ICoupledLowFreqManeuverSolver& slowSolver_;
    ISecondOrderLoadProvider& secondOrderProvider_;
    CoupledWaveEnvironment& waveEnv_;

    CoupledWriter writer_;
};