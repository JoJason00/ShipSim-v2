#pragma once
#include "IWindowSeakeepingSolver.h"
#include <Eigen/Dense>

class NullWindowSeakeepingSolver : public IWindowSeakeepingSolver
{
public:
    CoupledWindowResult solveWindow(
        const CoupledWindowRequest& req,
        const CoupledWaveEnvironment& env) override
    {
        (void)req;
        (void)env;

        CoupledWindowResult out;
        out.loads = CoupledWaveLoadsWindow{};
        out.etaMean = Eigen::VectorXd::Zero(3);
        out.etaRms = Eigen::VectorXd::Zero(3);
        out.fastLast.initialized = false;
        return out;
    }
};