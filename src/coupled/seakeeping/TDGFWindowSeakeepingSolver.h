#pragma once

#include "IWindowSeakeepingSolver.h"
#include "IWaveForceProvider.h"
#include "../hydro/CoupledRadiationKernelRepo.h"
#include "../../config/CaseConfig.h"
#include "../../config/SeakeepingConfig.h"
#include "../../config/CoupledConfig.h"
#include "../../seakeeping/CumminsTimeStepper.h"
#include "../../seakeeping/RadiationKernelCache.h"
#include "../../seakeeping/SeakeepingDOF.h"
#include "../../seakeeping/Element.h"

#include <Eigen/Dense>
#include <memory>
#include <string>

// TDGFWindowSeakeepingSolver — coupled fast-window solver using exactly the same
// Newmark-β + trapezoid Cummins time-stepper as pure seakeeping.
//
// Speed-aware kernel switching:
//   The radiation kernel A∞, B, C', K(t) depends on Fn. The stepper is fully rebuilt
//   only when |Fn - Fn_cached| / max(|Fn_cached|, 1e-3) exceeds refreshCfg.speedTolRatio.
//   Otherwise the cached stepper is reused across windows (the state q, v, a, vHist
//   carries over via setInitialState).
//
// The excitation source is abstracted behind IWaveForceProvider so the same solver
// supports both impulse-kernel (Qlag) and direct-pressure-FK wave force modes.

class TDGFWindowSeakeepingSolver : public IWindowSeakeepingSolver
{
public:
    TDGFWindowSeakeepingSolver(const ShipConfig& ship,
                               const SeakeepingConfig& skCfg,
                               const std::string& casePath,
                               const CoupledRadiationKernelRepo& radRepo,
                               const HydroRefreshConfig& refreshCfg,
                               IWaveForceProvider& forceProvider,
                               std::shared_ptr<Element> elementForHydrostatics,
                               bool enableRadiation = true);

    CoupledWindowResult solveWindow( const CoupledWindowRequest& req,
                                    const CoupledWaveEnvironment& env) override;

private:
    bool needRadRefresh(double Fn, double dtFast) const;

    void postprocessRadiationKernel(
        double Fn,
        RadiationKernelData& kernel) const;

    void refreshRadiation(double Fn, double dtFast);

    static CoupledSlowState3DOF interpSlow(const CoupledSlowState3DOF& a,
                                           const CoupledSlowState3DOF& b,
                                           double alpha);

private:
    ShipConfig                          ship_;
    SeakeepingConfig                    skCfg_;
    std::string                         casePath_;
    const CoupledRadiationKernelRepo&   radRepo_;
    HydroRefreshConfig                  refreshCfg_;
    IWaveForceProvider&                 forceProvider_;
    std::shared_ptr<Element>            element_;
    bool                                enableRadiation_ = true;

    // Hydrostatics — built once from the panel mesh (same path as pure seakeeping).
    HydrostaticsData                    hs_{};
    Eigen::MatrixXd                     Mphys_, Bphys_, Cphys_;
    RollViscDamping                     rollVisc_;
    int                                 rollIndex_ = -1;

    // Cached radiation kernel + stepper. Rebuilt only when |ΔFn| exceeds threshold.
    bool                                radValid_ = false;
    double                              FnRadCached_ = 0.0;
    double                              dtCached_    = 0.0;
    RadiationKernelData                 radCached_{};
    std::unique_ptr<cummins::CumminsTimeStepper> stepper_;
};
