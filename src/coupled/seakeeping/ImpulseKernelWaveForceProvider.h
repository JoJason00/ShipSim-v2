#pragma once

#include "IWaveForceProvider.h"
#include "../hydro/CoupledExcitationKernelRepo.h"
#include "../../config/CoupledConfig.h"
#include "../../config/SeakeepingConfig.h"
#include "../../seakeeping/TDGFExcitationKernelIO.h"

#include <vector>

class ImpulseKernelWaveForceProvider : public IWaveForceProvider
{
public:
    ImpulseKernelWaveForceProvider(const SeakeepingConfig& skCfg,
                                   const CoupledExcitationKernelRepo& excRepo,
                                   const HydroRefreshConfig& refreshCfg);

    void onWindowBegin(const CoupledSlowState3DOF& slow0,
                       const CoupledEncounterState& enc0,
                       double dtFast,
                       const CoupledFastWindowInfo& win = {}) override;

    Eigen::VectorXd evaluate(const CoupledSlowState3DOF& slow,
                             double t,
                             const CoupledEncounterState& enc,
                             Eigen::RowVectorXd* force6 = nullptr) override;

    void setEtaHistory(const std::vector<double>& hist) override;
    std::vector<double> etaHistory() const override { return etaHist_; }
    int memoryLength() const override { return memoryLen_; }

private:
    bool needRefresh(double Fn, double betaRel, double omegaInc, double dtFast) const;

private:
    SeakeepingConfig                  skCfg_;
    const CoupledExcitationKernelRepo& excRepo_;
    HydroRefreshConfig                refreshCfg_;

    bool                              cacheValid_ = false;
    double                            FnCached_    = 0.0;
    double                            betaCached_  = 0.0;
    double                            omegaCached_ = 0.0;
    double                            dtCached_    = 0.0;
    TDGFExcitationKernelData          excCached_{};

    std::vector<double>               etaHist_;
    int                               memoryLen_ = 0;
};
