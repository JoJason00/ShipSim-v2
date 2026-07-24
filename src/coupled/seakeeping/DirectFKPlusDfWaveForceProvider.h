#pragma once

// DirectFKPlusDfWaveForceProvider — hybrid first-order wave load for coupling:
//   F = F_FK^direct + F_D^impulse
//
//   F_FK: panel Gauss integration of incident pressure (DirectPressureFK).
//   F_D:  convolution of fkImpulse_*.csv diffraction (df) kernels with η_ref(t),
//         same grouping / trajectory / WaveField phase as FKImpulseWaveForceProvider.
//
// JSON:  "waveForceMode": "directFKdf"

#include "IWaveForceProvider.h"
#include "DirectPressureFKWaveForceProvider.h"
#include "FKImpulseWaveForceProvider.h"

#include <memory>

class DirectFKPlusDfWaveForceProvider : public IWaveForceProvider
{
public:
    DirectFKPlusDfWaveForceProvider(const ShipConfig& ship,
                                    const SeakeepingConfig& skCfg,
                                    const std::string& casePath,
                                    const std::string& shipName,
                                    const HydroRefreshConfig& refreshCfg,
                                    std::shared_ptr<Element> element);

    void onWindowBegin(const CoupledSlowState3DOF& slow0,
                       const CoupledEncounterState& enc0,
                       double dtFast,
                       const CoupledFastWindowInfo& win = {}) override;

    Eigen::VectorXd evaluate(const CoupledSlowState3DOF& slow,
                             double t,
                             const CoupledEncounterState& enc,
                             Eigen::RowVectorXd* force6 = nullptr) override;

    void setEtaHistory(const std::vector<double>& hist) override;
    std::vector<double> etaHistory() const override;
    int memoryLength() const override { return dfImpulse_.memoryLength(); }

private:
    DirectPressureFKWaveForceProvider fkDirect_;
    FKImpulseWaveForceProvider      dfImpulse_;
    std::string                     casePath_;  // 仅供临时诊断 dump 用
};
