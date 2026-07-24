#pragma once

#include "IWaveForceProvider.h"
#include "../../config/CaseConfig.h"
#include "../../config/SeakeepingConfig.h"
#include "../../seakeeping/Element.h"

#include <memory>

class DirectPressureFKWaveForceProvider : public IWaveForceProvider
{
public:
    DirectPressureFKWaveForceProvider(const ShipConfig& ship,
                                      const SeakeepingConfig& skCfg,
                                      std::shared_ptr<Element> element);

    void onWindowBegin(const CoupledSlowState3DOF& slow0,
                       const CoupledEncounterState& enc0,
                       double dtFast,
                       const CoupledFastWindowInfo& win = {}) override;

    Eigen::VectorXd evaluate(const CoupledSlowState3DOF& slow,
                             double t,
                             const CoupledEncounterState& enc,
                             Eigen::RowVectorXd* force6 = nullptr) override;

    int memoryLength() const override { return 0; }

private:
    ShipConfig               ship_;
    SeakeepingConfig         skCfg_;
    std::shared_ptr<Element> element_;
};
