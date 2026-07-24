#pragma once

#include "WaveBase.h"
#include "../config/WaveConfig.h"

class CrossWave : public WaveBase
{
public:
    CrossWave(const CrossWaveConfig& crosswave);

    double Eta(double t) const override;
    void Exciting(double tn, FKphi& fkphi) override;

    void loadData(const fkpData& Data) override;
    double getAmp() override;
    double getFreq() override;
    double direction() override;

    // Two crossing systems (each regular or irregular, each its own
    // direction). Consumed by the Step-4 multi-kernel assembler.
    const std::shared_ptr<WaveBase>& wave1() const { return config.wave1; }
    const std::shared_ptr<WaveBase>& wave2() const { return config.wave2; }

    // Ship start position [m] in the Earth frame (see CrossWaveConfig).
    double startX() const { return config.startX; }
    double startY() const { return config.startY; }

private:
    const CrossWaveConfig config;
};
