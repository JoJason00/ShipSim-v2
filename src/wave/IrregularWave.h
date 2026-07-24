#pragma once

#include "WaveBase.h"
#include "RegularWave.h"
#include "../config/WaveConfig.h"
#include <memory>

// Long-crested irregular wave at a single direction. The FK impulse response
// depends only on (direction, U) — not on the spectrum — so the kernel build
// is delegated to an internal RegularWave at the same direction. The spectrum
// only enters the simulation through the WaveField / incidentEtaAtShip sum
// built from spectrum() (see LinearCumminsTDGF). Eta() is therefore unused
// here and returns 0.
class IrregularWave : public WaveBase
{
public:
    IrregularWave(const IrregularWaveConfig& irregularwave);

    double Eta(double t) const override;             // unused (returns 0)
    void Exciting(double tn, FKphi& fkphi) override;  // delegate -> impulse reg
    void loadData(const fkpData& Data) override;      // forward U to delegate

    double getAmp() override;                         // significant amplitude
    double getFreq() override;                        // representative ω_inc
    double direction() override;

    const SpectrumParams& spectrum() const { return config.spectrum; }

private:
    const IrregularWaveConfig config;
    std::unique_ptr<RegularWave> impulseReg_;         // kernel delegate
    double repFreq_ = 0.0;                            // spectral peak ω
    double sigAmp_ = 0.0;                             // 2·√m0 (display scale)
};
