#include "IrregularWave.h"
#include "WaveSpectrum.h"

#include <cmath>

namespace
{
    // Spectral-peak (representative) incident frequency: argmax S over range.
    double peakFreq(const SpectrumParams& p)
    {
        const int M = std::max(2000, 20 * p.nComponents);
        const double dw = (p.omegaHi - p.omegaLo) / M;
        double wBest = p.omegaLo, sBest = -1.0;
        for (int i = 0; i <= M; ++i)
        {
            const double w = p.omegaLo + i * dw;
            const double s = wave_spectrum::density(p, w);
            if (s > sBest) { sBest = s; wBest = w; }
        }
        return wBest;
    }
}

IrregularWave::IrregularWave(const IrregularWaveConfig& irregularwave)
    : config(irregularwave)
{
    repFreq_ = peakFreq(config.spectrum);
    // Hs ≈ 4√m0  ->  characteristic amplitude Hs/2 = 2√m0 (display scaling).
    sigAmp_ = 2.0 * std::sqrt(std::max(wave_spectrum::m0(config.spectrum), 0.0));

    // Impulse-kernel delegate: same direction; H/ω irrelevant to Exciting
    // (the FK impulse response depends only on direction and loaded U).
    impulseReg_ = std::make_unique<RegularWave>(
        RegularWaveConfig(config.direction, 1.0, repFreq_));
}

double IrregularWave::Eta(double /*t*/) const
{
    return 0.0; // unused: incidentEtaAtShip sums the WaveField components
}

void IrregularWave::Exciting(double tn, FKphi& fkphi)
{
    impulseReg_->Exciting(tn, fkphi);
}

void IrregularWave::loadData(const fkpData& Data)
{
    impulseReg_->loadData(Data);
}

double IrregularWave::getAmp()  { return sigAmp_; }
double IrregularWave::getFreq() { return repFreq_; }
double IrregularWave::direction() { return config.direction; }
