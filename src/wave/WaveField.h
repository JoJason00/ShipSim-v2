#pragma once

#include "WaveComponent.h"
#include "WaveSpectrum.h"
#include <vector>

// Total incident sea as a flat list of WaveComponent. Regular = 1 component;
// long-crested irregular = N components at one direction; cross / directional
// = concatenation of several direction systems. The kernel/assembler layer
// only ever sees components(), so the three sea types share one code path.
//
// Phase convention matches the existing single-RegularWave path
// (CoupledWaveEnvironment): phase = -k(x cosθ + y sinθ) - ω t + ε,
// k = ω² / g (deep water). Keep this identical to how the impulse kernels
// were generated.
class WaveField
{
public:
    WaveField() = default;
    explicit WaveField(std::vector<WaveComponent> comps)
        : comps_(std::move(comps)) {}

    static WaveField makeRegular(double amp, double dirRad, double omega,
                                 double eps = 0.0);
    static WaveField makeIrregular(const SpectrumParams& sp, double dirRad);
    static WaveField makeCross(const WaveField& a, const WaveField& b);

    const std::vector<WaveComponent>& components() const { return comps_; }
    bool empty() const { return comps_.empty(); }

    // Incident free-surface elevation at point (x, y) and time t [m].
    double elevationAt(double x, double y, double t) const;

private:
    std::vector<WaveComponent> comps_;
};
