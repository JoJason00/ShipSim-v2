#pragma once

#include "WaveComponent.h"
#include "../config/SpectrumParams.h"
#include <vector>

// Spectral density S(ω) [m² s] for PM / ITTC / JONSWAP / Ochi-Hubble, plus an
// equal-energy (equal-area) discretiser. Equal-energy banding gives every
// component the same energy m0/N (uniform statistical weight, fast
// convergence); the representative frequency is drawn at random inside each
// band to break the deterministic 2π/Δω recurrence over long runs. The RNG is
// seeded for reproducibility.

namespace wave_spectrum
{
    double density(const SpectrumParams& p, double omega);

    double m0(const SpectrumParams& p, int nQuad = 4000);

    // theta = absolute wave propagation direction [rad], stamped on every
    // returned component.
    std::vector<WaveComponent> discretize(const SpectrumParams& p, double theta);
}
