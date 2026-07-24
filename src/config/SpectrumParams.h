#pragma once

#include <cstdint>

enum class SpectrumType
{
    JONSWAP,
    PM,
    ITTC,
    OH
};

// All possible spectrum parameters. The case file fills only the fields the
// chosen `type` needs (see WaveSpectrum.cpp for the exact formulae).
struct SpectrumParams
{
    SpectrumType type = SpectrumType::ITTC;

    // PM, eq.(4): S = a g² / ω⁵ · exp(-β (g/(U ω))⁴), a=0.0081, β=0.74.
    double Uwind = 0.0;           // wind speed at 19.5 m [m/s]

    // ITTC two-parameter, eq.(5): S = 173 Hs²/(T1⁴ ω⁵) · exp(-691/(T1⁴ ω⁴)).
    double Hs = 0.0;              // significant wave height [m]
    double T1 = 0.0;              // mean period [s]

    // JONSWAP, eq.(6): parameterised by Hs, Tp, gamma (ωm = 2π/Tp).
    double Tp = 0.0;              // peak period [s]
    double gamma = 3.3;           // peak enhancement factor

    // Ochi-Hubble bimodal, eq.(7): six parameters.
    double Hs1 = 0.0, Hs2 = 0.0;  // significant heights of the two systems [m]
    double wm1 = 0.0, wm2 = 0.0;  // peak frequencies [rad/s]
    double lam1 = 0.0, lam2 = 0.0;// shape parameters
    // 0 = both systems summed; 1 = system-1 (low-freq) only;
    // 2 = system-2 (high-freq) only. Used to feed the two OH peaks into two
    // separate directions for a "split-OH" cross sea.
    int ohMode = 0;

    // Equal-energy discretisation controls.
    double omegaLo = 0.2;
    double omegaHi = 3.0;
    int    nComponents = 100;
    std::uint64_t seed = 1;
    bool   randomFreqInBand = true;
};
