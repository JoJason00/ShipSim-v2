#pragma once

#include <memory>
#include <variant>
#include <wave/WaveBase.h>
#include "SpectrumParams.h"   // SpectrumType + SpectrumParams

enum class WaveType
{
    regularwave,
    irregularwave,
    crosswave
};

struct RegularWaveConfig
{
    double direction;
    double H;
    double w;
    double phase0 = 0.0;   // initial phase ε [rad]; default 0

    RegularWaveConfig(double dir, double H, double w, double phase0_ = 0.0)
        : direction(dir), H(H), w(w), phase0(phase0_) {
    }
};

struct IrregularWaveConfig
{
    double direction = 0.0;   // absolute propagation direction [rad]
    SpectrumParams spectrum;
};

//struct CrossWaveConfig
//{
//    std::shared_ptr<struct WaveConfig> wave1;
//    std::shared_ptr<struct WaveConfig> wave2;
//    double angle_between;
//};

struct CrossWaveConfig {
    std::shared_ptr<WaveBase> wave1;
    std::shared_ptr<WaveBase> wave2;

    // Ship start position in the Earth frame [m]. For a crossing sea the two
    // sub-waves have different (k, θ); the same (x0, y0) shifts each sub-wave's
    // phase by a different amount −k_s(x0 cosθ_s + y0 sinθ_s), so it changes
    // their *relative* phase at the ship and hence the combined excitation.
    // Default (0,0) keeps every existing case bit-identical.
    double startX = 0.0;
    double startY = 0.0;
};

struct WaveConfig
{
    WaveType type;

    std::variant<
        RegularWaveConfig,
        IrregularWaveConfig,
        CrossWaveConfig
    > config;
};



