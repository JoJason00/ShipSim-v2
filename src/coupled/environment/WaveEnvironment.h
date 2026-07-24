#pragma once

#include "../common/CoupledTypes.h"
#include "../../wave/WaveBase.h"
#include <memory>
#include <vector>

class CoupledWaveEnvironment
{
public:
    explicit CoupledWaveEnvironment(const std::vector<std::shared_ptr<WaveBase>>& waves);

    /// Which entry in `waves` is used for encounter / drift (default 0). Coupling sets this to `iWave` per run.
    void setActiveWaveIndex(int i);

    /// Soft-start ramp duration (seconds) for the incident wave amplitude.
    /// <= 0 disables the ramp. See CoupledTimeConfig::waveRampTimeS.
    void setWaveRampTime(double seconds);

    CoupledEncounterState evaluateEncounter(const CoupledSlowState3DOF& s, double t) const;
    double evaluateWaveElevationAtCG(const CoupledSlowState3DOF& s, double t) const;

    /// Incident wave scalars for the active `iWave` (deep-water k = ω²/g).
    struct ActiveWaveSpec
    {
        double omegaIncident = 0.0;
        double zetaA = 0.0;
        double k = 0.0;
        double directionRad = 0.0;  // absolute propagation direction [rad]
    };
    ActiveWaveSpec activeWaveSpec() const;

private:
    static double wrapAngle(double x);

    // Time-independent incident decomposition (a, ω, θ_abs, ε) of the active
    // wave: regular = 1 component; irregular = spectrum discretisation; cross
    // = sub-waves merged. Built once per active wave (spectrum discretisation
    // is expensive) and reused every fast step; only the ramp scales it.
    const std::vector<WaveComponent>& baseComponents() const;

private:
    std::vector<std::shared_ptr<WaveBase>> waves_;
    int activeWaveIndex_{ 0 };
    double waveRampTimeS_{ 0.0 };

    mutable int baseCompsWave_{ -1 };
    mutable std::vector<WaveComponent> baseComps_;
};
