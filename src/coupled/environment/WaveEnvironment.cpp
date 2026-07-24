#include "WaveEnvironment.h"
#include "../../wave/RegularWave.h"
#include "../../wave/IrregularWave.h"
#include "../../wave/CrossWave.h"
#include "../../wave/WaveSpectrum.h"
#include "../../const/Const.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

CoupledWaveEnvironment::CoupledWaveEnvironment(const std::vector<std::shared_ptr<WaveBase>>& waves)
    : waves_(waves)
{
}

void CoupledWaveEnvironment::setActiveWaveIndex(int i)
{
    if (waves_.empty())
    {
        activeWaveIndex_ = 0;
        return;
    }
    const int n = static_cast<int>(waves_.size());
    activeWaveIndex_ = std::clamp(i, 0, n - 1);
}

void CoupledWaveEnvironment::setWaveRampTime(double seconds)
{
    waveRampTimeS_ = seconds;
}

double CoupledWaveEnvironment::wrapAngle(double x)
{
    while (x <= -PI) x += 2.0 * PI;
    while (x > PI) x -= 2.0 * PI;
    return x;
}

const std::vector<WaveComponent>&
CoupledWaveEnvironment::baseComponents() const
{
    const int wi = std::clamp(activeWaveIndex_, 0,
                              static_cast<int>(waves_.size()) - 1);
    if (baseCompsWave_ == wi && !baseComps_.empty())
        return baseComps_;

    std::function<std::vector<WaveComponent>(const std::shared_ptr<WaveBase>&)>
        build = [&](const std::shared_ptr<WaveBase>& w)
        -> std::vector<WaveComponent>
    {
        if (auto rw = std::dynamic_pointer_cast<RegularWave>(w))
            return { WaveComponent{ rw->getAmp(), rw->getFreq(),
                                    rw->direction(), rw->initialPhase() } };
        if (auto ir = std::dynamic_pointer_cast<IrregularWave>(w))
            return wave_spectrum::discretize(ir->spectrum(), ir->direction());
        if (auto cw = std::dynamic_pointer_cast<CrossWave>(w))
        {
            auto a = build(cw->wave1());
            auto b = build(cw->wave2());
            a.insert(a.end(), b.begin(), b.end());
            return a;
        }
        throw std::runtime_error(
            "CoupledWaveEnvironment: unsupported wave type "
            "(expected Regular/Irregular/Cross).");
    };

    baseComps_ = build(waves_[static_cast<std::size_t>(wi)]);
    baseCompsWave_ = wi;
    return baseComps_;
}

CoupledEncounterState CoupledWaveEnvironment::evaluateEncounter(const CoupledSlowState3DOF& s, double t) const
{
    CoupledEncounterState info{};

    if (waves_.empty())
        return info;

    const int wi = std::clamp(activeWaveIndex_, 0, static_cast<int>(waves_.size()) - 1);
    const auto& activeWave = waves_[static_cast<std::size_t>(wi)];

    double rampFactor = 1.0;
    if (waveRampTimeS_ > 0.0)
    {
        if (t <= 0.0)
            rampFactor = 0.0;
        else if (t < waveRampTimeS_)
            rampFactor = 0.5 * (1.0 - std::cos(PI * t / waveRampTimeS_));
    }

    // Component decomposition (regular = 1; irregular = spectrum; cross =
    // merged), ramp applied to amplitudes. Built once, cached.
    const auto& base = baseComponents();
    info.comps.clear();
    info.comps.reserve(base.size());
    double sumA2 = 0.0;
    for (const auto& c : base)
    {
        const double a = c.a * rampFactor;
        info.comps.push_back(WaveComponent{ a, c.omega, c.theta, c.eps });
        sumA2 += a * a;
    }

    // Representative scalars for any not-yet-migrated consumer / diagnostics.
    // Regular wave: betaAbs/omega/waveAmp reduce exactly to the legacy values
    // (1 component, ε=0) -> bit-identical.
    const double betaAbs = activeWave->direction();
    const double omega   = activeWave->getFreq();
    const double betaRel = wrapAngle(betaAbs - s.psi);

    info.betaAbs = betaAbs;
    info.betaRel = betaRel;
    info.omega   = omega;
    info.waveAmp = std::sqrt(sumA2);   // = a·ramp for a single regular comp
    try { info.kin = calcEncounterInfo(omega, s.U, betaRel); }
    catch (...) { info.kin = EncounterInfo{}; }

    const double k = omega * omega / G;
    // Standard propagation convention: η = a·cos(k(x cosβ + y sinβ) − ωt),
    // so the phase rate along the trajectory is dθ/dt = −(ω − kU cosβ_rel),
    // i.e. |ωe| matches calcEncounterInfo (head sea ω+kU, following |ω−kU|).
    // The old "−k" spatial sign propagated the wave toward −β and swapped the
    // head/following encounter frequencies during a turn.
    info.phaseAtCG = k * (s.xe * std::cos(betaAbs) + s.ye * std::sin(betaAbs)) - omega * t;

    return info;
}

double CoupledWaveEnvironment::evaluateWaveElevationAtCG(const CoupledSlowState3DOF& s, double t) const
{
    // Sum over EVERY component so irregular and crossing seas show their true
    // (multi-frequency / multi-direction) elevation, not a single representative
    // cosine. Same phase convention as WaveField::elevationAt — the incident η
    // that actually drives the wave force — so this diagnostic matches it.
    // A single regular component (ε=0) reduces to enc.waveAmp·cos(enc.phaseAtCG),
    // i.e. the previous behaviour bit-for-bit.
    const auto enc = evaluateEncounter(s, t);
    double eta = 0.0;
    for (const auto& c : enc.comps)
    {
        const double k = c.omega * c.omega / G;   // deep water
        const double phase =
            k * (s.xe * std::cos(c.theta) + s.ye * std::sin(c.theta))
            - c.omega * t + c.eps;
        eta += c.a * std::cos(phase);
    }
    return eta;
}

CoupledWaveEnvironment::ActiveWaveSpec CoupledWaveEnvironment::activeWaveSpec() const
{
    ActiveWaveSpec spec{};
    if (waves_.empty())
        return spec;

    const int wi = std::clamp(activeWaveIndex_, 0, static_cast<int>(waves_.size()) - 1);
    const auto& w = waves_[static_cast<std::size_t>(wi)];
    if (!w)
        return spec;

    // Use the polymorphic WaveBase interface so EVERY wave type yields a
    // representative spec — not just RegularWave. This feeds the output
    // folder name (Fn..nada..fai..) and the writer's kA/ω_e reference.
    //   RegularWave  : exact (getFreq/getAmp/direction == legacy values).
    //   IrregularWave: spectral peak ω, significant amplitude, its heading.
    //   CrossWave    : delegates to wave1 (peak ω / amp / heading).
    // Without this, Cross/Irregular collapsed to nada0.0_fai0 and every such
    // run overwrote the previous one in the same folder.
    spec.omegaIncident = w->getFreq();
    spec.zetaA         = w->getAmp();
    spec.k             = spec.omegaIncident * spec.omegaIncident / G;  // deep water
    spec.directionRad  = w->direction();
    return spec;
}
