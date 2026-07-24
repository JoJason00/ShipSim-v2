#include "DirectPressureFKWaveForceProvider.h"
#include "../../seakeeping/DirectPressureFK.h"
#include "../../wave/Encounter.h"
#include "../../const/Const.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
    double wrapPiLocal(double x)
    {
        while (x <= -PI) x += 2.0 * PI;
        while (x >   PI) x -= 2.0 * PI;
        return x;
    }
}

DirectPressureFKWaveForceProvider::DirectPressureFKWaveForceProvider(
    const ShipConfig& ship,
    const SeakeepingConfig& skCfg,
    std::shared_ptr<Element> element)
    : ship_(ship), skCfg_(skCfg), element_(std::move(element))
{
    if (!element_)
        throw std::runtime_error("DirectPressureFKWaveForceProvider: element is null.");
}

void DirectPressureFKWaveForceProvider::onWindowBegin(
    const CoupledSlowState3DOF& /*slow0*/,
    const CoupledEncounterState& /*enc0*/,
    double /*dtFast*/,
    const CoupledFastWindowInfo& /*win*/)
{
    // Stateless — nothing to refresh.
}

Eigen::VectorXd DirectPressureFKWaveForceProvider::evaluate(
    const CoupledSlowState3DOF& slow,
    double t,
    const CoupledEncounterState& enc,
    Eigen::RowVectorXd* force6)
{
    // C3a: linear FK force = Σ over incident components (the integrand is
    // linear in wave amplitude). Regular wave = 1 component == the scalar
    // path -> bit-identical. Irregular/cross = many components summed.
    Eigen::RowVectorXd F6;

    // Body-frame FK integration. The trajectory enters only through Φ0, the
    // wave-field phase at the ship (body) origin evaluated with the real-time
    // earth position (slow.xe, slow.ye). DirectPressureFK::computeForce6 forms
    // the hull-point phase as Φ0 − k(x_b cosβ + y_b sinβ) — the NEGATED
    // standard phase (the pure-seakeeping caller passes Φ0 = +ω_e t + ε, see
    // LinearCumminsTDGF). So Φ0 here must be the negation of the true phase
    // k(xe cosθ + ye sinθ) − ω t + ε; then dΦ0/dt = ω − kU cosβ_rel = ω_e,
    // consistent with calcEncounterInfo. (The old form −k(...) − ωt swapped
    // the head/following encounter frequencies during a turn.)
    auto phi0 = [&](double a_omega, double thetaAbs, double eps) -> double
    {
        const double k = a_omega * a_omega / G;
        return -k * (slow.xe * std::cos(thetaAbs) + slow.ye * std::sin(thetaAbs))
               + a_omega * t - eps;
    };

    if (!enc.comps.empty())
    {
        for (const auto& w : enc.comps)
        {
            DirectPressureFKContext ctx;
            ctx.amp       = w.a;
            ctx.omega     = w.omega;
            ctx.beta      = wrapPiLocal(w.theta - slow.psi);
            ctx.phiOrigin = phi0(w.omega, w.theta, w.eps);

            Eigen::RowVectorXd f = DirectPressureFK::computeForce6(*element_, ctx, t);
            if (F6.size() == 0) F6 = f;
            else                F6 += f;
        }
    }
    else
    {
        DirectPressureFKContext ctx;
        ctx.amp       = enc.waveAmp;
        ctx.omega     = enc.omega;
        ctx.beta      = enc.betaRel;
        ctx.phiOrigin = phi0(enc.omega, enc.betaAbs, 0.0);
        F6 = DirectPressureFK::computeForce6(*element_, ctx, t);
    }

    const int D = skCfg_.DOF;
    Eigen::VectorXd F = Eigen::VectorXd::Zero(D);
    for (int i = 0; i < D; ++i)
    {
        const int mode = skCfg_.modes[i];
        if (mode >= 0 && mode < 6)
            F(i) = F6(mode);
    }

    if (force6)
        *force6 = F6;

    return F;
}
