#include "CoupledSjtuMmg3DOFCore.h"
#include "../../const/Const.h"

#include <algorithm>
#include <cmath>
#include <iostream>

CoupledSjtuMmg3DOFCore::CoupledSjtuMmg3DOFCore(const ShipConfig& ship, const MmgConfig& mmg)
    : ShipCfg(ship), s_(mmg.sjtuMmg)
{
    if (!s_.enabled)
        throw std::runtime_error("CoupledSjtuMmg3DOFCore: CoupledSjtuMmg3DOF.enabled must be true.");
    pi_ = PI;
    setupAddedMassAndInertia();
    U0_ = s_.Fn * std::sqrt(s_.g * s_.Lpp);
    lambdaR_ = s_.HR / std::max(1.0e-8, s_.BR);
    eta_ = s_.Dp / std::max(1.0e-8, s_.HR);

    std::cout << "[CoupledSjtuMmg3DOFCore] SJTU 3-DOF MMG: Lpp=" << s_.Lpp << " m, mass=" << s_.massKg
              << " kg, Fn=" << s_.Fn << ", np=" << s_.npRps << " rps, U0(Fn)=" << U0_ << " m/s\n";
}

void CoupledSjtuMmg3DOFCore::setupAddedMassAndInertia()
{
    Izz_ = s_.massKg * std::pow(s_.kzz_nd * s_.Lpp, 2.0);
    const double massScale = 0.5 * s_.rho * s_.Lpp * s_.Lpp * s_.d;
    const double yawScale = 0.5 * s_.rho * std::pow(s_.Lpp, 4.0) * s_.d;
    mx_ = s_.mx_nd * massScale;
    my_ = s_.my_nd * massScale;
    Jz_ = s_.Jz_nd * yawScale;
}

double CoupledSjtuMmg3DOFCore::ittcCf(double speedU) const
{
    const double uAbs = std::max(1.0e-6, std::abs(speedU));
    const double Re = std::max(1.0e6, s_.rho * uAbs * s_.Lwl / std::max(1.0e-12, s_.nuWater));
    const double lr = std::log10(Re) - 2.0;
    const double den = std::max(1.0e-6, lr * lr);
    return 0.075 / den;
}

double CoupledSjtuMmg3DOFCore::resistanceR(double speedU) const
{
    const double Cf = ittcCf(speedU);
    const double cd = s_.kf * Cf + s_.Cr;
    return cd * 0.5 * s_.rho * speedU * speedU * s_.S_wetted;
}

double CoupledSjtuMmg3DOFCore::referenceSpeedFromFn(double Fn) const
{
    if (std::abs(Fn - s_.Fn) < 1.0e-12)
        return U0_;
    return Fn * std::sqrt(s_.g * s_.Lpp);
}

double CoupledSjtuMmg3DOFCore::deg2rad(double deg) { return deg * PI / 180.0; }
double CoupledSjtuMmg3DOFCore::clamp(double x, double lo, double hi) { return std::max(lo, std::min(hi, x)); }

CoupledSjtuMmg3DOFCore::State CoupledSjtuMmg3DOFCore::toInner(const CoupledSlowState3DOF& s)
{
    State a;
    a.u = s.u;
    a.v = s.v;
    a.r = s.r;
    a.xe = s.xe;
    a.ye = s.ye;
    a.psi = s.psi;
    a.del = s.delta;
    return a;
}

CoupledSlowState3DOF CoupledSjtuMmg3DOFCore::toOuter(const State& s, double t, double Fn)
{
    CoupledSlowState3DOF o;
    o.t = t;
    o.u = s.u;
    o.v = s.v;
    o.r = s.r;
    o.xe = s.xe;
    o.ye = s.ye;
    o.psi = s.psi;
    o.delta = s.del;
    o.U = std::hypot(s.u, s.v);
    o.betaDrift = std::atan2(-s.v, std::max(1e-8, s.u));
    o.Fn = Fn;
    return o;
}

CoupledSjtuMmg3DOFCore::ForceBreakdown CoupledSjtuMmg3DOFCore::evalForces(const State& s) const
{
    ForceBreakdown fb{};
    const double V = std::max(1.0e-8, std::hypot(s.u, s.v));
    const double vp = s.v / V;
    const double rp = s.r * s_.Lpp / V;
    const double beta = std::atan2(-s.v, s.u);
    fb.V = V;
    fb.vp = vp;
    fb.rp = rp;
    fb.beta = beta;

    const double X0_nd = s_.Xvv * vp * vp + s_.Xvr * vp * rp + s_.Xrr * rp * rp;
    const double Y0_nd = s_.Yv * vp + s_.Yr * rp + s_.Yvvv * std::pow(vp, 3.0) + s_.Yvvr * vp * vp * rp
        + s_.Yvrr * vp * rp * rp + s_.Yrrr * std::pow(rp, 3.0);
    const double N0_nd = s_.Nv * vp + s_.Nr * rp + s_.Nvvv * std::pow(vp, 3.0) + s_.Nvvr * vp * vp * rp
        + s_.Nvrr * vp * rp * rp + s_.Nrrr * std::pow(rp, 3.0);

    const double qV = 0.5 * s_.rho * s_.Lpp * s_.Lpp * V * V;
    fb.X0 = qV * X0_nd;
    fb.Y0 = qV * Y0_nd;
    fb.N0 = qV * s_.Lpp * N0_nd;

    fb.Ru = resistanceR(s.u);

    const double betaP = beta - s_.xP_nd * rp;
    fb.wp = s_.wP0 * std::exp(-s_.wakeExpBetaP * betaP * betaP);
    fb.Jp = (s_.npRps > 1.0e-12) ? (s.u * (1.0 - fb.wp) / (s_.npRps * s_.Dp)) : 0.0;
    fb.KT = s_.j0 + s_.j1 * fb.Jp + s_.j2 * fb.Jp * fb.Jp;
    fb.XP = (1.0 - s_.tp) * s_.rho * s_.npRps * s_.npRps * std::pow(s_.Dp, 4.0) * fb.KT;

    const double uP = s.u * ((1.0 - fb.wp) + s_.tauPropInflow * betaP * betaP);

    // Same Jp→0 singularity as legacy MMG inflow; floor only for this subexpression.
    constexpr double kJpMagFloorForPropInflow = 0.02;
    double propTerm = 1.0;
    if (std::abs(fb.Jp) > 1.0e-10)
    {
        const double JpMag = std::max(std::abs(fb.Jp), kJpMagFloorForPropInflow);
        const double inside = std::max(0.0, 1.0 + 8.0 * fb.KT / (pi_ * JpMag * JpMag));
        propTerm = std::sqrt(inside);
    }
    fb.uR = uP * s_.epsilon * propTerm;

    const double gamma = (vp >= 0.0) ? s_.gammaR_plus : s_.gammaR_minus;
    const double vRp = gamma * vp + s_.C_Rv * vp + s_.C_Rrv * rp * vp + s_.C_Rrrv * rp * rp * vp;
    fb.vR = V * vRp;
    fb.UR = std::hypot(fb.uR, fb.vR);
    fb.aR = s.del - std::atan2(-fb.vR, fb.uR);

    const double fAlpha = 6.13 * lambdaR_ / (2.25 + lambdaR_);
    fb.FN = 0.5 * s_.rho * s_.AR * fb.UR * fb.UR * fAlpha * std::sin(fb.aR);

    const double xR = s_.xR_nd * s_.Lpp;
    const double xH = s_.xH_nd * s_.Lpp;
    fb.XR = -(1.0 - s_.tR) * fb.FN * std::sin(s.del);
    fb.YR = -(1.0 + s_.aH) * fb.FN * std::cos(s.del);
    fb.NR = -(xR + s_.aH * xH) * fb.FN * std::cos(s.del);
    return fb;
}

CoupledSjtuMmg3DOFCore::Deriv CoupledSjtuMmg3DOFCore::rhs(
    const State& s,
    const CoupledExternalLoads3DOF& ext,
    double deltaCmdRad,
    double localDt) const
{
    const ForceBreakdown fb = evalForces(s);
    Deriv k{};
    const double Xh = fb.X0 - fb.Ru;
    k.du = ((s_.massKg + my_) * s.v * s.r + Xh + fb.XP + fb.XR + ext.X) / (s_.massKg + mx_);
    k.dv = (-(s_.massKg + mx_) * s.u * s.r + fb.Y0 + fb.YR + ext.Y) / (s_.massKg + my_);
    k.dr = (fb.N0 + fb.NR + ext.N) / (Izz_ + Jz_);
    k.dxe = s.u * std::cos(s.psi) - s.v * std::sin(s.psi);
    k.dye = s.u * std::sin(s.psi) + s.v * std::cos(s.psi);
    k.dpsi = s.r;

    const double deltaMaxRad = deg2rad(s_.maxDeltaDeg);
    const double deltaRateRad = deg2rad(s_.deltaRateDegPerS);
    const double cmd = clamp(deltaCmdRad, -deltaMaxRad, deltaMaxRad);
    const double diff = cmd - s.del;
    k.ddel = (localDt <= 1.0e-14) ? 0.0 : clamp(diff / localDt, -deltaRateRad, +deltaRateRad);
    return k;
}

CoupledSjtuMmg3DOFCore::State CoupledSjtuMmg3DOFCore::rk4Step(
    const State& s,
    const CoupledExternalLoads3DOF& ext,
    double dt,
    double deltaCmdRad) const
{
    auto addScaled = [](const State& a, const Deriv& k, double h) -> State
    {
        State out = a;
        out.u += h * k.du;
        out.v += h * k.dv;
        out.r += h * k.dr;
        out.xe += h * k.dxe;
        out.ye += h * k.dye;
        out.psi += h * k.dpsi;
        out.del += h * k.ddel;
        return out;
    };

    const Deriv k1 = rhs(s, ext, deltaCmdRad, dt);
    const Deriv k2 = rhs(addScaled(s, k1, 0.5 * dt), ext, deltaCmdRad, 0.5 * dt);
    const Deriv k3 = rhs(addScaled(s, k2, 0.5 * dt), ext, deltaCmdRad, 0.5 * dt);
    const Deriv k4 = rhs(addScaled(s, k3, dt), ext, deltaCmdRad, dt);

    State out = s;
    out.u += dt * (k1.du + 2.0 * k2.du + 2.0 * k3.du + k4.du) / 6.0;
    out.v += dt * (k1.dv + 2.0 * k2.dv + 2.0 * k3.dv + k4.dv) / 6.0;
    out.r += dt * (k1.dr + 2.0 * k2.dr + 2.0 * k3.dr + k4.dr) / 6.0;
    out.xe += dt * (k1.dxe + 2.0 * k2.dxe + 2.0 * k3.dxe + k4.dxe) / 6.0;
    out.ye += dt * (k1.dye + 2.0 * k2.dye + 2.0 * k3.dye + k4.dye) / 6.0;
    out.psi += dt * (k1.dpsi + 2.0 * k2.dpsi + 2.0 * k3.dpsi + k4.dpsi) / 6.0;
    out.del += dt * (k1.ddel + 2.0 * k2.ddel + 2.0 * k3.ddel + k4.ddel) / 6.0;

    const double deltaMaxRad = deg2rad(s_.maxDeltaDeg);
    out.del = clamp(out.del, -deltaMaxRad, +deltaMaxRad);
    return out;
}

CoupledSlowState3DOF CoupledSjtuMmg3DOFCore::step(
    const CoupledSlowState3DOF& s,
    const CoupledExternalLoads3DOF& ext,
    double deltaCmd,
    double dt)
{
    const State next = rk4Step(toInner(s), ext, dt, deltaCmd);
    CoupledSlowState3DOF out = toOuter(next, s.t + dt, s.Fn);
    const Deriv dPost = rhs(next, ext, deltaCmd, dt);
    out.rDotRadPerS2 = dPost.dr;
    return out;
}
