#include "CoupledMmg3DOFCore.h"
#include "../../const/Const.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

/// Inner sqrt(1 + ...) in Eq. (3.25) prop-race bracket: 1 + kappa * (sqrt(radicandKt) - 1).
enum class RudderURadPropRaceMode : int
{
    Fig325_KtOverPiJp = 0,         // reference figure: sqrt(1 + K_T / (pi J_P))
    MomentumDisc_8KtOverPiJp2 = 1, // common MMG / actuator disc: sqrt(1 + 8 K_T / (pi J_P^2))
};

// ---------------------------------------------------------------------------
// *** Edit here: pick model + tweak numerators (denominator always pi * J_p^k). ***
constexpr RudderURadPropRaceMode kRudderURadPropRaceMode =
    RudderURadPropRaceMode::MomentumDisc_8KtOverPiJp2;

// Floor on |J_P| when forming denominators (avoid |J|->0 blow-up).
constexpr double kRudderURJpMagFloor = 1e-10;

// Fig325 branch: radicandKt = 1 + kFig325KtNumer * K_T / (pi * J_p^kFig325JDenPower)
constexpr double kFig325KtNumer = 1.0;
constexpr int kFig325JDenPower = 1; // 1 -> .../(pi J_p); 2 -> .../(pi J_p^2)

// Momentum-disc branch: radicandKt = 1 + kMomDiscKtNumer * K_T / (pi * J_p^kMomDiscJDenPower)
constexpr double kMomDiscKtNumer = 8.0;
constexpr int kMomDiscJDenPower = 2;
// ---------------------------------------------------------------------------

inline double rudderURadDenomPiJp(double jpMag, int jDenPower)
{
    double d = PI * jpMag;
    for (int k = 1; k < jDenPower; ++k)
        d *= jpMag;
    return d;
}

} // namespace

void CoupledMmg3DOFCore::initFromCoupledMmg3DOFSettings(const MmgConfig& mmg)
{
    const CoupledMmg3DOFSettings& s = mmg.coupledMmg3DOF;
    pi_ = PI;
    rho_ = s.rho;
    g_ = s.g;
    Lpp_ = s.Lpp;
    Lwl_ = s.Lwl;
    B_ = s.B;
    d_ = s.d;
    Cb_ = s.Cb;
    m_ = s.massKg;
    xG_ = s.xG_m;
    kzz_ = s.kzz_nd;
    Dp_ = s.Dp;
    HR_ = s.HR;
    BR_ = s.BR;
    AR_ = s.AR;
    xP_nd_ = s.xP_nd;
    xR_nd_ = s.xR_nd;
    xH_nd_ = s.xH_nd;
    lR_nd_ = s.lR_nd;
    epsilon_ = s.epsilon;
    kappa_ = s.kappa;
    tR_ = s.tR;
    aH_ = s.aH;
    gammaR_plus_ = s.gammaR_plus;
    gammaR_minus_ = s.gammaR_minus;
    deltaRateDeg_ = s.deltaRateDegPerS;
    maxDeltaDeg_ = s.maxDeltaDeg;
    Fn_ = s.Fn;
    U0_ = Fn_ * std::sqrt(g_ * Lpp_);
    np_rps_ = s.npRps;
    mx_nd_ = s.mx_nd;
    my_nd_ = s.my_nd;
    Jz_nd_ = s.Jz_nd;

    if (s.reuseManoeuvringNondimHullProp)
    {
        Xbb_ = mmg.Hull.Xvv;
        Xbr_ = mmg.Hull.Xvr;
        Xrr_ = mmg.Hull.Xrr;
        Xuu_ = mmg.Hull.Xuu;
        Yb_ = mmg.Hull.Yv;
        Yr_ = mmg.Hull.Yr;
        Ybbb_ = mmg.Hull.Yvvv;
        Ybbr_ = mmg.Hull.Yvvr;
        Ybrr_ = mmg.Hull.Yvrr;
        Yrrr_ = mmg.Hull.Yrrr;
        Nb_ = mmg.Hull.Nv;
        Nr_ = mmg.Hull.Nr;
        Nbbb_ = mmg.Hull.Nvvv;
        Nbbr_ = mmg.Hull.Nvvr;
        Nbrr_ = mmg.Hull.Nvrr;
        Nrrr_ = mmg.Hull.Nrrr;
        wp0_ = mmg.Propeller.wP0;
        C1_ = mmg.Propeller.C1;
        tp_ = mmg.Propeller.tp;
        J0_ = mmg.Propeller.k0;
        J1_ = mmg.Propeller.k1;
        J2_ = mmg.Propeller.k2;
    }
    else
    {
        throw std::runtime_error(
            "CoupledMmg3DOF: reuseManoeuvringNondimHullProp=false is not implemented yet.");
    }

    // Legacy scenario knobs (unused by TDGF coupling today); kept for parity with old header.
    dt_ = 10.0;
    tMax_ = 212100.3203435596424;
    approachTime_ = 0.0;
    stopPsiDeg_ = 540.0;

    std::cout << "[CoupledMmg3DOFCore] CoupledMmg3DOF from case.json: Lpp=" << Lpp_
              << " m, mass=" << m_ << " kg, Fn=" << Fn_ << ", np=" << np_rps_
              << " rps, U0(Fn)=" << U0_ << " m/s\n";
}

void CoupledMmg3DOFCore::initLegacyHardcodedFullScaleS175()
{
    // -------------------------------------------------------------------------
    // Archived defaults: S175 full-scale equivalent (Lpp = 175 m, Froude-scaled
    // time/rudder-rate from a 1:50 tank model). Superseded by
    // Manoeuvring.CoupledMmg3DOF in case.json when enabled=true.
    // -------------------------------------------------------------------------
    pi_ = 3.14159265358979323846;
    rho_ = 1025.0;
    g_ = 9.81;

    Lpp_ = 175.0;
    Lwl_ = 178.2;
    B_ = 25.4;
    d_ = 9.5;
    Cb_ = 0.5720;

    m_ = 24196250.0;
    xG_ = -2.55;
    kzz_ = 0.2690;

    Dp_ = 6.505;
    HR_ = 7.7;
    BR_ = 4.215;
    AR_ = 32.5;

    lambdaR_ = 0.0;
    eta_ = 0.0;

    xP_nd_ = 0.0;
    xR_nd_ = -0.50;
    xH_nd_ = -0.48;
    lR_nd_ = -1.0;

    Xbb_ = -0.0711;
    Xbr_ = 0.0573;
    Xrr_ = 0.0037;

    Yb_ = 0.2137;
    Yr_ = 0.0446;
    Ybbb_ = 2.0080;
    Ybbr_ = 0.3942;
    Ybrr_ = 0.7461;
    Yrrr_ = 0.0326;

    Nb_ = 0.0710;
    Nr_ = -0.0409;
    Nbbb_ = -0.0275;
    Nbbr_ = -0.7811;
    Nbrr_ = -0.0287;
    Nrrr_ = -0.0422;

    mx_nd_ = 0.0044;
    my_nd_ = 0.1299;
    Jz_nd_ = 0.0077;
    Xuu_ = 0.01563;

    mx_ = 0.0;
    my_ = 0.0;
    Jz_ = 0.0;

    epsilon_ = 0.921;
    kappa_ = 0.631;
    tR_ = 0.29;
    aH_ = 0.237;
    gammaR_plus_ = 0.139;
    gammaR_minus_ = 0.088;

    deltaRateDeg_ = 1.697056274847714;

    C1_ = -8.0;
    wp0_ = 0.1684;
    tp_ = 0.175;
    J0_ = 0.2932;
    J1_ = -0.1971;
    J2_ = -0.0481;

    Fn_ = 0.15;
    U0_ = 6.215468606629753;
    np_rps_ = 1.4212846301849607;

    dt_ = 10.0;
    tMax_ = 212100.3203435596424;
    approachTime_ = 0.0;
    stopPsiDeg_ = 540.0;
    maxDeltaDeg_ = 35.0;

    std::cout << "[CoupledMmg3DOFCore] legacy hard-coded S175 full-scale defaults (Lpp=175 m).\n";
}

CoupledMmg3DOFCore::CoupledMmg3DOFCore(const ShipConfig& ship, const MmgConfig& mmg)
    : ShipCfg(ship), mmg_(mmg)
{
    if (!mmg.addedCoeff_defined) compute_addedCoeff();
    if (!mmg.hullCoeff_defined) compute_hullCoeff();
    if (!mmg.propellerCoeff_defined) compute_propellerCoeff();
    if (!mmg.rudderCoeff_defined) compute_rudderCoeff();

    if (mmg_.coupledMmg3DOF.enabled)
        initFromCoupledMmg3DOFSettings(mmg_);
    else
        initLegacyHardcodedFullScaleS175();

    setupModelScaleS175();

    // Optional 4-DOF roll (slow heel). Only via the case.json-driven path; the
    // legacy full-scale path leaves rollEnabled_ = false (3-DOF, unchanged).
    if (mmg_.coupledMmg3DOF.enabled && mmg_.coupledMmg3DOF.maneuverDOF == 4)
    {
        const CoupledMmg3DOFSettings& sc = mmg_.coupledMmg3DOF;
        rollEnabled_ = true;
        GM_        = ShipCfg.Mass.GM;
        I44_       = ShipCfg.Mass.Ixx + sc.rollJxxAdd;
        C44roll_   = m_ * g_ * GM_;
        rollBlin_  = sc.rollB44_lin;
        rollBquad_ = sc.rollB44_quad;
        zH_        = sc.zH_over_d * d_;
        zR_        = sc.zR_over_d * d_;
        zG_        = sc.zG_over_d * d_;
        const double Tn = (C44roll_ > 0.0) ? 2.0 * PI * std::sqrt(I44_ / C44roll_) : 0.0;
        std::cout << "[CoupledMmg3DOFCore] 4-DOF roll ON: I44=" << I44_
                  << " C44=" << C44roll_ << " Tn=" << Tn << " s GM=" << GM_
                  << " zH=" << zH_ << " zR=" << zR_ << " zG=" << zG_
                  << " B44(lin,quad)=(" << rollBlin_ << "," << rollBquad_ << ")\n";
    }
}

double CoupledMmg3DOFCore::referenceSpeedFromFn(double Fn) const
{
    if (std::abs(Fn - Fn_) < 1.0e-12)
        return U0_;
    return Fn * std::sqrt(g_ * Lpp_);
}

void CoupledMmg3DOFCore::setupModelScaleS175()
{
    Izz_ = m_ * std::pow(kzz_ * Lpp_, 2.0);
    lambdaR_ = HR_ / BR_;
    eta_ = Dp_ / HR_;

    const double massScale = 0.5 * rho_ * Lpp_ * Lpp_ * d_;
    const double yawScale = 0.5 * rho_ * std::pow(Lpp_, 4.0) * d_;

    mx_ = mx_nd_ * massScale;
    my_ = my_nd_ * massScale;
    Jz_ = Jz_nd_ * yawScale;
}

double CoupledMmg3DOFCore::deg2rad(double deg) { return deg * PI / 180.0; }
double CoupledMmg3DOFCore::clamp(double x, double lo, double hi) { return std::max(lo, std::min(hi, x)); }

CoupledMmg3DOFCore::State CoupledMmg3DOFCore::toInner(const CoupledSlowState3DOF& s)
{
    State a; a.u=s.u; a.v=s.v; a.r=s.r; a.xe=s.xe; a.ye=s.ye; a.psi=s.psi; a.del=s.delta; return a;
}

CoupledSlowState3DOF CoupledMmg3DOFCore::toOuter(const State& s, double t, double Fn)
{
    CoupledSlowState3DOF o;
    o.t=t; o.u=s.u; o.v=s.v; o.r=s.r; o.xe=s.xe; o.ye=s.ye; o.psi=s.psi; o.delta=s.del;
    o.U=std::hypot(s.u,s.v); o.betaDrift=std::atan2(-s.v, std::max(1e-8,s.u)); o.Fn=Fn; return o;
}

CoupledMmg3DOFCore::ForceBreakdown CoupledMmg3DOFCore::evalForces(const State& s) const
{
    ForceBreakdown fb{};
    const double U = std::max(1.0e-8, std::hypot(s.u, s.v));
    const double beta = std::atan2(-s.v, s.u);
    const double r_nd = s.r * Lpp_ / U;
    fb.U = U; fb.beta = beta; fb.r_nd = r_nd;

    const double q = 0.5 * rho_ * Lpp_ * d_ * U * U;
    const double qN = q * Lpp_;

    const double XH_nd = Xbb_ * beta * beta + Xbr_ * beta * r_nd + Xrr_ * r_nd * r_nd;
    const double YH_nd = Yb_ * beta + Yr_ * r_nd + Ybbb_ * beta * beta * beta + Ybbr_ * beta * beta * r_nd + Ybrr_ * beta * r_nd * r_nd + Yrrr_ * r_nd * r_nd * r_nd;
    const double NH_nd = Nb_ * beta + Nr_ * r_nd + Nbbb_ * beta * beta * beta + Nbbr_ * beta * beta * r_nd + Nbrr_ * beta * r_nd * r_nd + Nrrr_ * r_nd * r_nd * r_nd;

    fb.Ru = 0.5 * rho_ * Lpp_ * d_ * s.u * s.u * Xuu_;
    fb.XH = q * XH_nd - fb.Ru;
    fb.YH = q * YH_nd;
    fb.NH = qN * NH_nd;

    const double betaP = beta - xP_nd_ * r_nd;
    fb.wp = wp0_ * std::exp(C1_ * betaP * betaP);
    // Advance coefficient: literature (MMG) J = (1-w_p) U / (n_p D_p) with U = |u+iv|,
    // not longitudinal u alone — matters once |v| grows during turning.
    fb.Jp = (np_rps_ > 1.0e-12) ? (s.u * (1.0 - fb.wp) / (np_rps_ * Dp_)) : 0.0;
    fb.KT = J0_ + J1_ * fb.Jp + J2_ * fb.Jp * fb.Jp;
    fb.XP = (1.0 - tp_) * rho_ * np_rps_ * np_rps_ * std::pow(Dp_, 4.0) * fb.KT;

    // Rudder axial inflow (Eq. 3.25): u_R = epsilon * u_P * sqrt( eta * bracket^2 + (1-eta) ),
    // bracket = 1 + kappa * (sqrt(radicandKt) - 1). Select radicandKt model in file-top anonymous namespace.
    const double uP = s.u * (1.0 - fb.wp);

    double bracket = 1.0;
    if (std::abs(fb.Jp) > 1.0e-10)
    {
        const double JpMag = std::max(std::abs(fb.Jp), kRudderURJpMagFloor);
        double radicandKt = 1.0;
        switch (kRudderURadPropRaceMode)
        {
        case RudderURadPropRaceMode::Fig325_KtOverPiJp:
            radicandKt = std::max(
                0.0,
                1.0 + fb.KT / rudderURadDenomPiJp(JpMag, kFig325JDenPower));
            
            break;
        case RudderURadPropRaceMode::MomentumDisc_8KtOverPiJp2:
            radicandKt = std::max(
                0.0,
                1.0 + 8.0 * fb.KT / rudderURadDenomPiJp(JpMag, kMomDiscJDenPower));
            break;
        default:
            radicandKt = 1.0;
            break;
        }
        bracket = 1.0 + kappa_ * (std::sqrt(radicandKt) - 1.0);
    }

    const double blendUnderSqrt =
        std::max(0.0, eta_ * bracket * bracket + (1.0 - eta_));
    fb.uR = epsilon_ * uP * std::sqrt(blendUnderSqrt);

    const double betaR = beta - lR_nd_ * r_nd;
    //const double gamma = (betaR >= 0.0) ? gammaR_plus_ : gammaR_minus_;
    const double gamma = (s.v <= 0.0) ? gammaR_plus_ : gammaR_minus_;
    fb.vR = -U * gamma * betaR;
    fb.UR = std::hypot(fb.uR, fb.vR);
    fb.aR = s.del - std::atan2(-fb.vR, fb.uR);

    const double fAlpha = 6.13 * lambdaR_ / (2.25 + lambdaR_);
    fb.FN = 0.5 * rho_ * AR_ * fb.UR * fb.UR * fAlpha * std::sin(fb.aR);

    const double xR = xR_nd_ * Lpp_;
    const double xH = xH_nd_ * Lpp_;
    fb.XR = -(1.0 - tR_) * fb.FN * std::sin(s.del);
    fb.YR = -(1.0 + aH_) * fb.FN * std::cos(s.del);
    fb.NR = -(xR + aH_ * xH) * fb.FN * std::cos(s.del);
    return fb;
}

CoupledMmg3DOFCore::Deriv CoupledMmg3DOFCore::rhs(
    const State& s,
    const CoupledExternalLoads3DOF& ext,
    double deltaCmdRad,
    double localDt) const
{
    const ForceBreakdown fb = evalForces(s);
    Deriv k{};
    k.du = ((m_ + my_) * s.v * s.r + fb.XH + fb.XP + fb.XR + ext.X) / (m_ + mx_);
    k.dv = (-(m_ + mx_) * s.u * s.r + fb.YH + fb.YR + ext.Y) / (m_ + my_);
    k.dr = (fb.NH + fb.NR + ext.N) / (Izz_ + Jz_);
    k.dxe = s.u * std::cos(s.psi) - s.v * std::sin(s.psi);
    k.dye = s.u * std::sin(s.psi) + s.v * std::cos(s.psi);
    k.dpsi = s.r;

    const double deltaMaxRad = deg2rad(maxDeltaDeg_);
    const double deltaRateRad = deg2rad(deltaRateDeg_);
    const double cmd = clamp(deltaCmdRad, -deltaMaxRad, deltaMaxRad);
    const double diff = cmd - s.del;
    k.ddel = (localDt <= 1.0e-14) ? 0.0 : clamp(diff / localDt, -deltaRateRad, +deltaRateRad);
    return k;
}

CoupledMmg3DOFCore::State CoupledMmg3DOFCore::rk4Step(
    const State& s,
    const CoupledExternalLoads3DOF& ext,
    double dt,
    double deltaCmdRad) const
{
    auto addScaled = [](const State& a, const Deriv& k, double h) -> State
    {
        State out = a;
        out.u += h * k.du; out.v += h * k.dv; out.r += h * k.dr;
        out.xe += h * k.dxe; out.ye += h * k.dye; out.psi += h * k.dpsi; out.del += h * k.ddel;
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

    const double deltaMaxRad = deg2rad(maxDeltaDeg_);
    out.del = clamp(out.del, -deltaMaxRad, +deltaMaxRad);
    return out;
}

void CoupledMmg3DOFCore::integrateRoll(double& phi, double& p,
    double YH, double YR, double vdot, double u, double r,
    double Kext, double dt) const
{
    if (!rollEnabled_ || I44_ <= 0.0 || dt <= 0.0) return;

    // Heeling moment (frozen over the slow step). Eq.(4)+(5)+(16) + mean wave heel:
    //   I44*phidd = zG*(my*vdot + mx*u*r) - zH*YH - zG*YR + Kext - C44*sin(phi) - damping
    // Kext = K2 (二阶平均横倾矩). 注意:K2 来自面元解(耐波系),与离心项 -zH*YH
    // 的 MMG 横摇符号约定若不一致,横/斜浪首测时若heel方向反了,把这里的 +Kext
    // 改成 -Kext 即可（单字符）。迎/随浪 K2≈0,不影响。
    const double Mext = zG_ * (my_ * vdot + mx_ * u * r) - zH_ * YH - zG_ * YR + Kext;

    const double Tn = (C44roll_ > 0.0) ? 2.0 * PI * std::sqrt(I44_ / C44roll_) : dt;
    const int nSub = std::max(1,
        static_cast<int>(std::ceil(40.0 * dt / std::max(Tn, 1.0e-9))));
    const double h = dt / nSub;

    auto pdot = [&](double ph, double pp) -> double
    {
        const double Mdamp = -(rollBlin_ * pp + rollBquad_ * std::abs(pp) * pp);
        return (Mext + Mdamp - C44roll_ * std::sin(ph)) / I44_;
    };

    for (int i = 0; i < nSub; ++i)
    {
        const double k1p = p,             k1v = pdot(phi, p);
        const double k2p = p + 0.5*h*k1v, k2v = pdot(phi + 0.5*h*k1p, p + 0.5*h*k1v);
        const double k3p = p + 0.5*h*k2v, k3v = pdot(phi + 0.5*h*k2p, p + 0.5*h*k2v);
        const double k4p = p + h*k3v,     k4v = pdot(phi + h*k3p,     p + h*k3v);
        phi += h * (k1p + 2.0*k2p + 2.0*k3p + k4p) / 6.0;
        p   += h * (k1v + 2.0*k2v + 2.0*k3v + k4v) / 6.0;
    }
}

CoupledSlowState3DOF CoupledMmg3DOFCore::step(
    const CoupledSlowState3DOF& s,
    const CoupledExternalLoads3DOF& ext,
    double deltaCmd,
    double dt)
{
    const State inner = toInner(s);
    const State next = rk4Step(inner, ext, dt, deltaCmd);
    CoupledSlowState3DOF out = toOuter(next, s.t + dt, s.Fn);
    const Deriv dPost = rhs(next, ext, deltaCmd, dt);
    out.rDotRadPerS2 = dPost.dr;

    // 4-DOF: integrate the slow heel driven by the (start-state) hull/rudder
    // lateral loads. Roll does not feed back to u/v/r, so the 3-DOF trajectory
    // above is untouched. Total physical roll = out.phi + seakeeping fast roll.
    if (rollEnabled_)
    {
        const ForceBreakdown fb0 = evalForces(inner);
        const Deriv d0 = rhs(inner, ext, deltaCmd, dt);
        double phi = s.phi, p = s.p;
        integrateRoll(phi, p, fb0.YH, fb0.YR, d0.dv, inner.u, inner.r, ext.K2, dt);
        out.phi = phi;
        out.p   = p;
    }
    return out;
}
