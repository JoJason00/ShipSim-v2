// mmg.cpp
#include "mmg.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

Mmg::Mmg(const ShipConfig& Ship, std::string casePath, const MmgConfig& MmgCfg)
    : ShipCfg(Ship),
    filePath(std::move(casePath)),
    waves(MmgCfg.waves),
    Time(MmgCfg.Time),
    Hull(MmgCfg.Hull),
    Propeller(MmgCfg.Propeller),
    Rudder(MmgCfg.Rudder),
    Added(MmgCfg.Added),
    Fns(MmgCfg.Fn),
    Turning(MmgCfg.Turning),
    Zigzag(MmgCfg.Zigzag),
    turn_on_turningcase(MmgCfg.turn_on_turningcase),
    turn_on_zigzagcase(MmgCfg.turn_on_zigzagcase)
{
    // 旧版经验公式，当前论文验证版本不依赖这些经验系数
    if (!MmgCfg.addedCoeff_defined)      compute_addedCoeff();
    if (!MmgCfg.hullCoeff_defined)       compute_hullCoeff();
    if (!MmgCfg.propellerCoeff_defined)  compute_propellerCoeff();
    if (!MmgCfg.rudderCoeff_defined)     compute_rudderCoeff();

    setupModelScaleS175();

    if (Time.dt <= 0.0)
        Time.dt = dt_;

    if (Fns.empty())
        Fns.push_back(Fn_);

    if (std::abs(Turning.rudder_angle) < 1.0e-12)
        Turning.rudder_angle = deg2rad(maxDeltaDeg_);

    if (Turning.TurningCaseCicle <= 0.0)
        Turning.TurningCaseCicle = stopPsiDeg_ / 360.0;

    if (std::abs(Zigzag.rudder_angle) < 1.0e-12)
        Zigzag.rudder_angle = deg2rad(10.0);

    if (Zigzag.ZigzagCaseCicle <= 0)
        Zigzag.ZigzagCaseCicle = 4;
}

void Mmg::setupModelScaleS175()
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

double Mmg::deg2rad(double deg)
{
    return deg * 3.14159265358979323846 / 180.0;
}

double Mmg::rad2deg(double rad)
{
    return rad * 180.0 / 3.14159265358979323846;
}

double Mmg::clamp(double x, double lo, double hi)
{
    return std::max(lo, std::min(hi, x));
}

void Mmg::compute_hullCoeff() {}
void Mmg::compute_propellerCoeff() {}
void Mmg::compute_rudderCoeff() {}
void Mmg::compute_addedCoeff() {}

CaseContext_mmg Mmg::buildCaseContext(const int i_wave, const double Fn)
{
    CaseContext_mmg ctx{};

    ctx.i_wave = i_wave;

    if (!waves.empty())
        ctx.wave = waves.at(i_wave);
    else
        ctx.wave = nullptr;

    ctx.reg = std::dynamic_pointer_cast<RegularWave>(ctx.wave);

    // 当前验证版继续用当前这版数值时间步
    ctx.dt = Time.dt > 0.0 ? Time.dt : dt_;

    ctx.Fn = Fn > 0.0 ? Fn : Fn_;

    // 对 Fn=0.15 保持论文给定的 U0=0.879
    if (std::abs(ctx.Fn - Fn_) < 1.0e-12)
        ctx.U = U0_;
    else
        ctx.U = ctx.Fn * std::sqrt(g_ * Lpp_);

    ctx.A.setIdentity();
    ctx.A_lu.compute(ctx.A);

    return ctx;
}

Mmg::ForceBreakdown Mmg::evalForces(double /*t*/, const State& s, const CaseContext_mmg& /*ctx*/) const
{
    ForceBreakdown fb{};

    const double U = std::max(1.0e-8, std::hypot(s.u, s.v));
    const double beta = std::atan2(-s.v, s.u); // beta = atan(-v/u)
    const double r_nd = s.r * Lpp_ / U;        // r' = rL/U

    fb.U = U;
    fb.beta = beta;
    fb.r_nd = r_nd;

    // -------------------------
    // Hull
    // -------------------------
    const double q = 0.5 * rho_ * Lpp_ * d_ * U * U;
    const double qN = q * Lpp_;

    const double XH_nd =
        Xbb_ * beta * beta +
        Xbr_ * beta * r_nd +
        Xrr_ * r_nd * r_nd;

    const double YH_nd =
        Yb_ * beta +
        Yr_ * r_nd +
        Ybbb_ * beta * beta * beta +
        Ybbr_ * beta * beta * r_nd +
        Ybrr_ * beta * r_nd * r_nd +
        Yrrr_ * r_nd * r_nd * r_nd;

    const double NH_nd =
        Nb_ * beta +
        Nr_ * r_nd +
        Nbbb_ * beta * beta * beta +
        Nbbr_ * beta * beta * r_nd +
        Nbrr_ * beta * r_nd * r_nd +
        Nrrr_ * r_nd * r_nd * r_nd;

    fb.Ru = 0.5 * rho_ * Lpp_ * d_ * s.u * s.u * Xuu_;

    fb.XH = q * XH_nd - fb.Ru;
    fb.YH = q * YH_nd;
    fb.NH = qN * NH_nd;

    // -------------------------
    // Propeller
    // -------------------------
    const double betaP = beta - xP_nd_ * r_nd;
    fb.wp = wp0_ * std::exp(C1_ * betaP * betaP);

    fb.Jp = (np_rps_ > 1.0e-12) ? (s.u * (1.0 - fb.wp) / (np_rps_ * Dp_)) : 0.0;
    fb.KT = J0_ + J1_ * fb.Jp + J2_ * fb.Jp * fb.Jp;
    fb.XP = (1.0 - tp_) * rho_ * np_rps_ * np_rps_ * std::pow(Dp_, 4.0) * fb.KT;

    // -------------------------
    // Rudder inflow
    // -------------------------
    const double wR = 1.0 - epsilon_ * (1.0 - fb.wp);
    const double U0p = s.u * (1.0 - wR);

    // Fujii–Nomoto inflow factor diverges as |Jp|→0; floor |Jp| for this subexpression only.
    constexpr double kJpMagFloorForPropInflow = 0.02;
    double propTerm = 1.0;
    if (std::abs(fb.Jp) > 1.0e-10)
    {
        const double JpMag = std::max(std::abs(fb.Jp), kJpMagFloorForPropInflow);
        const double inside = std::max(0.0, 1.0 + 8.0 * fb.KT / (pi_ * JpMag * JpMag));
        propTerm = 1.0 + kappa_ * (std::sqrt(inside) - 1.0);
    }

    const double URP = (1.0 - fb.wp) * s.u * epsilon_ * propTerm;
    fb.uR = std::sqrt(std::max(0.0, eta_ * URP * URP + (1.0 - eta_) * U0p * U0p));

    const double betaR = beta - lR_nd_ * r_nd;
    const double gamma = (betaR >= 0.0) ? gammaR_plus_ : gammaR_minus_;
    fb.vR = U * gamma * betaR;

    fb.UR = std::hypot(fb.uR, fb.vR);
    fb.aR = s.del - std::atan2(-fb.vR, fb.uR);

    const double fAlpha = 6.13 * lambdaR_ / (2.25 + lambdaR_);
    fb.FN = 0.5 * rho_ * AR_ * fb.UR * fb.UR * fAlpha * std::sin(fb.aR);

    // -------------------------
    // Rudder force / moment
    // -------------------------
    const double xR = xR_nd_ * Lpp_;
    const double xH = xH_nd_ * Lpp_;

    fb.XR = -(1.0 - tR_) * fb.FN * std::sin(s.del);
    fb.YR = -(1.0 + aH_) * fb.FN * std::cos(s.del);
    fb.NR = -(xR + aH_ * xH) * fb.FN * std::cos(s.del);


    return fb;
}

Mmg::Deriv Mmg::rhs(double t, const State& s, const CaseContext_mmg& ctx,
    double deltaCmdRad, double localDt,
    ForceBreakdown* fbPtr) const
{
    ForceBreakdown fb = evalForces(t, s, ctx);
    if (fbPtr) *fbPtr = fb;

    Deriv k{};

    k.du = ((m_ + my_) * s.v * s.r + fb.XH + fb.XP + fb.XR) / (m_ + mx_);
    k.dv = (-(m_ + mx_) * s.u * s.r + fb.YH + fb.YR) / (m_ + my_);
    k.dr = (fb.NH + fb.NR) / (Izz_ + Jz_);

    k.dxe = s.u * std::cos(s.psi) - s.v * std::sin(s.psi);
    k.dye = s.u * std::sin(s.psi) + s.v * std::cos(s.psi);
    k.dpsi = s.r;

    const double deltaMaxRad = deg2rad(maxDeltaDeg_);
    const double deltaRateRad = deg2rad(deltaRateDeg_);
    const double cmd = clamp(deltaCmdRad, -deltaMaxRad, deltaMaxRad);
    const double diff = cmd - s.del;

    if (localDt <= 1.0e-14)
    {
        k.ddel = 0.0;
    }
    else
    {
        k.ddel = clamp(diff / localDt, -deltaRateRad, +deltaRateRad);
    }

    return k;
}

Mmg::State Mmg::rk4Step(double t, const State& s, double dt,
    const CaseContext_mmg& ctx, double deltaCmdRad) const
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

    const Deriv k1 = rhs(t, s, ctx, deltaCmdRad, dt, nullptr);
    const Deriv k2 = rhs(t + 0.5 * dt, addScaled(s, k1, 0.5 * dt), ctx, deltaCmdRad, 0.5 * dt, nullptr);
    const Deriv k3 = rhs(t + 0.5 * dt, addScaled(s, k2, 0.5 * dt), ctx, deltaCmdRad, 0.5 * dt, nullptr);
    const Deriv k4 = rhs(t + dt, addScaled(s, k3, dt), ctx, deltaCmdRad, dt, nullptr);

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

void Mmg::writeCsvHeader(std::ofstream& ofs) const
{
    ofs << "t,u,v,r,x,y,psi,delta\n";
}

void Mmg::writeCsvRow(std::ofstream& ofs, double t, const State& s,
    const ForceBreakdown& /*fb*/) const
{
    ofs << std::fixed << std::setprecision(10)
        << t << ","
        << s.u << ","
        << s.v << ","
        << s.r << ","
        << s.xe / Lpp_ << ","
        << s.ye / Lpp_ << ","
        << s.psi << ","
        << s.del
        << "\n";
}

void Mmg::run()
{
    const int n_cases = waves.empty() ? 1 : static_cast<int>(waves.size());

    const bool turn_left = true;
    const bool turn_right = false;

    for (int iFn = 0; iFn < static_cast<int>(Fns.size()); ++iFn)
    {
        const double Fn = Fns[iFn];
        for (int i_case = 0; i_case < n_cases; ++i_case)
        {
            CaseContext_mmg ctx = buildCaseContext(i_case, Fn);
            ctx.i_wave = i_case;

            if (turn_on_turningcase)
            {
                runTurningCase(ctx, turn_left);
                runTurningCase(ctx, turn_right);
            }

            if (turn_on_zigzagcase)
                runZigzagCase(ctx);
        }
    }
}

void Mmg::runTurningCase(const CaseContext_mmg& ctx, const bool direction)
{
    State s{};
    s.u = ctx.U;

    double t = 0.0;

    double psi_acc = 0.0;
    double psi_prev = s.psi;

    const double delta_mag = (std::abs(Turning.rudder_angle) > 1.0e-12)
        ? std::abs(Turning.rudder_angle)
        : deg2rad(maxDeltaDeg_);

    const double delta_target = direction ? delta_mag : -delta_mag;

    const double turningCircle = (Turning.TurningCaseCicle > 0.0)
        ? Turning.TurningCaseCicle
        : (stopPsiDeg_ / 360.0);

    const double turn_angle = turningCircle * 2.0 * pi_;

    std::filesystem::create_directories(filePath + "/Mmg_turning");
    std::ostringstream oss;
    oss << filePath << "/Mmg_turning/"
        << std::fixed << std::setprecision(3) << delta_target
        << "_Fn_" << ctx.Fn << "_wave_" << ctx.i_wave
        << "_turning.csv";

    std::ofstream fout(oss.str());
    if (!fout)
        throw std::runtime_error("Cannot open output file: " + oss.str());

    fout << std::fixed << std::setprecision(10);
    writeCsvHeader(fout);

    bool maneuverStarted = false;
    double t0 = 0.0;
    double x0 = 0.0;
    double y0 = 0.0;
    double psi0 = 0.0;

    while (t <= tMax_)
    {
        if (!maneuverStarted && t >= approachTime_)
        {
            maneuverStarted = true;
            t0 = t;
            x0 = s.xe;
            y0 = s.ye;
            psi0 = s.psi;
            psi_prev = s.psi;
        }

        const double cmd = maneuverStarted ? delta_target : 0.0;

        s = rk4Step(t, s, ctx.dt, ctx, cmd);
        t += ctx.dt;

        if (maneuverStarted)
        {
            State sout = s;
            sout.xe -= x0;
            sout.ye -= y0;
            sout.psi -= psi0;

            ForceBreakdown fb = evalForces(t, s, ctx);
            writeCsvRow(fout, t - t0, sout, fb);

            const double dpsi = std::atan2(std::sin(s.psi - psi_prev), std::cos(s.psi - psi_prev));
            psi_acc += dpsi;
            psi_prev = s.psi;

            if (std::abs(psi_acc) >= turn_angle)
                break;
        }
    }
}

void Mmg::runZigzagCase(const CaseContext_mmg& ctx)
{
    State s{};
    s.u = ctx.U;

    const int zig_cycle = (Zigzag.ZigzagCaseCicle > 0) ? Zigzag.ZigzagCaseCicle : 4;
    const double zig_angle = (std::abs(Zigzag.rudder_angle) > 1.0e-12)
        ? std::abs(Zigzag.rudder_angle)
        : deg2rad(10.0);

    int switch_count = 0;
    double t = 0.0;
    double delta_cmd = 0.0;

    std::filesystem::create_directories(filePath + "/Mmg_zigzag");
    std::ostringstream oss;
    oss << filePath << "/Mmg_zigzag/"
        << std::fixed << std::setprecision(3) << zig_angle
        << "_Fn_" << ctx.Fn << "_wave_" << ctx.i_wave
        << "_zigzag.csv";

    std::ofstream fout(oss.str());
    if (!fout)
        throw std::runtime_error("Cannot open output file: " + oss.str());

    fout << std::fixed << std::setprecision(10);
    writeCsvHeader(fout);

    bool maneuverStarted = false;
    double t0 = 0.0;
    double x0 = 0.0;
    double y0 = 0.0;
    double psi0 = 0.0;

    while (t <= tMax_)
    {
        if (!maneuverStarted && t >= approachTime_)
        {
            maneuverStarted = true;
            t0 = t;
            x0 = s.xe;
            y0 = s.ye;
            psi0 = s.psi;
            delta_cmd = +zig_angle;
        }

        const double cmd = maneuverStarted ? delta_cmd : 0.0;

        s = rk4Step(t, s, ctx.dt, ctx, cmd);
        t += ctx.dt;

        if (maneuverStarted)
        {
            State sout = s;
            sout.xe -= x0;
            sout.ye -= y0;
            sout.psi -= psi0;

            ForceBreakdown fb = evalForces(t, s, ctx);
            writeCsvRow(fout, t - t0, sout, fb);

            if (delta_cmd > 0.0 && sout.psi >= +zig_angle)
            {
                delta_cmd = -zig_angle;
                ++switch_count;
            }
            else if (delta_cmd < 0.0 && sout.psi <= -zig_angle)
            {
                delta_cmd = +zig_angle;
                ++switch_count;
            }

            if (switch_count >= zig_cycle)
                break;
        }
    }
}