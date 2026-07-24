//#include "mmg.h"
//#include "const/Const.h"
//
//#include <fstream>
//#include <iomanip>
//#include <filesystem>
//#include <algorithm>
//#include <sstream>
//#include <cmath>
//#include <iostream>
//
//Mmg::Mmg(const ShipConfig& Ship, std::string casePath, const MmgConfig& Mmg)
//    : ShipCfg(Ship),
//    filePath(std::move(casePath)),
//    waves(Mmg.waves),
//    Time(Mmg.Time),
//    Hull(Mmg.Hull),
//    Propeller(Mmg.Propeller),
//    Rudder(Mmg.Rudder),
//    Added(Mmg.Added),
//    Fns(Mmg.Fn),
//    Turning(Mmg.Turning),
//    Zigzag(Mmg.Zigzag),
//    turn_on_turningcase(Mmg.turn_on_turningcase),
//    turn_on_zigzagcase(Mmg.turn_on_zigzagcase)
//{
//    // 先按你原接口保留逻辑，但静水操纵核心改成 Fortran 路线
//    if (!Mmg.addedCoeff_defined)
//        compute_addedCoeff();
//
//    compute_normMass();
//
//    if (!Mmg.hullCoeff_defined)
//        compute_hullCoeff();
//
//    // 这两个即便不严格需要，保留接口，避免以后扩展时断掉
//    if (!Mmg.propellerCoeff_defined)
//        compute_propellerCoeff();
//
//    if (!Mmg.rudderCoeff_defined)
//        compute_rudderCoeff();
//}
//
//void Mmg::run()
//{
//    const int n_Fns = static_cast<int>(Fns.size());
//    const int n_cases = waves.empty() ? 1 : static_cast<int>(waves.size());
//
//    const bool turn_left = true;
//    const bool turn_right = false;
//
//    for (int iFn = 0; iFn < n_Fns; ++iFn)
//    {
//        const double Fn = Fns[iFn];
//
//        for (int i_case = 0; i_case < n_cases; ++i_case)
//        {
//            CaseContext_mmg ctx = buildCaseContext(i_case, Fn);
//            ctx.i_wave = i_case;
//
//            if (turn_on_turningcase)
//            {
//                runTurningCase(ctx, turn_left);
//                runTurningCase(ctx, turn_right);
//            }
//
//            if (turn_on_zigzagcase)
//                runZigzagCase(ctx);
//        }
//    }
//}
//
//void Mmg::runTurningCase(const CaseContext_mmg& ctx, const bool direction)
//{
//    Vec6 Y{};
//    Y.setZero();
//
//    // 状态定义：
//    // Y(0)=u, Y(1)=v, Y(2)=r, Y(3)=x, Y(4)=y, Y(5)=psi
//    Y(0) = ctx.U;
//
//    double psi_acc = 0.0;
//    double t = 0.0;
//    double psi_prev = Y(5);
//    double delta = 0.0;
//
//    const double delta_target = direction ? Turning.rudder_angle : -Turning.rudder_angle;
//    const double delta_rate = Turning.rudder_rate;
//
//    std::filesystem::create_directories(filePath + "/Mmg_turning");
//
//    std::ostringstream oss;
//    oss << filePath << "/Mmg_turning/"
//        << std::fixed << std::setprecision(3) << delta_target
//        << "_Fn_" << ctx.Fn
//        << "_wave_" << ctx.i_wave
//        << "_turning.csv";
//
//    std::ofstream fout(oss.str());
//    fout << std::fixed << std::setprecision(10);
//    fout << "t,u,v,r,x,y,psi,delta\n";
//
//    const double turn_angle = Turning.TurningCaseCicle * 2.0 * PI;
//    const double L = ShipCfg.Geometry.Length;
//
//    while (std::abs(psi_acc) < turn_angle)
//    {
//        if (delta < delta_target)
//            delta = std::min(delta + delta_rate * ctx.dt, delta_target);
//        else if (delta > delta_target)
//            delta = std::max(delta - delta_rate * ctx.dt, delta_target);
//
//        stepRK4Turning(Y, t, ctx.dt, ctx, delta);
//        t += ctx.dt;
//
//        fout << t
//            << "," << Y(0)
//            << "," << Y(1)
//            << "," << Y(2)
//            << "," << Y(3) / L
//            << "," << Y(4) / L
//            << "," << Y(5)
//            << "," << delta
//            << "\n";
//
//        const double psi_now = Y(5);
//        const double dpsi = std::atan2(std::sin(psi_now - psi_prev), std::cos(psi_now - psi_prev));
//        psi_acc += dpsi;
//        psi_prev = psi_now;
//    }
//}
//
//void Mmg::runZigzagCase(const CaseContext_mmg& ctx)
//{
//    Vec6 Y{};
//    Y.setZero();
//    Y(0) = ctx.U;
//
//    const int zig_cycle = static_cast<int>(Zigzag.ZigzagCaseCicle);
//    const double zig_angle = Zigzag.rudder_angle;
//    const double zig_rate = Zigzag.rudder_rate;
//
//    int switch_count = 0;
//    double delta = 0.0;
//    double delta_cmd = zig_angle;
//    double psi_trigger = zig_angle;
//    double t = 0.0;
//
//    const double L = ShipCfg.Geometry.Length;
//
//    std::filesystem::create_directories(filePath + "/Mmg_zigzag");
//
//    std::ostringstream oss;
//    oss << filePath << "/Mmg_zigzag/"
//        << std::fixed << std::setprecision(3) << zig_angle
//        << "_Fn_" << ctx.Fn
//        << "_wave_" << ctx.i_wave
//        << "_zigzag.csv";
//
//    std::ofstream fout(oss.str());
//    fout << std::fixed << std::setprecision(10);
//    fout << "t,u,v,r,x,y,psi,delta\n";
//
//    while (switch_count < zig_cycle)
//    {
//        if (delta < delta_cmd)
//            delta = std::min(delta + zig_rate * ctx.dt, delta_cmd);
//        else if (delta > delta_cmd)
//            delta = std::max(delta - zig_rate * ctx.dt, delta_cmd);
//
//        stepRK4Turning(Y, t, ctx.dt, ctx, delta);
//        t += ctx.dt;
//
//        fout << t
//            << "," << Y(0)
//            << "," << Y(1)
//            << "," << Y(2)
//            << "," << Y(3) / L
//            << "," << Y(4) / L
//            << "," << Y(5)
//            << "," << delta
//            << "\n";
//
//        const bool reached_trigger =
//            (psi_trigger > 0.0 && Y(5) >= psi_trigger) ||
//            (psi_trigger < 0.0 && Y(5) <= psi_trigger);
//
//        if (reached_trigger)
//        {
//            delta_cmd = -delta_cmd;
//            psi_trigger = -psi_trigger;
//            switch_count++;
//        }
//    }
//}
//
//void Mmg::compute_normMass()
//{
//    const double L = ShipCfg.Geometry.Length;
//    const double d = ShipCfg.Geometry.Draft;
//
//    const double denom_force = 0.5 * rho * L * L * d;
//    const double denom_moment = 0.5 * rho * d * std::pow(L, 4);
//
//    MassND.Mm1 = ShipCfg.Mass.Mass / denom_force;
//    MassND.Mm4 = ShipCfg.Mass.Izz / denom_moment;
//
//    MassND.mx = Added.mx / denom_force;
//    MassND.my = Added.my / denom_force;
//    MassND.Jz = Added.Jz / denom_moment;
//}
//
//void Mmg::compute_hullCoeff()
//{
//    // 这里尽量按 Fortran 的 F_coef 来
//    const double L = ShipCfg.Geometry.Length;
//    const double B = ShipCfg.Geometry.Breadth;
//    const double d = ShipCfg.Geometry.Draft;
//    const double CB = ShipCfg.Geometry.CB;
//    const double trm = ShipCfg.Geometry.Trim;
//
//    const double kad = 2.0 * d / L;
//
//    Hull.Xvv = 1.15 * CB / (L / B) - 0.18;
//    Hull.Xvr = -1.91 * CB / (L / B) + 0.08 + MassND.my;
//    Hull.Xrr = -0.085 * CB / (L / B) + 0.008;
//    Hull.Xuu = 0.01563;
//
//    Hull.Yv = (0.5 * PI * kad + 1.4 * CB / (L / B)) * (1.0 + 0.54 * std::pow(-trm / d, 2.0));
//    Hull.Yr = 0.5 * CB / (L / B) * (1.0 + 1.82 * std::pow(-trm / d, 2.0)) + MassND.mx;
//    Hull.Nv = kad * (1.0 - 0.85 * (-trm / d));
//    Hull.Nr = (-0.54 * kad + kad * kad) * (1.0 + 0.33 * (-trm / d));
//
//    Hull.Yvvv = 0.185 * (L / B) + 0.48;
//    Hull.Yrrr = 0.069 * (-trm / d) - 0.051;
//    Hull.Yvrr = 0.26 * (1.0 - CB) * (L / B) + 0.11;
//    Hull.Yvvr = 0.97 * (-trm / d) / CB - 0.75;
//
//    Hull.Nvvv = -0.69 * CB + 0.66;
//    Hull.Nrrr = 0.25 * CB / (L / B) - 0.056;
//    Hull.Nvrr = 0.075 * (1.0 - CB) * (L / B) - 0.098;
//    Hull.Nvvr = 1.55 * CB / (L / B) - 0.76;
//}
//
//void Mmg::compute_propellerCoeff()
//{
//    // 保留接口。
//    // Fortran 注释里的经验式：
//    // wp_0=0.7*CP-0.18
//    // t_thrust=0.5*CP-0.12
//    // 这里只在“未定义”时给一个经验估计，避免空函数完全失效。
//    const double CP = 0.605; // 你当前结构里没有单独存 CP，先用 Fortran S175 的值
//    Propeller.wP0 = 0.7 * CP - 0.18;
//    Propeller.tp = 0.5 * CP - 0.12;
//}
//
//void Mmg::compute_rudderCoeff()
//{
//    // 保留接口。
//    // Fortran 注释里的经验式：
//    // aH=(0.6784-1.3374*CB+1.8891*CB**2)
//    // xH=-(0.4+0.1*CB)*Lpp
//    // t_R=1-(0.7382-0.0539*CB+0.1755*CB**2)
//    // epR=-156.2*(CB*Bread/Lpp)^2+41.6*(CB*Bread/Lpp)-1.76
//    // kP=0.55/epR
//    const double CB = ShipCfg.Geometry.CB;
//    const double B = ShipCfg.Geometry.Breadth;
//    const double L = ShipCfg.Geometry.Length;
//
//    Rudder.aH = 0.6784 - 1.3374 * CB + 1.8891 * CB * CB;
//    Rudder.xH = -(0.4 + 0.1 * CB) * L;
//    Rudder.tR = 1.0 - (0.7382 - 0.0539 * CB + 0.1755 * CB * CB);
//
//    Rudder.epsilonP = -156.2 * std::pow(CB * B / L, 2.0) + 41.6 * (CB * B / L) - 1.76;
//    if (std::abs(Rudder.epsilonP) > 1.0e-12)
//        Rudder.kP = 0.55 / Rudder.epsilonP;
//}
//
//void Mmg::compute_addedCoeff()
//{
//    // 维持你原来的经验式，得到有量纲附加质量，再在 compute_normMass 里无量纲化
//    const double L = ShipCfg.Geometry.Length;
//    const double B = ShipCfg.Geometry.Breadth;
//    const double d = ShipCfg.Geometry.Draft;
//    const double CB = ShipCfg.Geometry.CB;
//    const double M = ShipCfg.Mass.Mass;
//
//    Added.mx = 0.01 * M * (0.398 + 11.97 * CB * (1.0 + 3.73 * d / B)
//        - 2.89 * CB * (L / B) * (1.0 + 1.13 * d / B)
//        + 0.175 * CB * (L / B) * (L / B) * (1.0 + 0.541 * d / B)
//        - 1.107 * (L / B) * (d / B));
//
//    Added.my = M * (0.882 - 0.54 * CB * (1.0 - 1.6 * d / B)
//        - 0.156 * (L / B) * (1.0 - 0.673 * CB)
//        + 0.826 * (d / B) * (L / B) * (1.0 - 0.678 * d / B)
//        - 0.638 * CB * (d / B) * (L / B) * (1.0 - 0.669 * d / B));
//
//    Added.Jz = std::pow(0.01 * (33.0 - 76.85 * CB * (1.0 - 0.784 * CB)
//        + 3.43 * (L / B) * (1.0 - 0.63 * CB)), 2.0) * M * L * L;
//}
//
//CaseContext_mmg Mmg::buildCaseContext(const int i_wave, const double Fn)
//{
//    CaseContext_mmg ctx{};
//
//    if (!waves.empty())
//    {
//        ctx.wave = waves.at(i_wave);
//        ctx.reg = std::dynamic_pointer_cast<RegularWave>(ctx.wave);
//    }
//    else
//    {
//        ctx.wave = nullptr;
//        ctx.reg = nullptr;
//    }
//
//    ctx.dt = Time.dt;
//    ctx.Fn = Fn;
//    ctx.U = Fn * std::sqrt(G * ShipCfg.Geometry.Length);
//
//    return ctx;
//}
//
//void Mmg::stepRK4Turning(Vec6& Y, double t, double dt, const CaseContext_mmg& ctx, const double rudder_angle)
//{
//    const Vec6 k1 = stateDerivativeTurning(t, Y, ctx, rudder_angle);
//    const Vec6 k2 = stateDerivativeTurning(t + 0.5 * dt, Y + 0.5 * dt * k1, ctx, rudder_angle);
//    const Vec6 k3 = stateDerivativeTurning(t + 0.5 * dt, Y + 0.5 * dt * k2, ctx, rudder_angle);
//    const Vec6 k4 = stateDerivativeTurning(t + dt, Y + dt * k3, ctx, rudder_angle);
//
//    Y += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
//}
//
//HullForceND Mmg::HullForceFortranND(const Vec6& Y) const
//{
//    HullForceND force{};
//
//    const double u = Y(0);
//    const double v = Y(1);
//    const double r = Y(2);
//
//    const double L = ShipCfg.Geometry.Length;
//    const double U = std::hypot(u, v);
//    const double eps = 1.0e-10;
//
//    if (U < eps)
//        return force;
//
//    const double beta = -std::asin(std::clamp(v / U, -1.0, 1.0));
//    const double r_nd = r * L / U;
//
//    // 按 Fortran 的 F_hull
//    force.X =
//        MassND.Mm1 * r * v * L / (U * U)
//        + Hull.Xvv * beta * beta
//        + Hull.Xrr * r_nd * r_nd
//        + Hull.Xvr * beta * r_nd;
//
//    force.Y =
//        Hull.Yv * beta
//        + Hull.Yr * r_nd
//        - MassND.Mm1 * r * u * L / (U * U)
//        + Hull.Yvvv * std::pow(beta, 3.0)
//        + Hull.Yvvr * r_nd * beta * beta
//        + Hull.Yrrr * std::pow(r_nd, 3.0)
//        + Hull.Yvrr * beta * r_nd * r_nd;
//
//    force.N =
//        Hull.Nv * beta
//        + Hull.Nr * r_nd
//        + Hull.Nvvv * std::pow(beta, 3.0)
//        + Hull.Nvvr * r_nd * beta * beta
//        + Hull.Nvrr * beta * r_nd * r_nd
//        + Hull.Nrrr * std::pow(r_nd, 3.0);
//
//    return force;
//}
//
//double Mmg::ResistanceForce(const Vec6& Y) const
//{
//    const double L = ShipCfg.Geometry.Length;
//    const double d = ShipCfg.Geometry.Draft;
//    const double U = std::hypot(Y(0), Y(1));
//
//    const double eps = 1.0e-10;
//    if (U < eps)
//        return 0.0;
//
//    // 按 Fortran 的 X_R
//    const double cf =
//        0.075 / std::pow(std::log10(L * U * 1.0e6) - 2.03, 2.0)
//        + 0.0001
//        + 0.00054;
//
//    const double Xuu = 9225.0 * cf / (L * d);
//    return Xuu * 0.5 * rho * L * d * U * U;
//}
//
//double Mmg::PropellerForce(const Vec6& Y, double& wP, double& Jp, double& KT) const
//{
//    const double u = Y(0);
//    const double v = Y(1);
//    const double r = Y(2);
//
//    const double L = ShipCfg.Geometry.Length;
//    const double U = std::hypot(u, v);
//    const double eps = 1.0e-10;
//
//    wP = Propeller.wP0;
//    Jp = 0.0;
//    KT = 0.0;
//
//    if (U < eps || std::abs(Propeller.np) < eps || std::abs(Propeller.Dp) < eps)
//        return 0.0;
//
//    const double beta = -std::asin(std::clamp(v / U, -1.0, 1.0));
//    const double r_nd = r * L / U;
//
//    // 按 Fortran 的 X_P
//    wP = Propeller.wP0 * std::exp(-8.0 * std::pow(beta + r_nd, 2.0));
//    const double uP = u * (1.0 - wP);
//
//    Jp = uP / (Propeller.np * Propeller.Dp);
//    KT = Propeller.k0 + Propeller.k1 * Jp + Propeller.k2 * Jp * Jp;
//
//    return (1.0 - Propeller.tp) * rho * std::pow(Propeller.Dp, 4) * KT * Propeller.np * Propeller.np;
//}
//
//RudderForceDim Mmg::RudderForce(const Vec6& Y, const double rudder_angle, const double Jp, const double KT) const
//{
//    RudderForceDim force{};
//
//    const double u = Y(0);
//    const double v = Y(1);
//    const double r = Y(2);
//
//    const double L = ShipCfg.Geometry.Length;
//    const double eps = 1.0e-10;
//
//    if (std::abs(Jp) < eps)
//        return force;
//
//    // 按 Fortran rud_forcea:
//    // ur = epR * u * (1 - wp_0) * sqrt(1 + kP*(8*KT0/(pi*J_p^2)))
//    double ur_term = 1.0 + Rudder.kP * (8.0 * KT / (PI * Jp * Jp));
//    ur_term = std::max(ur_term, 0.0);
//
//    const double ur = Rudder.epsilonP * u * (1.0 - Propeller.wP0) * std::sqrt(ur_term);
//
//    const double gammaR = (v <= 0.0) ? Rudder.gammaR_minus : Rudder.gammaR_plus;
//    const double vr = gammaR * (v + r * Rudder.lR * L);
//
//    const double alphaR = rudder_angle - std::atan2(-vr, ur);
//
//    // 当前输入直接给了 f_alpha = 2.747，对应 Fortran 里 fa=6.13*la/(2.25+la)
//    const double fa = Rudder.f_alpha;
//    const double FN = 0.5 * rho * Rudder.AR * fa * (ur * ur + vr * vr) * std::sin(alphaR);
//
//    force.X = -(1.0 - Rudder.tR) * FN * std::sin(rudder_angle);
//    force.Y = -(1.0 + Rudder.aH) * FN * std::cos(rudder_angle);
//    force.N = -(Rudder.xR + Rudder.aH * Rudder.xH) * FN * std::cos(rudder_angle);
//
//    return force;
//}
//
//Vec6 Mmg::stateDerivativeTurning(double t, const Vec6& Y, const CaseContext_mmg& ctx, const double rudder_angle)
//{
//    (void)t;
//    (void)ctx;
//
//    Vec6 dY{};
//    dY.setZero();
//
//    // 状态定义：
//    // Y(0)=u, Y(1)=v, Y(2)=r, Y(3)=x, Y(4)=y, Y(5)=psi
//    const double u = Y(0);
//    const double v = Y(1);
//    const double r = Y(2);
//    const double psi = Y(5);
//
//    const double L = ShipCfg.Geometry.Length;
//    const double d = ShipCfg.Geometry.Draft;
//
//    // 运动学
//    dY(3) = u * std::cos(psi) - v * std::sin(psi);
//    dY(4) = u * std::sin(psi) + v * std::cos(psi);
//    dY(5) = r;
//
//    const double U = std::hypot(u, v);
//    const double eps = 1.0e-10;
//    if (U < eps)
//        return dY;
//
//    const double q = 0.5 * rho * L * d * U * U;
//
//    // 1) 船体无量纲力
//    const HullForceND FH = HullForceFortranND(Y);
//
//    // 2) 螺旋桨推力、静水阻力、舵力
//    double wP = 0.0;
//    double Jp = 0.0;
//    double KT = 0.0;
//
//    const double XP = PropellerForce(Y, wP, Jp, KT);
//    const double XR = ResistanceForce(Y);
//    const RudderForceDim FR = RudderForce(Y, rudder_angle, Jp, KT);
//
//    // 3) 按 Fortran FCNM 的对角质量阵
//    const double M11 = (MassND.Mm1 + MassND.mx) * L / (U * U);
//    const double M22 = (MassND.Mm1 + MassND.my) * L / (U * U);
//    const double M66 = (MassND.Mm4 + MassND.Jz) * L * L / (U * U);
//
//    // 4) 连续方程（返回纯导数，不乘 dt）
//    dY(0) = (FH.X + (XP - XR + FR.X) / q) / M11;
//    dY(1) = (FH.Y + FR.Y / q) / M22;
//    dY(2) = (FH.N + FR.N / (q * L)) / M66;
//
//    return dY;
//}