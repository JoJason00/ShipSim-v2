// mmg.h
#pragma once

#include "../config/CaseConfig.h"
#include <string>
#include "../wave/WaveBase.h"
#include <memory>
#include "Eigen/Dense"
#include <wave/RegularWave.h>
#include <vector>
#include <fstream>

using Vec6 = Eigen::Matrix<double, 6, 1>;
using Mat6 = Eigen::Matrix<double, 6, 6>;
using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix<double, 3, 3>;

struct CaseContext_mmg
{
    int i_wave{};

    std::shared_ptr<WaveBase> wave;
    std::shared_ptr<RegularWave> reg;

    double dt{};
    double Fn{};
    double U{};

    // 先保留旧接口，当前静水版本未实际使用
    Mat3 A;
    Eigen::PartialPivLU<Mat3> A_lu;
};

class Mmg
{
public:
    Mmg(const ShipConfig& Ship, std::string casePath, const MmgConfig& Mmg);
    void run();
    void runTurningCase(const CaseContext_mmg& ctx, const bool direction);
    void runZigzagCase(const CaseContext_mmg& ctx);

private:
    // 经验公式版本
    void compute_hullCoeff();
    void compute_propellerCoeff();
    void compute_rudderCoeff();
    void compute_addedCoeff();

    CaseContext_mmg buildCaseContext(const int i_wave, const double Fn);

private:
    struct State
    {
        double u = 0.0;   // surge [m/s]
        double v = 0.0;   // sway  [m/s]
        double r = 0.0;   // yaw rate [rad/s]
        double xe = 0.0;  // earth-fixed x [m]
        double ye = 0.0;  // earth-fixed y [m]
        double psi = 0.0; // heading [rad]
        double del = 0.0; // rudder angle [rad]
    };

    struct Deriv
    {
        double du = 0.0;
        double dv = 0.0;
        double dr = 0.0;
        double dxe = 0.0;
        double dye = 0.0;
        double dpsi = 0.0;
        double ddel = 0.0;
    };

    struct ForceBreakdown
    {
        double U = 0.0;
        double beta = 0.0;
        double r_nd = 0.0;

        double XH = 0.0;
        double YH = 0.0;
        double NH = 0.0;
        double Ru = 0.0;

        double wp = 0.0;
        double Jp = 0.0;
        double KT = 0.0;
        double XP = 0.0;

        double uR = 0.0;
        double vR = 0.0;
        double UR = 0.0;
        double aR = 0.0;
        double FN = 0.0;

        double XR = 0.0;
        double YR = 0.0;
        double NR = 0.0;
    };

private:
    void setupModelScaleS175();

    ForceBreakdown evalForces(double t, const State& s, const CaseContext_mmg& ctx) const;
    Deriv rhs(double t, const State& s, const CaseContext_mmg& ctx,
        double deltaCmdRad, double localDt,
        ForceBreakdown* fb = nullptr) const;
    State rk4Step(double t, const State& s, double dt,
        const CaseContext_mmg& ctx, double deltaCmdRad) const;

    void writeCsvHeader(std::ofstream& ofs) const;
    void writeCsvRow(std::ofstream& ofs, double t, const State& s,
        const ForceBreakdown& fb) const;

    static double deg2rad(double deg);
    static double rad2deg(double rad);
    static double clamp(double x, double lo, double hi);

private:
    std::string         filePath;
    ShipConfig          ShipCfg;

    std::vector<std::shared_ptr<WaveBase>> waves;
    MmgTimeConfig Time;
    hullCoeff Hull;
    propellerCoeff Propeller;
    rudderCoeff Rudder;
    addedCoeff Added;

    std::vector<double> Fns;

    TuringCfg Turning;
    ZigzagCfg Zigzag;

    bool turn_on_turningcase{ false };
    bool turn_on_zigzagcase{ false };

private:
    // =========================
    // 当前验证版：S175 模型尺度参数
    // =========================
    double rho_ = 1000.0;
    double g_ = 9.81;
    double pi_ = 3.14159265358979323846;

    // Geometry / mass
    double Lpp_ = 3.5000;
    double Lwl_ = 3.5640;
    double B_ = 0.5080;
    double d_ = 0.1900;
    double Cb_ = 0.5720;

    double m_ = 193.57;    // kg
    double xG_ = -0.0510;
    double kzz_ = 0.2690;
    double Izz_ = 0.0;

    // Propeller / rudder geometry
    double Dp_ = 0.1301;
    double HR_ = 0.1540;
    double BR_ = 0.0843;
    double AR_ = 0.0130;

    double lambdaR_ = 0.0;
    double eta_ = 0.0;

    // Positions
    double xP_nd_ = 0.0;
    double xR_nd_ = -0.50;
    double xH_nd_ = -0.48;
    double lR_nd_ = -1.0;

    // Hydrodynamic derivatives
    double Xbb_ = -0.0711;
    double Xbr_ = 0.0573;
    double Xrr_ = 0.0037;

    double Yb_ = 0.2137;
    double Yr_ = 0.0446;
    double Ybbb_ = 2.0080;
    double Ybbr_ = 0.3942;
    double Ybrr_ = 0.7461;
    double Yrrr_ = 0.0326;

    double Nb_ = 0.0710;
    double Nr_ = -0.0409;
    double Nbbb_ = -0.0275;
    double Nbbr_ = -0.7811;
    double Nbrr_ = -0.0287;
    double Nrrr_ = -0.0422;

    double mx_nd_ = 0.0044;
    double my_nd_ = 0.1299;
    double Jz_nd_ = 0.0077;
    double Xuu_ = 0.01563;

    double mx_ = 0.0;
    double my_ = 0.0;
    double Jz_ = 0.0;

    // Rudder coefficients
    double epsilon_ = 0.921;
    double kappa_ = 0.631;
    double tR_ = 0.29;
    double aH_ = 0.237;
    double gammaR_plus_ = 0.139;
    double gammaR_minus_ = 0.088;
    double deltaRateDeg_ = 12.0;   

    // Propeller coefficients
    double C1_ = -8.0;
    double wp0_ = 0.1684;
    double tp_ = 0.175;
    double J0_ = 0.2932;
    double J1_ = -0.1971;
    double J2_ = -0.0481;

    // Operating condition
    double Fn_ = 0.15;
    double U0_ = 0.879;       // 论文工况
    double np_rps_ = 10.05;   // 固定转速
    double dt_ = 0.01;
    double tMax_ = 300.0;
    double approachTime_ = 0.0;  
    double stopPsiDeg_ = 540.0;
    double maxDeltaDeg_ = 35.0;
};