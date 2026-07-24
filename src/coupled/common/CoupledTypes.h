#pragma once

#include <Eigen/Dense>
#include "../../wave/Encounter.h"
#include "../../wave/WaveComponent.h"
#include <vector>

// ===== Coupled-only types =====
// 全部加 Coupled 前缀，避免和纯耐波/纯操纵侧的全局类型重名。

struct CoupledSlowState3DOF
{
    double t = 0.0;

    double u = 0.0;
    double v = 0.0;
    double r = 0.0;

    double xe = 0.0;
    double ye = 0.0;
    double psi = 0.0;

    double delta = 0.0;

    double U = 0.0;
    double betaDrift = 0.0;
    double Fn = 0.0;

    /// Yaw angular acceleration dr/dt from MMG RHS at the integrated state (rad/s^2).
    /// Written to CSV to avoid post-processing derivative on coarse dtSlow.
    double rDotRadPerS2 = 0.0;

    /// Manoeuvring (slow-time) roll: heel angle [rad] and roll rate [rad/s].
    /// Only integrated when the 4-DOF MMG is enabled (maneuverDOF==4); otherwise
    /// stays 0 and every other field is bit-identical to the 3-DOF run. The total
    /// physical roll is this slow heel PLUS the wave-frequency roll from the
    /// seakeeping fast solver (linear superposition).
    double phi = 0.0;
    double p = 0.0;
};

struct CoupledExternalLoads3DOF
{
    double X = 0.0;
    double Y = 0.0;
    double N = 0.0;

    // 二阶力/漂移力预留接口
    double X2 = 0.0;
    double Y2 = 0.0;
    double N2 = 0.0;
    // 二阶平均横摇(横倾)矩。仅在 4-DOF 横摇里使用：迎/随浪≈0，横/斜浪给出
    // 稳态波浪致横倾。3-DOF 与现有路径不读它，默认 0 -> 行为不变。
    double K2 = 0.0;
    bool hasSecondOrder = false;
};

struct CoupledWaveLoadsWindow
{
    // 当前版本：快窗耐波不反推 MMG 低频载荷。
    // 这些量保留为调试/后续扩展接口。
    double X_mean = 0.0;
    double Y_mean = 0.0;
    double N_mean = 0.0;

    double X_lf = 0.0;
    double Y_lf = 0.0;
    double N_lf = 0.0;
};

struct CoupledEncounterState
{
    EncounterInfo kin;
    double betaAbs = 0.0;
    double betaRel = 0.0;
    double omega = 0.0;
    double waveAmp = 0.0;
    double phaseAtCG = 0.0;

    // Full incident-sea decomposition (a, ω, θ_abs, ε). Regular = 1 component
    // equal to the scalars above; irregular/cross = many. Providers that have
    // been migrated sum over this; legacy ones keep using the scalars, so a
    // regular-wave run is bit-identical until each provider opts in.
    std::vector<WaveComponent> comps;
};

struct CoupledFastState
{
    Eigen::VectorXd q;   // 当前窗口起点位移
    Eigen::VectorXd v;   // 当前窗口起点速度
    Eigen::VectorXd a;   // 当前窗口起点加速度

    std::vector<Eigen::VectorXd> vMemHist; // 辐射卷积历史
    std::vector<double> etaHist;           // 激励卷积历史

    bool initialized = false;
};

/// Optional kinematics passed to IWaveForceProvider::onWindowBegin so providers
/// can evaluate η_ref(t−τ) with the ship position at t−τ (not only at t).
struct CoupledFastWindowInfo
{
    const CoupledSlowState3DOF* slowEnd = nullptr;
    double t0 = 0.0;
    double t1 = 0.0;
};

struct CoupledWindowRequest
{
    CoupledSlowState3DOF slow0;
    CoupledSlowState3DOF slow1_pred;

    double t0 = 0.0;
    double dtFast = 0.0;
    double dtSlow = 0.0;
    int nFast = 0;

    CoupledFastState fast0;
};

struct CoupledWindowResult
{
    CoupledWaveLoadsWindow loads;
    Eigen::VectorXd etaMean;
    Eigen::VectorXd etaRms;

    Eigen::VectorXd tHist;
    Eigen::MatrixXd qHist; // nFast x DOF
    Eigen::MatrixXd vHist; // nFast x DOF

    // 新增：真正进入耐波运动方程的波浪力
    // forceHist: nFast x DOF，顺序与 skCfg.modes 一致
    Eigen::MatrixXd forceHist;

    // 新增：完整 6DOF 波浪力，固定顺序 F0,F1,F2,F3,F4,F5
    Eigen::MatrixXd force6Hist;

    // 新增：调试波浪相位、相对浪向、插值浪向等
    Eigen::VectorXd betaRelDegHist;       // 当前真实相对浪向，0~360 deg
    Eigen::VectorXd betaQueryDegHist;     // 查表用浪向，0~180 deg
    Eigen::VectorXd mirroredHist;         // 是否由 180~360 镜像而来
    Eigen::VectorXd thetaHist;            // 用于 cos(theta + phase) 的相位
    Eigen::VectorXd omegaEncounterHist;   // 遭遇频率
    Eigen::VectorXd waveAmpHist;          // 当前波幅
    Eigen::VectorXd etaCGHist;            // 当前船舶重心处波面升高

    CoupledFastState fastLast;
};


