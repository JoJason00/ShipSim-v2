#pragma once

#include <string>

struct CoupledTimeConfig
{
    double dtSlow = 0.1;
    double dtFast = 0.01;
    double totalTime = 120.0;
    int predictorCorrectorIters = 1;

    /// Soft-start ramp for the incident wave amplitude (seconds). Over the first
    /// waveRampTimeS the wave amplitude is faded in with a raised-cosine window
    /// 0.5·(1−cos(π·t/T)), so a near-undamped roll mode is not hit by a step at
    /// t=0. Applied at the single encounter source, so the fast-loop excitation
    /// (∝ amp), the second-order drift loads (∝ amp²) and the CG wave elevation
    /// all ramp consistently. 0 disables the ramp (legacy step behaviour). A
    /// good value is ≈ 2 roll natural periods.
    double waveRampTimeS = 0.0;
};

struct HydroRefreshConfig
{
    double betaTolDeg = 5.0;
    double omegaTolRatio = 0.03;
    double speedTolRatio = 0.03;

    bool enableExcitationInterpolation = true;
    bool enableRadiationInterpolation = true;

    /// When true, the Fn / wave-heading refresh tolerances are derived from the
    /// actual excitation-kernel library grid spacing at the current incident
    /// frequency instead of the fixed betaTolDeg / speedTolRatio above. A denser
    /// library then refreshes more often (and more accurately). Falls back to the
    /// fixed values for any axis with fewer than two distinct library nodes.
    bool adaptiveRefresh = false;

    /// Refresh once the state has moved more than this fraction of one library
    /// grid cell (only used when adaptiveRefresh is true). 0.5 ≈ half a cell.
    double adaptiveRefreshFactor = 0.5;
};

struct CoupledOutputConfig
{
    bool saveSlowStates = true;
    bool saveWindowLoads = true;
    bool saveFastSummary = true;
    bool saveExternalLoads = true;

    /// 若 >0 且小于实船 \p Lpp，在 *_slow.csv 末尾附加弗劳德换算到该「文献船模」尺度的
    /// u_mdl_mps, v_mdl_mps, U_mdl_mps, t_mdl_s（与实船 u,v,U,t 对应，U 为合速度）。
    /// 设为 0 可关闭附加列。
    double referenceModelLppM = 3.5;
};

struct CoupledPhysicsConfig
{
    bool enableFastSeakeeping = true;
    bool enableSecondOrderLoads = true;

    // "impulseKernel" — convolve precomputed Qlag with eta history (TDGF kernel)
    // "directFK"      — integrate incident pressure over panels at every step
    // "directFKdf"    — directFK + diffraction impulse (fkImpulse df columns)
    // "fkImpulse"     — fk+df impulse from cases/{case}/fkImpulse/
    // "fitted"        — legacy lookup table (FittedWaveForceWindowSolver)
    std::string waveForceMode = "impulseKernel";

    // Toggle for the radiation memory term in the Cummins equation.
    // When false, the stepper drops A_inf / B / C' / K(t) and reduces to
    // M_phys · q̈ + b_visc · q̇_roll + C_phys · q = F_exc.
    // Useful for bring-up before the per-Fn radiation kernel database is complete.
    bool enableRadiation = true;
};

struct CoupledConfig
{
    CoupledTimeConfig time;
    HydroRefreshConfig refresh;
    CoupledOutputConfig output;
    CoupledPhysicsConfig physics;
};
