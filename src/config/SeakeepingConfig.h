#pragma once

#include <vector>
#include <string>
#include "wave/WaveBase.h"
#include "../seakeeping/RollDamping.h"

struct PanelSettings
{
    std::string NEType;
    int NE{};
};

struct TimeSettings
{
    double dt{};
    int TimeStep{};
    int PreStep{};
    int GreenStep{};
    double TimeCircle{};
    double PreCircle{};
    double GreenCircle{};
};

struct FreeRollDecayConfig
{
    bool enabled = false;
    double phi0_deg = 10.0;
    double phidot0_deg_s = 0.0;
    double dt = 0.02;
    double duration = 120.0;
};

enum class FKForceMethod
{
    PrescribedAmp,
    TDGFImpulse,    // FK + diffraction both via the impulse kernel
    DirectPressure  // FK by direct panel-pressure integration (body frame),
                    // diffraction via the df part of the impulse kernel.
                    // (Formerly also exposed as "DirectFKImpulseDf"; the two
                    // were numerically identical and have been merged.)
};

struct TDGFConfig
{
    // Exact:    ʵʱ���� Greenf::GreenFunctionCal
    // Interp:   ʹ������ GreenTable::eval
    // ODETable: ���߽ڵ�� + ���� RK4 �ƽ�
    std::string Method = "Exact";

    // Interp Ĭ���� green_table.bin
    // ODETable Ĭ���� tdgf_ode_table.bin
    std::string TablePath = "green_table.bin";

    // ODETable ���� RK4 С����
    double RK4Step = 0.001;
};

struct FKImpulseKernelConfig
{
    // Method:
    //   "Auto"   — try library first; if no matching file found, compute and save.
    //   "Direct" — always recompute and overwrite the library file.
    //   "Load"   — only read from library; error out when nothing matches.
    std::string Method = "Auto";

    // Include the ship name in the cached file name.
    bool IncludeShipName = true;

    // ----- Physics-driven (dt, tMot) auto-derivation, used when building kernels.
    //
    // Lu et al. (2024, J. Ship Mech.) — diffraction/FK impulse with 3D TDFGF + Beck-type
    // free-surface memory: following / stern-quartering needs a long τ convolution window;
    // head seas decay faster. Nondimensional dumps go to {case}/fkImpulse/fkImpulse_nd_*.csv
    // (τ_nd = τ·√(g/L)); half-width ≈ (tMot·dt)/√(L/g) once Te-based mem hits MemoryCapNd).
    //
    // Practical tuning (region-adaptive ON):
    //   • Following F1/F2: MemoryCapNdFollowing (+ optional MemoryCapNdF1) and MemoryCapNdF2
    //     set the |τ_nd| half-width cap (τ_nd = τ·√(g/L)); use larger windows in bands 1–2.
    //   • Following F3 (ω_inc ≥ ω_c2): MemoryCapNdF3 defaults to ~10 (τ_nd ∈ [−10,+10] half-span
    //     is usually enough); F1/F2 are not limited by this.
    //   • F2 (near-critical |ω_e|): raise MemoryCapNdF2 / MemoryEncounterPeriodsF2 first.
    //   • Head (H): MemoryCapNdHead ≈ 6 keeps |τ_nd| ≲ 6 (compact head seas); raise SamplesPerWavePeriodHead
    //     only if numerical noise appears near τ=0. Impulse dt is also capped at 0.015·√(L/g)
    //     so head-seas FK/diffraction kernels match the legacy fixed grid (coarser dt was a regression).
    //   • Time step: dt = min(T_inc,T_e)/Samples*, clamped to DtCapNd*√(L/g). Aim ≥24–32
    //     samples per min(T_inc,T_e) in following seas; coarser caps speed runs but can
    //     undersmooth the impulse near τ=0.
    //   • MaxFKImpulseRows: if 2*tMot+1 would exceed this, dt is enlarged (up to DtCap) so
    //     long following-sea windows stay tractable; 0 disables. Head seas skip this relaxation.
    //
    // dt = min(T_inc, T_e) / SamplesPerWavePeriod, clamped to [DtFloorNd, DtCapNd] × √(L/g).
    double SamplesPerWavePeriod = 20.0;
    double DtFloorNd            = 0.005;
    double DtCapNd              = 0.10;

    // tMot·dt = MemoryEncounterPeriods · T_e, clamped to [MemoryFloorNd, MemoryCapNd] × √(L/g).
    // In F2 (near critical encounter frequency) T_e blows up — the cap keeps the
    // step count from exploding; raise MemoryCapNd if you suspect the tail is
    // still active.
    double MemoryEncounterPeriods = 8.0;
    double MemoryFloorNd          = 5.0;
    double MemoryCapNd            = 30.0;

    // ----- Region-adaptive grid (Lu et al. 2024 following-seas bands; see WaveForceRegion.h)
    //
    // When true, head/beam/bow-quartering (region H) and following (F1/F2/F3) use the
    // *Head / *Following / F2* fields below so that: finer dt + shorter memory in head
    // seas, and longer √(L/g)-scaled memory caps in following seas (especially F2) where
    // a single global cap truncates the diffraction / FK impulse tail too aggressively.
    bool RegionAdaptiveKernelGrid = true;

    // Head seas: (tMot·dt)/√(L/g) ≤ MemoryCapNdHead — default ~6 matches τ_nd ∈ [−6, +6] half-span.
    double SamplesPerWavePeriodHead = 40.0;
    double DtCapNdHead             = 0.05;
    double MemoryEncounterPeriodsHead = 3.5;
    double MemoryFloorNdHead       = 2.0;
    double MemoryCapNdHead        = 6.0;

    // Following / F1 & F3: |τ_nd| half-width cap (τ_nd = τ·√(g/L)); 150 → test window ≈ ±150.
    double SamplesPerWavePeriodFollowing = 30.0;
    double DtCapNdFollowing         = 0.10;
    double MemoryEncounterPeriodsFollowing = 12.5;
    double MemoryFloorNdFollowing   = 6.0;
    double MemoryCapNdFollowing     = 150.0;

    // F1 (ω_inc < ω_c1): optional extra cap vs Following; 0 = use MemoryCapNdFollowing only.
    double MemoryCapNdF1            = 0.0;

    // F2 (ω_c1 ≤ ω_inc < ω_c2): |ω_e| small — extend encounter-based window and cap.
    double MemoryEncounterPeriodsF2 = 15.0;
    double MemoryCapNdF2            = 150.0;

    // F3 (ω_inc ≥ ω_c2) following: shorter impulse window in τ_nd (default ±10).
    double MemoryCapNdF3            = 10.0;

    // Lower bound on memory span as a multiple of T_inc in following seas (adds to Te-based
    // window) so slow encounter modulation is not resolved with an overly short τ window.
    double FollowingMemoryMinIncidentPeriods = 3.5;

    // If the impulse table would have more than MaxFKImpulseRows (2*tMot+1), increase dt
    // up to the active DtCap so tMot drops — balances long τ windows vs CPU. 0 = off.
    int MaxFKImpulseRows = 2200;

    // Adaptive truncation (reserved): trim trailing rows below `TailTolRel · peak`.
    double TailTolRel = 1.0e-3;
};

struct KradOverrideConfig
{
    bool enabled = false;
    std::string file;   // path relative to the case folder
};

// 非均匀 chi 时间网格（仅离线 K 计算用，运动方程仍是等距 dt_online）。
// 每一段 (stride, count) 含义：
//   - stride: 步长 = stride * dt_fine（必须是正整数，确保所有节点都落在 dt_fine 细网格上）
//   - count : 该段的区间数（节点数贡献 = count；首点 t=0 由网格生成器自动加）
//
// 总节点数 = 1 + Σ count。最后一个节点的时间 = Σ (stride·count) · dt_fine。
//
// 例（你之前确认的方案）：
//   dt_fine = 0.05
//   blocks = [ (1,100), (4,125), (20,170) ]
//   memory_cutoff_lag = 800
// → 396 chi 节点，从 0 到 200s；chi 卷积里 lag > 800 (= 40s) 的项被截断。
//
// 当 blocks 为空或 enabled=false 时，走老的等距 chi build 路径。
struct ChiTimeGridBlock
{
    int stride = 1;     // 步长（dt_fine 倍数），必须 >= 1
    int count  = 0;     // 区间数，必须 >= 1
};

struct ChiTimeGridConfig
{
    bool   enabled  = false;
    double dt_fine  = 0.0;
    std::vector<ChiTimeGridBlock> blocks;
    int    memory_cutoff_lag = 0;   // <= 0 表示不截断（保留全长卷积）
};

struct SeakeepingConfig
{
    std::vector<std::shared_ptr<WaveBase>> waves;

    std::vector<double> Fn;
    int DOF;
    std::vector<int> modes;

    PanelSettings Panel;
    TimeSettings Time;

    std::string Solver;

    RollDampingConfig RollDamping;
    FreeRollDecayConfig FreeRollDecay;

    std::string ResponseMethod = "LegacyAB";
    FKForceMethod FKMethod = FKForceMethod::PrescribedAmp;

    TDGFConfig GreenFunction;

    FKImpulseKernelConfig FKImpulseKernel;

    KradOverrideConfig KradOverride;

    // 非均匀 K 网格（默认 enabled=false，走老的等距路径，结果与之前一致）
    ChiTimeGridConfig ChiTimeGrid;
};