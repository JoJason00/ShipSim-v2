#pragma once

#include "../common/CoupledTypes.h"
#include "../../const/Const.h"

#include <fstream>
#include <string>
#include <vector>

// Coupled output writer (slim).
//   每个工况输出 3 个 CSV：
//     <fileTag>_slow.csv  — 慢回路（操纵）状态时历
//     <fileTag>_ext.csv   — 加给慢回路的外部力 / 慢回路看到的二阶载荷
//     <fileTag>_fast.csv  — 快回路（耐波）广义位移时历 (q0,q1,q2 = heave/roll/pitch
//                            按 Seakeeping.modes 顺序排列)
//   轨迹列顺序：ye (横轴) 在 xe (纵轴) 之前，与可视化“Y 横、X 纵”的地图视角一致。
class CoupledWriter
{
public:
    explicit CoupledWriter(
        const std::string& casePath,
        double refLppForNondim = 175.0);
    ~CoupledWriter();

    void openCase(const std::string& folderTag, const std::string& fileTag);

    /// 每个耦合算例开始时设置：U0 初速、入射波 ω（仅用于在 slow 末尾留底，不再做大量派生列）。
    void setScenarioReference(
        double U0_ref_mps,
        double omega_incident_rad_s,
        double zeta_a_m,
        const std::vector<int>& seakeepingModes);

    void closeCase();

    void writeSlowState(const CoupledSlowState3DOF& s);
    void writeExternalLoads(double t, const CoupledExternalLoads3DOF& ext);
    /// `phiSlow` = manoeuvring slow heel [rad] at this window; written as the
    /// extra "total roll" columns (wave roll + slow heel), dim + nondim.
    void writeFastSummary(const CoupledWindowResult& win, double phiSlow = 0.0);

private:
    std::string casePath_;
    double refLppForNondim_{ 175.0 };

    double U0Ref_{ 1.0 };
    double omegaIncRef_{ 0.0 };
    double zetaARef_{ 1.0 };
    double timeNondimScale_{ 0.0 };

    std::vector<int> modes_;  // 耐波 modes 列顺序（MODE_HEAVE/ROLL/PITCH 等），供 fast 列头注释

    std::ofstream slowOfs_;
    std::ofstream extOfs_;
    std::ofstream fastOfs_;
    std::ofstream waveForceOfs_;   // 耐波 6 自由度波浪激励力时历 (<fileTag>_waveforce.csv)
};
