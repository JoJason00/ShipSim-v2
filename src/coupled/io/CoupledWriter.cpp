#include "CoupledWriter.h"

#include "../../seakeeping/SeakeepingDOF.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <stdexcept>

namespace
{
    const char* modeShortName(int mode)
    {
        switch (mode)
        {
        case MODE_SURGE: return "surge";
        case MODE_SWAY:  return "sway";
        case MODE_HEAVE: return "heave";
        case MODE_ROLL:  return "roll";
        case MODE_PITCH: return "pitch";
        case MODE_YAW:   return "yaw";
        default:         return "mode";
        }
    }
}

CoupledWriter::CoupledWriter(
    const std::string& casePath,
    double refLppForNondim)
    : casePath_(casePath),
      refLppForNondim_(refLppForNondim > 0.0 ? refLppForNondim : 175.0)
{
}

CoupledWriter::~CoupledWriter()
{
    closeCase();
}

void CoupledWriter::setScenarioReference(
    double U0_ref_mps,
    double omega_incident_rad_s,
    double zeta_a_m,
    const std::vector<int>& seakeepingModes)
{
    U0Ref_ = std::max(U0_ref_mps, 1.0e-6);
    omegaIncRef_ = omega_incident_rad_s;
    zetaARef_ = std::max(zeta_a_m, 1.0e-12);
    timeNondimScale_ = U0Ref_ / refLppForNondim_;
    modes_ = seakeepingModes;
}

void CoupledWriter::closeCase()
{
    if (slowOfs_.is_open()) slowOfs_.close();
    if (extOfs_.is_open()) extOfs_.close();
    if (fastOfs_.is_open()) fastOfs_.close();
    if (waveForceOfs_.is_open()) waveForceOfs_.close();
}

void CoupledWriter::openCase(const std::string& folderTag, const std::string& fileTag)
{
    closeCase();

    const std::filesystem::path dir =
        std::filesystem::path(casePath_) / folderTag;

    std::filesystem::create_directories(dir);

    slowOfs_.open((dir / (fileTag + "_slow.csv")).string());
    extOfs_.open((dir / (fileTag + "_ext.csv")).string());
    fastOfs_.open((dir / (fileTag + "_fast.csv")).string());
    waveForceOfs_.open((dir / (fileTag + "_waveforce.csv")).string());

    if (!slowOfs_ || !extOfs_ || !fastOfs_ || !waveForceOfs_)
        throw std::runtime_error("CoupledWriter: failed to open output files.");

    // 耐波 6 自由度波浪激励力时历（有量纲：力 [N]，力矩 [N·m]）。
    // 列序与 force6Hist 一致：0 surge,1 sway,2 heave,3 roll,4 pitch,5 yaw。
    // Mx_roll(横摇力矩)=第4列、My_pitch(纵摇力矩)=第5列，供后处理对比"横浪横摇激励大/迎浪纵摇激励大"。
    waveForceOfs_ << "t,Fx_surge,Fy_sway,Fz_heave,Mx_roll,My_pitch,Mz_yaw\n";

    // Slow state CSV — Y (横轴) 在 X (纵轴) 之前。
    slowOfs_
        << "t,t_star,"
        << "ye,xe,ye_Lpp,xe_Lpp,"
        << "u,v,r,psi,delta,U,betaDrift,Fn,Fn_inst,"
        << "r_degps,rdot_rad_s2,"
        << "phi_roll_deg,phi_roll_rad,"   // slow heel (4-DOF MMG); 0 in 3-DOF
        // 艏向角缠绕到 [-180,180]°（便于读"什么浪向"）。注意：psi(弧度,累积/未缠绕)
        // 仍保留在前面——回转参数(进距/横距/战术直径)靠 psi 累积过 90/180°，不能用缠绕值。
        << "psi_deg"
        << "\n";

    extOfs_ << "t,X,Y,N,X2,Y2,N2,hasSecondOrder\n";

    // Fast (seakeeping) CSV — 左侧为有量纲量 (t [s], q [m] 或 [rad])，
    // 右侧补出对应的无量纲量 (t* = t·U0/Lpp, q* = q/ζa 或 q/(k·ζa))。
    fastOfs_ << "t";
    const std::size_t nModes = modes_.size();
    for (std::size_t i = 0; i < nModes; ++i)
        fastOfs_ << ",q" << i << "_" << modeShortName(modes_[i]);
    if (nModes == 0)
        fastOfs_ << ",q0,q1,q2";

    fastOfs_ << ",t_star";
    for (std::size_t i = 0; i < nModes; ++i)
        fastOfs_ << ",q" << i << "_star_" << modeShortName(modes_[i]);
    if (nModes == 0)
        fastOfs_ << ",q0_star,q1_star,q2_star";
    // Total roll = wave-frequency roll (seakeeping) + slow heel (manoeuvring),
    // dimensional [rad] and non-dimensional. 0 + wave-roll if 3-DOF (no heel).
    fastOfs_ << ",roll_total,roll_total_star";
    fastOfs_ << "\n";
}

void CoupledWriter::writeSlowState(const CoupledSlowState3DOF& s)
{
    if (!slowOfs_.is_open())
        return;

    const double Lpp = refLppForNondim_;
    const double Fn_inst = s.U / std::sqrt(G * Lpp);
    const double tStar = s.t * timeNondimScale_;

    slowOfs_ << std::fixed << std::setprecision(10)
        << s.t << ","
        << tStar << ","
        << s.ye << "," << s.xe << ","
        << s.ye / Lpp << "," << s.xe / Lpp << ","
        << s.u << "," << s.v << "," << s.r << ","
        << s.psi << "," << s.delta << "," << s.U << "," << s.betaDrift << ","
        << s.Fn << "," << Fn_inst << ","
        << s.r * (180.0 / PI) << ","
        << s.rDotRadPerS2 << ","
        << s.phi * (180.0 / PI) << "," << s.phi;

    // psi_deg: 艏向角缠绕到 [-180,180]°
    double psiDeg = std::fmod(s.psi * (180.0 / PI) + 180.0, 360.0);
    if (psiDeg < 0.0) psiDeg += 360.0;
    psiDeg -= 180.0;
    slowOfs_ << "," << psiDeg << "\n";
}

void CoupledWriter::writeExternalLoads(double t, const CoupledExternalLoads3DOF& ext)
{
    if (!extOfs_.is_open())
        return;

    extOfs_ << std::fixed << std::setprecision(10)
        << t << ","
        << ext.X << "," << ext.Y << "," << ext.N << ","
        << ext.X2 << "," << ext.Y2 << "," << ext.N2 << ","
        << (ext.hasSecondOrder ? 1 : 0)
        << "\n";
}

void CoupledWriter::writeFastSummary(const CoupledWindowResult& win, double phiSlow)
{
    if (!fastOfs_.is_open())
        return;

    const int N = static_cast<int>(win.tHist.size());
    const int nCols = static_cast<int>(win.qHist.cols());

    // Roll column index in qHist (modes order), so we can add the slow heel.
    int jRoll = -1;
    for (std::size_t j = 0; j < modes_.size(); ++j)
        if (modes_[j] == MODE_ROLL) { jRoll = static_cast<int>(j); break; }

    // 无量纲化尺度：
    //   q*  = q / ζa            (平动自由度: surge/sway/heave)
    //   q*  = q / (k·ζa)         (转动自由度: roll/pitch/yaw), k = ω²/g
    //   t*  = t · U0/Lpp
    // 没有波（ω≈0 或 ζa≈0）时退化为直接输出 q。
    const double k = (omegaIncRef_ * omegaIncRef_) / G;
    const double zetaA = zetaARef_;
    const bool waveValid = (omegaIncRef_ > 1.0e-9) && (zetaA > 1.0e-9);

    auto isRotational = [](int mode) {
        return mode == MODE_ROLL || mode == MODE_PITCH || mode == MODE_YAW;
    };

    for (int i = 0; i < N; ++i)
    {
        fastOfs_ << std::fixed << std::setprecision(10) << win.tHist(i);
        for (int j = 0; j < nCols; ++j)
        {
            if (isRotational(modes_[j]))
                fastOfs_ << "," << win.qHist(i, j)*180.0/3.1416;
            else
                fastOfs_ << "," << win.qHist(i, j);
        }

        fastOfs_ << "," << (win.tHist(i) * timeNondimScale_);
        for (int j = 0; j < nCols; ++j)
        {
            double qs = win.qHist(i, j);
            if (waveValid)
            {
                const int mode = (j < static_cast<int>(modes_.size())) ? modes_[j] : -1;
                const double scale = isRotational(mode) ? (k * zetaA) : zetaA;
                qs = win.qHist(i, j) / scale;
            }
            fastOfs_ << "," << qs;
        }

        // Total roll = wave roll (this row) + slow heel (constant over window).
        const double rollWave = (jRoll >= 0) ? win.qHist(i, jRoll) : 0.0;
        const double rollTotal = rollWave + phiSlow;
        const double rollTotalStar =
            waveValid ? (rollTotal / (k * zetaA)) : rollTotal;
        fastOfs_ << "," << rollTotal << "," << rollTotalStar;

        fastOfs_ << "\n";
    }

    // 6 自由度波浪激励力时历 -> 单独的 <fileTag>_waveforce.csv。
    // 仅当求解器填了 force6Hist（N×6）时写；Null/缺省求解器为空时自动跳过。
    if (waveForceOfs_.is_open()
        && win.force6Hist.rows() == N
        && win.force6Hist.cols() == 6)
    {
        for (int i = 0; i < N; ++i)
        {
            waveForceOfs_ << std::fixed << std::setprecision(10) << win.tHist(i);
            for (int c = 0; c < 6; ++c)
                waveForceOfs_ << "," << win.force6Hist(i, c);
            waveForceOfs_ << "\n";
        }
    }
}
