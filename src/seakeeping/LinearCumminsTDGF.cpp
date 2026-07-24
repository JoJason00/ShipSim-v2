#include "LinearCumminsTDGF.h"
#include "RadiationKernelTailSplit.h"

#include <functional>

#include "../const/Const.h"
#include "../tool/Fit.h"
#include "../io/Write.h"
#include "../wave/IrregularWave.h"
#include "../wave/CrossWave.h"
#include "../wave/Wave.h"
#include "Element.h"

#include <omp.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <sstream>
#include <filesystem>
#include <array>
#include <tool/ChiHistorySmooth.h>
#include "../tool/ImpulseKernelSG.h"
#include "../tool/ParallelGuard.h"
#include "CumminsTimeStepper.h"
#include "../seakeeping/WaveForceRegion.h"
#include "../io/FKImpulseKernelIO.h"

namespace
{
    constexpr bool kRadiationKernelSgSmooth = true;
    constexpr bool kRadiationKernelSgCompareCsv = true;

    // Derive (dt, tMot, tMot·dt) for the FK impulse kernel from physics
    // (T_inc, T_e) and the user-tunable knobs in FKImpulseKernelConfig.
    // With RegionAdaptiveKernelGrid, head seas use a finer dt / shorter memory span;
    // following seas (esp. F2) relax the √(L/g) memory cap so the τ tail is not cut
    // where Lu et al. (2024) impulse formulation still carries free-surface memory.
    // Output τ_nd = τ·√(g/L) in fkImpulse.csv: half-range ≈ (tMot·dt)/√(L/g) once capped.
    struct KernelGrid { double dt; int tMot; double memory; };

    KernelGrid deriveKernelGrid(const CaseContext& ctx,
        double shipLength,
        const FKImpulseKernelConfig& cfg)
    {
        const double tDim = std::sqrt(shipLength / G);
        const double tInc = 2.0 * PI / std::max(ctx.W, 1e-6);
        const double tEnc = 2.0 * PI / std::max(std::abs(ctx.we), 1e-6);
        const double tShort = std::min(tInc, tEnc);

        const WaveForceRegion region =
            wave_force_region::classify(ctx.U, ctx.dirRad, ctx.W);

        double samp = std::max(cfg.SamplesPerWavePeriod, 4.0);
        double dtCapNd = cfg.DtCapNd;
        double memPeriods = cfg.MemoryEncounterPeriods;
        double memFloorNdUse = cfg.MemoryFloorNd;
        double memCapNdUse = cfg.MemoryCapNd;

        if (cfg.RegionAdaptiveKernelGrid)
        {
            if (region == WaveForceRegion::Head)
            {
                samp = std::max(cfg.SamplesPerWavePeriodHead, 4.0);
                dtCapNd = cfg.DtCapNdHead;
                memPeriods = cfg.MemoryEncounterPeriodsHead;
                memFloorNdUse = cfg.MemoryFloorNdHead;
                memCapNdUse = cfg.MemoryCapNdHead;
            }
            else
            {
                samp = std::max(cfg.SamplesPerWavePeriodFollowing, 4.0);
                dtCapNd = cfg.DtCapNdFollowing;
                memPeriods = cfg.MemoryEncounterPeriodsFollowing;
                memFloorNdUse = cfg.MemoryFloorNdFollowing;
                memCapNdUse = cfg.MemoryCapNdFollowing;

                if (region == WaveForceRegion::F2)
                {
                    memPeriods = std::max(memPeriods, cfg.MemoryEncounterPeriodsF2);
                    memCapNdUse = std::max(memCapNdUse, cfg.MemoryCapNdF2);
                }
                else if (region == WaveForceRegion::F3)
                {
                    // Third following band: compact τ window (Lu-type impulse decays faster).
                    memCapNdUse = cfg.MemoryCapNdF3;
                }
                else if (region == WaveForceRegion::F1
                    && cfg.MemoryCapNdF1 > 0.0)
                {
                    memCapNdUse = std::max(memCapNdUse, cfg.MemoryCapNdF1);
                }
            }
        }

        double dt = tShort / samp;
        dt = std::max(dt, cfg.DtFloorNd * tDim);
        dt = std::min(dt, dtCapNd * tDim);

        // Head seas: older builds used a fixed nondim step dt_nd ≈ 0.015 on the impulse
        // grid. Purely physics-based sampling (SamplesPerWavePeriodHead / DtCapNdHead) can
        // land noticeably coarser than that, which changes the stored FK+diffraction kernel
        // and downstream convolution. Never use a coarser impulse dt than the legacy cap.
        if (region == WaveForceRegion::Head)
            dt = std::min(dt, 0.015 * tDim);

        double mem = memPeriods * tEnc;

        if (cfg.RegionAdaptiveKernelGrid
            && region != WaveForceRegion::Head
            && cfg.FollowingMemoryMinIncidentPeriods > 0.0)
        {
            mem = std::max(mem, cfg.FollowingMemoryMinIncidentPeriods * tInc);
        }

        mem = std::max(mem, memFloorNdUse * tDim);
        mem = std::min(mem, memCapNdUse * tDim);

        int tMot = std::max(1, static_cast<int>(std::round(mem / dt)));

        // Row-cap relaxation is for long following-sea windows. Head impulses stay compact
        // (MemoryCapNdHead); coarsening dt here was a regression vs the legacy fixed grid.
        if (cfg.MaxFKImpulseRows > 2 && region != WaveForceRegion::Head)
        {
            const int rowCap = cfg.MaxFKImpulseRows;
            for (int iter = 0; iter < 8 && 2 * tMot + 1 > rowCap; ++iter)
            {
                const int tMotTarget = std::max(1, (rowCap - 1) / 2);
                const double dtNeed = mem / static_cast<double>(tMotTarget);
                const double dtCapAbs = dtCapNd * tDim;
                const double dtNew = std::min(
                    std::max(dtNeed, cfg.DtFloorNd * tDim),
                    dtCapAbs);
                if (dtNew <= dt * (1.0 + 1e-12))
                    break;
                dt = dtNew;
                tMot = std::max(1, static_cast<int>(std::round(mem / dt)));
            }
        }

        return { dt, tMot, tMot * dt };
    }

    FKImpulseKernelIO::Params makeFKImpulseParams(
        const std::string& filePath,
        const ShipConfig& shipCfg,
        const CaseContext& ctx)
    {
        FKImpulseKernelIO::Params p;

        p.casePath = filePath;
        p.shipName = shipCfg.Name;

        p.Fn = ctx.Fn;
        p.dirRad = ctx.dirRad;
        p.dt = ctx.dt_const;
        p.tMot = ctx.tMot;

        // In following seas the kernel structure depends on which side of
        // the critical encounter frequencies (ω_c1, ω_c2) we are on, so
        // tag each saved file with its region; head seas always tag "H".
        p.region = wave_force_region::classify(ctx.U, ctx.dirRad, ctx.W);

        p.includeShipName = true;

        return p;
    }


    // Fixed-precision number with trailing zeros trimmed: 0.70 -> "0.7".
    // Used for Fn / ω_e in the run tag; the per-system descriptors come from
    // Wave::conditionTag (shared with the coupled output-folder naming).
    static std::string shortNum(double x, int prec = 2)
    {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(prec) << x;
        std::string s = ss.str();
        if (s.find('.') != std::string::npos)
        {
            s.erase(s.find_last_not_of('0') + 1);
            if (!s.empty() && s.back() == '.') s.pop_back();
        }
        return s;
    }


    static Eigen::VectorXd smooth5Binomial(
        const Eigen::VectorXd& y,
        int passes = 2)
    {
        const int N = static_cast<int>(y.size());

        if (N <= 0)
            return y;

        Eigen::VectorXd a = y;
        Eigen::VectorXd b = y;

        for (int pass = 0; pass < passes; ++pass)
        {
            b = a;

            if (N >= 5)
            {
                for (int i = 2; i <= N - 3; ++i)
                {
                    b(i) =
                        (a(i - 2)
                            + 4.0 * a(i - 1)
                            + 6.0 * a(i)
                            + 4.0 * a(i + 1)
                            + a(i + 2)) / 16.0;
                }
            }

            if (N >= 3)
            {
                b(1) =
                    (a(0) + 2.0 * a(1) + a(2)) / 4.0;

                b(N - 2) =
                    (a(N - 3) + 2.0 * a(N - 2) + a(N - 1)) / 4.0;
            }

            b(0) = a(0);
            b(N - 1) = a(N - 1);

            a = b;
        }

        return a;
    }

    static Eigen::VectorXd differentiateCentral(
        const Eigen::VectorXd& y,
        double dt)
    {
        const int N = static_cast<int>(y.size());
        Eigen::VectorXd dy = Eigen::VectorXd::Zero(N);

        if (N < 2 || dt <= 0.0)
            return dy;

        if (N == 2)
        {
            const double v = (y(1) - y(0)) / dt;
            dy(0) = v;
            dy(1) = v;
            return dy;
        }

        // 二阶单边端点
        dy(0) =
            (-3.0 * y(0) + 4.0 * y(1) - y(2)) / (2.0 * dt);

        dy(N - 1) =
            (3.0 * y(N - 1) - 4.0 * y(N - 2) + y(N - 3)) / (2.0 * dt);

        // 中间点：优先四阶中心差分
        for (int i = 1; i <= N - 2; ++i)
        {
            if (i >= 2 && i <= N - 3)
            {
                dy(i) =
                    (-y(i + 2)
                        + 8.0 * y(i + 1)
                        - 8.0 * y(i - 1)
                        + y(i - 2)) / (12.0 * dt);
            }
            else
            {
                dy(i) =
                    (y(i + 1) - y(i - 1)) / (2.0 * dt);
            }
        }
        return dy;
    }

    struct SimpleAmpPhase
    {
        double amp = 0.0;
        double phase = 0.0; // y = amp * cos(w t + phase)
    };

    static int findModeColInCfg(const SeakeepingConfig& cfg, int mode)
    {
        for (int i = 0; i < cfg.DOF; ++i)
        {
            if (cfg.modes[i] == mode)
                return i;
        }
        return -1;
    }

    static double wrap2Pi(double x)
    {
        while (x < 0.0) x += 2.0 * PI;
        while (x >= 2.0 * PI) x -= 2.0 * PI;
        return x;
    }

    static SimpleAmpPhase fitOneCosineWithMean(
        const Eigen::VectorXd& y,
        double dt,
        double w,
        int keepPeriods = 4)
    {
        SimpleAmpPhase out;

        const int N = static_cast<int>(y.size());
        if (N < 20 || std::abs(w) < 1.0e-12 || dt <= 0.0)
            return out;

        const int periodPoints =
            std::max(12, static_cast<int>(std::round(2.0 * PI / (std::abs(w) * dt))));

        const int keepPoints = keepPeriods * periodPoints;
        const int i0 = std::max(0, N - keepPoints);

        // y = a*sin(wt) + b*cos(wt) + c
        Eigen::Matrix3d ATA = Eigen::Matrix3d::Zero();
        Eigen::Vector3d ATy = Eigen::Vector3d::Zero();

        for (int i = i0; i < N; ++i)
        {
            const double t = i * dt;
            Eigen::Vector3d p;
            p << std::sin(w * t), std::cos(w * t), 1.0;

            ATA.noalias() += p * p.transpose();
            ATy.noalias() += p * y(i);
        }

        Eigen::Vector3d coef = ATA.fullPivLu().solve(ATy);

        const double aSin = coef(0);
        const double bCos = coef(1);

        out.amp = std::sqrt(aSin * aSin + bCos * bCos);

        // a*sin(wt)+b*cos(wt) = Amp*cos(wt+phase)
        // b = Amp*cos(phase), a = -Amp*sin(phase)
        out.phase = std::atan2(-aSin, bCos);
        out.phase = wrap2Pi(out.phase);

        return out;
    }

    static std::array<SimpleAmpPhase, 6> fitAmpPhase6Modes(
        const Eigen::MatrixXd& hist,
        const SeakeepingConfig& cfg,
        double dt,
        double fitOmega,
        int keepPeriods = 4)
    {
        std::array<SimpleAmpPhase, 6> fp;
        for (int mode = 0; mode < 6; ++mode)
        {
            const int col = findModeColInCfg(cfg, mode);
            if (col >= 0 && col < hist.cols())
            {
                fp[mode] = fitOneCosineWithMean(
                    hist.col(col), dt, fitOmega, keepPeriods);
            }
        }
        return fp;
    }

    static void writeAmpPhase6Modes(
        std::ofstream& out,
        const std::array<SimpleAmpPhase, 6>& fp,
        const CaseContext& ctx)
    {
        for (int mode = 0; mode < 6; ++mode)
        {
            double ampNd = 0.0;
            const double denom = ctx.a_scale[static_cast<std::size_t>(mode)];
            if (std::abs(denom) > 1.0e-30)
                ampNd = fp[mode].amp / denom;
            else
                ampNd = fp[mode].amp;
            out << "," << ampNd << "," << fp[mode].phase;
        }
    }

    static void writeFirstOrderWaveForceAmpPhaseTable(
        const std::string& filePath,
        const CaseContext& ctx,
        const SeakeepingConfig& cfg,
        double shipLengthM)
    {
        namespace fs = std::filesystem;

        fs::create_directories(fs::path(filePath) / "ExcitingForce");

        const std::string outFile =
            filePath + "/ExcitingForce/FirstOrderWaveForceAmpPhase.csv";

        const double fitOmega =
            (std::abs(ctx.we) > 1.0e-12) ? std::abs(ctx.we) : std::abs(ctx.W);

        // sqrt(L/nada): nada = incident wavelength λ (deep water), same ω–λ as CaseLoader
        // (ω = sqrt(2πg/λ) ⇒ λ = 2πg/ω²). Lu et al. use sqrt(L/λ) on the horizontal axis.
        double sqrt_L_over_nada = 0.0;
        const double Wabs = std::abs(ctx.W);
        if (shipLengthM > 0.0 && Wabs > 1.0e-12)
        {
            const double lambda = 2.0 * PI * G / (Wabs * Wabs);
            if (lambda > 1.0e-12)
                sqrt_L_over_nada = std::sqrt(shipLengthM / lambda);
        }

        const std::array<SimpleAmpPhase, 6> fpTotal =
            fitAmpPhase6Modes(ctx.ExcitingForce, cfg, ctx.dt, fitOmega, 4);
        const std::array<SimpleAmpPhase, 6> fpFK =
            fitAmpPhase6Modes(ctx.FKForceHist, cfg, ctx.dt, fitOmega, 4);
        const std::array<SimpleAmpPhase, 6> fpDF =
            fitAmpPhase6Modes(ctx.DiffForceHist, cfg, ctx.dt, fitOmega, 4);

        const bool needHeader = !fs::exists(outFile) || fs::file_size(outFile) == 0;

        std::ofstream out(outFile, std::ios::app);
        if (!out.is_open())
            throw std::runtime_error("cannot create " + outFile);

        if (needHeader)
        {
            out << "Fn,betaDeg,omegaIncident,omegaEncounter,waveAmp,\"sqrt(L/nada)\"";
            for (int mode = 0; mode < 6; ++mode)
                out << ",A" << mode << "_nd,P" << mode;
            for (int mode = 0; mode < 6; ++mode)
                out << ",A" << mode << "_FK_nd,P" << mode << "_FK";
            for (int mode = 0; mode < 6; ++mode)
                out << ",A" << mode << "_DF_nd,P" << mode << "_DF";
            out << "\n";
        }

        out << std::fixed << std::setprecision(12)
            << ctx.Fn << ","
            << ctx.dirRad * 180.0 / PI << ","
            << ctx.W << ","
            << ctx.we << ","
            << ctx.Amp << ","
            << sqrt_L_over_nada;

        writeAmpPhase6Modes(out, fpTotal, ctx);
        writeAmpPhase6Modes(out, fpFK, ctx);
        writeAmpPhase6Modes(out, fpDF, ctx);

        out << "\n";

        std::cout << "[FirstOrderFit] write " << outFile << "\n";
        std::cout << "  Fn=" << ctx.Fn
            << " betaDeg=" << ctx.dirRad * 180.0 / PI
            << " omegaIncident=" << ctx.W
            << " omegaEncounter=" << ctx.we << "\n";

        for (int mode = 0; mode < 6; ++mode)
        {
            if (findModeColInCfg(cfg, mode) < 0)
                continue;
            const double denom = ctx.a_scale[static_cast<std::size_t>(mode)];
            const auto ampNd = [&](const SimpleAmpPhase& p) -> double
                {
                    return std::abs(denom) > 1.0e-30 ? p.amp / denom : p.amp;
                };
            std::cout << "  mode " << mode
                << " total Amp_nd=" << ampNd(fpTotal[mode])
                << " FK Amp_nd=" << ampNd(fpFK[mode])
                << " DF Amp_nd=" << ampNd(fpDF[mode])
                << " (Ptot=" << fpTotal[mode].phase
                << " PFK=" << fpFK[mode].phase
                << " PDF=" << fpDF[mode].phase << " rad)\n";
        }
    }
}


LinearCumminsTDGF::LinearCumminsTDGF(
    const ShipConfig& shipCfg,
    std::string casePath,
    const SeakeepingConfig& seakeepingCfg)
    : filePath(std::move(casePath)),
    ShipCfg(shipCfg),
    SeakeepingCfg(seakeepingCfg),
    ak(seakeepingCfg.Panel.NE),
    rVn(seakeepingCfg.Panel.NE),
    Solver(seakeepingCfg.Solver)
{
    cfgTimeDtInput_ = seakeepingCfg.Time.dt;
    cfgGreenStepInput_ = seakeepingCfg.Time.GreenStep;
}

void LinearCumminsTDGF::assertRigidModesOnly() const
{
    for (int m : SeakeepingCfg.modes)
    {
        if (m < 0 || m > 5)
        {
            throw std::runtime_error(
                "LinearCumminsTDGF: this implementation follows the paper rigid-body formulas; "
                "modes must be a subset of [0..5].");
        }
    }
}

void LinearCumminsTDGF::run()
{
    assertRigidModesOnly();

    processElement();
    hs = SeakeepingDOF::hydrostaticsFromWaterline(ShipCfg, SeakeepingCfg, *element);
    //initGreenTable();

    //根据工况选择 Exact / Interp / ODETable
    initTDGFProvider();

    if (SeakeepingCfg.FreeRollDecay.enabled)
    {
        runFreeRollDecay(filePath);
        return;
    }

    rollVisc_ = RollDampingBuilder::build(
        filePath,
        SeakeepingCfg.RollDamping,
        ShipCfg.Mass.Mass,
        ShipCfg.Mass.GM);

    std::cout << "[LinearCumminsTDGF] roll damping:\n"
        << "  B44_lin  = " << rollVisc_.B44_lin << "\n"
        << "  B44_quad = " << rollVisc_.B44_quad << "\n"
        << "  B44_cube = " << rollVisc_.B44_cube << "\n";

    solve();
}

void LinearCumminsTDGF::processElement()
{
    const std::string file = filePath + ShipCfg.Name + ".element";

    auto elementData = CaseLoader::loadelement(
        file, SeakeepingCfg.Panel.NEType, SeakeepingCfg.Panel.NE);

    element = std::make_shared<Element>(
        SeakeepingCfg.Solver, SeakeepingCfg.Panel.NE, std::move(elementData));

    Eigen::Vector3d cg;
    cg << ShipCfg.Mass.CG.at(0), ShipCfg.Mass.CG.at(1), ShipCfg.Mass.CG.at(2);

    element->Geometry(cg);
    element->RankineSource2();

    const int nWL = element->n_WL;
    N0sq.resize(nWL);
    PotL_idx.resize(nWL);
    wlTx_.resize(nWL);
    wlTy_.resize(nWL);
    wlN1_.resize(nWL);
    wlN2_.resize(nWL);
    wlSegLen_.resize(nWL);
    wlDn1n2Dl_.resize(nWL);
    wlOrder_.assign(nWL, 0);

    struct WLSeg
    {
        int orig = 0;
        int panel = 0;
        double x0 = 0.0, y0 = 0.0;
        double x1 = 0.0, y1 = 0.0;
        double mx = 0.0, my = 0.0;
        double n1 = 0.0, n2 = 0.0;
        double L = 1.0;
        double ang = 0.0;
    };

    std::vector<WLSeg> segs(static_cast<std::size_t>(nWL));
    double cx = 0.0, cy = 0.0;

    for (int k = 0; k < nWL; ++k)
    {
        WLSeg s;
        s.orig = k;
        s.panel = element->PotL[k];
        s.x0 = element->xpl(k, 0);
        s.y0 = element->ypl(k, 0);
        s.x1 = element->xpl(k, 1);
        s.y1 = element->ypl(k, 1);
        s.mx = 0.5 * (s.x0 + s.x1);
        s.my = 0.5 * (s.y0 + s.y1);
        s.n1 = -element->Nvec(s.panel, 0);
        s.n2 = -element->Nvec(s.panel, 1);
        const double dx = s.x1 - s.x0;
        const double dy = s.y1 - s.y0;
        s.L = std::sqrt(dx * dx + dy * dy);
        cx += s.mx;
        cy += s.my;
        segs[static_cast<std::size_t>(k)] = s;
    }


    //调整水线段顺序，使其按极角排序，逆时针方向
    //if (nWL > 0)
    //{
    //    const double TWO_PI = 2.0 * PI;
    //    const double EPS = 1e-12;

    //    auto wrapTo2Pi = [&](double a) -> double
    //    {
    //        while (a < 0.0) a += TWO_PI;
    //        while (a >= TWO_PI) a -= TWO_PI;
    //        return a;
    //    };

    //    auto flipSeg = [&](WLSeg& s)
    //    {
    //        std::swap(s.x0, s.x1);
    //        std::swap(s.y0, s.y1);
    //    };

    //    cx /= static_cast<double>(nWL);
    //    cy /= static_cast<double>(nWL);

    //    for (auto& s : segs)
    //    {
    //        // 中点
    //        s.mx = 0.5 * (s.x0 + s.x1);
    //        s.my = 0.5 * (s.y0 + s.y1);

    //        // 排序角：只看中点极角
    //        s.ang = wrapTo2Pi(std::atan2(s.my - cy, s.mx - cx));

    //        // 方向：切向量相对中心必须是逆时针
    //        const double rx = s.mx - cx;
    //        const double ry = s.my - cy;
    //        const double tx = s.x1 - s.x0;
    //        const double ty = s.y1 - s.y0;

    //        const double cross = rx * ty - ry * tx;

    //        if (cross < -EPS)
    //        {
    //            flipSeg(s);
    //        }
    //    }

    //    std::sort(segs.begin(), segs.end(), [](const WLSeg& a, const WLSeg& b)
    //    {
    //        return a.ang < b.ang;
    //    });

    //    // 可选：把序列起点旋到最右上，便于看结果
    //    int start = 0;
    //    for (int i = 1; i < nWL; ++i)
    //    {
    //        const auto& a = segs[static_cast<std::size_t>(i)];
    //        const auto& b = segs[static_cast<std::size_t>(start)];

    //        if (a.mx > b.mx + EPS ||
    //            (std::fabs(a.mx - b.mx) <= EPS && a.my > b.my))
    //        {
    //            start = i;
    //        }
    //    }
    //    std::rotate(segs.begin(), segs.begin() + start, segs.end());
    //}



    if (nWL > 0)
    {
        const double TWO_PI = 2.0 * PI;
        const double EPS = 1e-12;

        auto wrapTo2Pi = [&](double a) -> double
            {
                while (a < 0.0) a += TWO_PI;
                while (a >= TWO_PI) a -= TWO_PI;
                return a;
            };

        auto flipSeg = [&](WLSeg& s)
            {
                std::swap(s.x0, s.x1);
                std::swap(s.y0, s.y1);
            };

        auto segHand = [&](const WLSeg& s) -> double
            {
                if (s.L <= 1e-14) return 0.0;
                const double tx = (s.x1 - s.x0) / s.L;
                const double ty = (s.y1 - s.y0) / s.L;
                // (-ty, tx) · (n1, n2)
                return (-ty) * s.n1 + tx * s.n2;
            };

        cx /= static_cast<double>(nWL);
        cy /= static_cast<double>(nWL);

        for (auto& s : segs)
        {
            s.mx = 0.5 * (s.x0 + s.x1);
            s.my = 0.5 * (s.y0 + s.y1);
            s.ang = wrapTo2Pi(std::atan2(s.my - cy, s.mx - cx));

            const double rx = s.mx - cx;
            const double ry = s.my - cy;
            const double tx = s.x1 - s.x0;
            const double ty = s.y1 - s.y0;
            const double cross = rx * ty - ry * tx;

            if (cross < -EPS)
            {
                flipSeg(s);
            }
        }

        std::sort(segs.begin(), segs.end(), [](const WLSeg& a, const WLSeg& b)
            {
                return a.ang < b.ang;
            });

        int start = 0;
        for (int i = 1; i < nWL; ++i)
        {
            const auto& a = segs[static_cast<std::size_t>(i)];
            const auto& b = segs[static_cast<std::size_t>(start)];

            if (a.mx > b.mx + EPS ||
                (std::fabs(a.mx - b.mx) <= EPS && a.my > b.my))
            {
                start = i;
            }
        }
        std::rotate(segs.begin(), segs.begin() + start, segs.end());

        // -------- 按首尾连接关系再修正每段方向 --------
        for (int k = 1; k < nWL; ++k)
        {
            auto& prev = segs[static_cast<std::size_t>(k - 1)];
            auto& cur = segs[static_cast<std::size_t>(k)];

            const double d_keep = std::hypot(cur.x0 - prev.x1, cur.y0 - prev.y1);
            const double d_flip = std::hypot(cur.x1 - prev.x1, cur.y1 - prev.y1);

            if (d_flip < d_keep)
            {
                flipSeg(cur);
            }
        }

        // 闭合处也检查一下最后一段
        {
            const auto& first = segs.front();
            auto& last = segs.back();

            const double d_keep = std::hypot(last.x1 - first.x0, last.y1 - first.y0);
            const double d_flip = std::hypot(last.x0 - first.x0, last.y0 - first.y0);

            if (d_flip < d_keep)
            {
                flipSeg(last);
            }
        }

        // -------- 统一整条水线的手性 --------
        double refHand = 0.0;
        bool foundRef = false;
        for (const auto& s : segs)
        {
            const double h = segHand(s);
            if (std::fabs(h) > 1e-10)
            {
                refHand = h;
                foundRef = true;
                break;
            }
        }

        // 如果当前整条链的手性是反的，就整体 reverse + flip
        if (foundRef && refHand < 0.0)
        {
            std::reverse(segs.begin(), segs.end());
            for (auto& s : segs)
            {
                flipSeg(s);
            }
        }

        // -------- 调试输出 --------
        const double meanL = [&]()
            {
                double sumL = 0.0;
                for (const auto& s : segs) sumL += s.L;
                return (nWL > 0) ? sumL / static_cast<double>(nWL) : 1.0;
            }();
        const double connTol = std::max(1e-8, 1e-3 * meanL);

        for (int k = 0; k < nWL; ++k)
        {
            const auto& cur = segs[static_cast<std::size_t>(k)];
            const auto& nxt = segs[static_cast<std::size_t>((k + 1) % nWL)];

            const double gap = std::hypot(cur.x1 - nxt.x0, cur.y1 - nxt.y0);
            const double hand = segHand(cur);

            if (gap > connTol)
            {
                std::cout << "[WL bad connect] k=" << k
                    << " gap=" << gap
                    << "  (" << cur.x0 << "," << cur.y0 << ")->(" << cur.x1 << "," << cur.y1 << ")"
                    << "  next=(" << nxt.x0 << "," << nxt.y0 << ")->(" << nxt.x1 << "," << nxt.y1 << ")\n";
            }

            if (foundRef && hand * refHand < -1e-10)
            {
                std::cout << "[WL bad hand] k=" << k
                    << " hand=" << hand
                    << " ref=" << refHand << "\n";
            }
        }


        // 把当前逆时针链改成顺时针
        //std::reverse(segs.begin(), segs.end());
        //for (auto& s : segs)
        //{
        //    flipSeg(s);
        //}

        ////  如果 GreenCal_WL 内部按 element->xpl/ypl 的原始端点方向积分，
        ////    这里把端点方向也同步写回去
        //for (const auto& s : segs)
        //{
        //    element->xpl(s.orig, 0) = s.x0;
        //    element->ypl(s.orig, 0) = s.y0;
        //    element->xpl(s.orig, 1) = s.x1;
        //    element->ypl(s.orig, 1) = s.y1;
        //}
    }


    for (int k = 0; k < nWL; ++k)
    {
        const auto& s = segs[static_cast<std::size_t>(k)];

        //std::cout << "line segment\t" << k << ":\tx0:\t" << s.x0 << "\ty0:\t " << s.y0
        //    << " \tx1:\t" << s.x1 << " \ty1:\t" << s.y1
        //    << "\tn1:\t " << s.n1 << " \tn2:\t" << s.n2 << "\tPotL_idx:\t"<< s.panel <<"\n";

        wlOrder_[static_cast<std::size_t>(k)] = s.orig;
        PotL_idx(k) = s.panel;
        wlN1_(k) = s.n1;
        wlN2_(k) = s.n2;
        N0sq(k) = static_cast<GScalar>(s.n1 * s.n1);
        wlSegLen_(k) = s.L;

        const double dx = s.x1 - s.x0;
        const double dy = s.y1 - s.y0;
        if (s.L > 1e-14)
        {
            wlTx_(k) = dx / s.L;
            wlTy_(k) = dy / s.L;
        }
        else
        {
            wlTx_(k) = 1.0;
            wlTy_(k) = 0.0;
        }
    }

    //for (int k = 0; k < nWL; ++k)
    //{
    //    const int km1 = (k - 1 + nWL) % nWL;
    //    const int kp1 = (k + 1) % nWL;

    //    const double fkm1 = wlN1_(km1) * wlN2_(km1);
    //    const double fkp1 = wlN1_(kp1) * wlN2_(kp1);

    //    //const double ds = 0.5 * (wlSegLen_(km1) + wlSegLen_(kp1));
    //    //wlDn1n2Dl_(k) = (fkp1 - fkm1) / std::max(1e-12, ds);

    //    const double ds =
    //        0.5 * wlSegLen_(km1) + wlSegLen_(k) + 0.5 * wlSegLen_(kp1);

    //    wlDn1n2Dl_(k) = (fkp1 - fkm1) / std::max(1e-12, ds);
    //}


    wlDn1n2Dl_.setZero();

    std::vector<int> cornerAfter(static_cast<std::size_t>(nWL), 0);
    std::vector<double> tanJump(static_cast<std::size_t>(nWL), 0.0);
    std::vector<double> norJump(static_cast<std::size_t>(nWL), 0.0);
    std::vector<double> fRaw(static_cast<std::size_t>(nWL), 0.0);
    std::vector<double> fUse(static_cast<std::size_t>(nWL), 0.0);

    auto prevIdx = [&](int k) -> int
        {
            return (k - 1 + nWL) % nWL;
        };

    auto nextIdx = [&](int k) -> int
        {
            return (k + 1) % nWL;
        };

    auto clampUnit = [](double x) -> double
        {
            return std::max(-1.0, std::min(1.0, x));
        };

    auto angleBetween2D = [&](double ax, double ay, double bx, double by) -> double
        {
            const double na = std::hypot(ax, ay);
            const double nb = std::hypot(bx, by);
            if (na <= 1e-14 || nb <= 1e-14) return 0.0;
            const double c = clampUnit((ax * bx + ay * by) / (na * nb));
            return std::acos(c);
        };

    auto medianPositive = [](std::vector<double> v) -> double
        {
            v.erase(std::remove_if(v.begin(), v.end(),
                [](double x) { return x <= 1e-12; }),
                v.end());

            if (v.empty()) return 0.0;

            std::sort(v.begin(), v.end());
            const std::size_t n = v.size();
            if (n % 2 == 1) return v[n / 2];
            return 0.5 * (v[n / 2 - 1] + v[n / 2]);
        };

    // --------  计算相邻线段的切向/法向跳变 --------
    for (int k = 0; k < nWL; ++k)
    {
        const int kp1 = nextIdx(k);

        tanJump[static_cast<std::size_t>(k)] =
            angleBetween2D(
                wlTx_(k), wlTy_(k),
                wlTx_(kp1), wlTy_(kp1));

        norJump[static_cast<std::size_t>(k)] =
            angleBetween2D(
                wlN1_(k), wlN2_(k),
                wlN1_(kp1), wlN2_(kp1));

        fRaw[static_cast<std::size_t>(k)] = wlN1_(k) * wlN2_(k);
    }

    // -------- 自适应阈值：识别明显折角 --------
    const double deg = PI / 180.0;
    const double medTan = medianPositive(tanJump);
    const double medNor = medianPositive(norJump);

    // 下限防止把普通平滑曲率误判；上限防止阈值过大
    const double thrTan = std::min(60.0 * deg, std::max(25.0 * deg, 4.0 * medTan));
    const double thrNor = std::min(60.0 * deg, std::max(20.0 * deg, 4.0 * medNor));

    int nCorner = 0;
    for (int k = 0; k < nWL; ++k)
    {
        if (tanJump[static_cast<std::size_t>(k)] > thrTan ||
            norJump[static_cast<std::size_t>(k)] > thrNor)
        {
            cornerAfter[static_cast<std::size_t>(k)] = 1; // 折角位于 k 和 k+1 之间
            ++nCorner;
        }
    }

    std::cout << "[WL corner detect] count=" << nCorner
        << "  thrTanDeg=" << thrTan / deg
        << "  thrNorDeg=" << thrNor / deg << "\n";

    for (int k = 0; k < nWL; ++k)
    {
        if (cornerAfter[static_cast<std::size_t>(k)] != 0)
        {
            const int kp1 = nextIdx(k);
            std::cout << "[WL corner] between seg " << k << " and " << kp1
                << "  tanJumpDeg=" << tanJump[static_cast<std::size_t>(k)] / deg
                << "  norJumpDeg=" << norJump[static_cast<std::size_t>(k)] / deg
                << "\n";
        }
    }

    // 邻接中心点之间的弧长距离：segment-center -> next-segment-center
    auto centerStepRight = [&](int k) -> double
        {
            const int kp1 = nextIdx(k);
            return 0.5 * wlSegLen_(k) + 0.5 * wlSegLen_(kp1);
        };

    // 是否允许跨界走到右边/左边（不能跨折角）
    auto canStepRight = [&](int k) -> bool
        {
            // k 与 k+1 之间若是折角，则不能跨
            return cornerAfter[static_cast<std::size_t>(k)] == 0;
        };

    auto canStepLeft = [&](int k) -> bool
        {
            // k-1 与 k 之间若是折角，则不能跨
            return cornerAfter[static_cast<std::size_t>(prevIdx(k))] == 0;
        };

    // -------- 先对 f = n1*n2 做轻微局部平滑，但绝不跨折角 --------
    // 权重: self=4, first neighbors=2, second neighbors=1
    for (int k = 0; k < nWL; ++k)
    {
        double sum = 4.0 * fRaw[static_cast<std::size_t>(k)];
        double wsum = 4.0;

        if (canStepLeft(k))
        {
            const int km1 = prevIdx(k);
            sum += 2.0 * fRaw[static_cast<std::size_t>(km1)];
            wsum += 2.0;

            if (canStepLeft(km1))
            {
                const int km2 = prevIdx(km1);
                sum += 1.0 * fRaw[static_cast<std::size_t>(km2)];
                wsum += 1.0;
            }
        }

        if (canStepRight(k))
        {
            const int kp1 = nextIdx(k);
            sum += 2.0 * fRaw[static_cast<std::size_t>(kp1)];
            wsum += 2.0;

            if (canStepRight(kp1))
            {
                const int kp2 = nextIdx(kp1);
                sum += 1.0 * fRaw[static_cast<std::size_t>(kp2)];
                wsum += 1.0;
            }
        }

        fUse[static_cast<std::size_t>(k)] = sum / std::max(1e-12, wsum);
    }

    // --------  在同一光滑链内，用局部二次最小二乘重构导数 --------
    auto fitDerivativeLocal = [&](int k) -> double
        {
            struct Sample
            {
                double s; // 相对第 k 段中心的弧长坐标
                double y; // fUse
                double w; // 拟合权重
            };

            std::vector<Sample> samples;
            samples.reserve(5);

            samples.push_back({ 0.0, fUse[static_cast<std::size_t>(k)], 1.0 });

            int left = k;
            int right = k;
            double sLeft = 0.0;
            double sRight = 0.0;

            while (static_cast<int>(samples.size()) < 5)
            {
                bool added = false;

                if (canStepLeft(left))
                {
                    const int km1 = prevIdx(left);
                    sLeft -= centerStepRight(km1);
                    left = km1;

                    const double w = 1.0 / std::max(1e-12, std::abs(sLeft) + 0.25 * wlSegLen_(k));
                    samples.push_back({ sLeft, fUse[static_cast<std::size_t>(left)], w });
                    added = true;

                    if (static_cast<int>(samples.size()) >= 5)
                        break;
                }

                if (canStepRight(right))
                {
                    sRight += centerStepRight(right);
                    right = nextIdx(right);

                    const double w = 1.0 / std::max(1e-12, std::abs(sRight) + 0.25 * wlSegLen_(k));
                    samples.push_back({ sRight, fUse[static_cast<std::size_t>(right)], w });
                    added = true;
                }

                if (!added)
                    break;
            }

            // 少于 3 点，不做拟合，导数置零
            if (static_cast<int>(samples.size()) < 3)
                return 0.0;

            // 拟合 y(s) = a0 + a1 s + a2 s^2
            // 则导数 y'(0) = a1
            Eigen::Matrix3d ATA = Eigen::Matrix3d::Zero();
            Eigen::Vector3d ATy = Eigen::Vector3d::Zero();

            for (const auto& q : samples)
            {
                const Eigen::Vector3d r(1.0, q.s, q.s * q.s);
                ATA.noalias() += q.w * (r * r.transpose());
                ATy.noalias() += q.w * r * q.y;
            }

            Eigen::Vector3d a = ATA.ldlt().solve(ATy);

            if (!a.allFinite())
            {
                // fallback：尽量用最近的单边或双边差分
                const bool hasL = canStepLeft(k);
                const bool hasR = canStepRight(k);

                if (hasL && hasR)
                {
                    const int km1 = prevIdx(k);
                    const int kp1 = nextIdx(k);
                    const double ds =
                        0.5 * wlSegLen_(km1) + wlSegLen_(k) + 0.5 * wlSegLen_(kp1);

                    return (fUse[static_cast<std::size_t>(kp1)] -
                        fUse[static_cast<std::size_t>(km1)]) / std::max(1e-12, ds);
                }
                else if (hasR)
                {
                    const int kp1 = nextIdx(k);
                    const double ds = centerStepRight(k);
                    return (fUse[static_cast<std::size_t>(kp1)] -
                        fUse[static_cast<std::size_t>(k)]) / std::max(1e-12, ds);
                }
                else if (hasL)
                {
                    const int km1 = prevIdx(k);
                    const double ds = centerStepRight(km1);
                    return (fUse[static_cast<std::size_t>(k)] -
                        fUse[static_cast<std::size_t>(km1)]) / std::max(1e-12, ds);
                }
                return 0.0;
            }

            return a(1);
        };

    for (int k = 0; k < nWL; ++k)
    {
        wlDn1n2Dl_(k) = fitDerivativeLocal(k);
    }


    // -------- 最后，凡是折角相邻的段，导数都置零 --------
    auto touchesCorner = [&](int k) -> bool
        {
            // 左边界：k-1 与 k 之间
            const bool leftCorner =
                cornerAfter[static_cast<std::size_t>(prevIdx(k))] != 0;

            // 右边界：k 与 k+1 之间
            const bool rightCorner =
                cornerAfter[static_cast<std::size_t>(k)] != 0;

            return leftCorner || rightCorner;
        };

    for (int k = 0; k < nWL; ++k)
    {
        // 只要该段相邻两侧任一处是折角，
        // 就认为这里的 d(n1*n2)/dl 不可可靠定义，直接置零
        if (touchesCorner(k))
        {
            wlDn1n2Dl_(k) = 0.0;
            std::cout << "[WL dn1n2/dl zeroed at corner-adjacent seg] k=" << k << "\n";
            continue;
        }

        //std::cout << "wlDn1n2Dl_(" << k << "):\t" << wlDn1n2Dl_(k) << std::endl;
    }

    // ====== 临时诊断：强制整条水线导数项为 0 ======
    //wlDn1n2Dl_.setZero();
    //std::cout << "[DEBUG] force all wlDn1n2Dl_ = 0\n";
}

RAO4DTable LinearCumminsTDGF::initRAO4D()
{
    RAO4DTable tab;

    tab.meta.L = ShipCfg.Geometry.Length;
    tab.meta.dirUnit = RAO4DMeta::DirUnit::Rad;
    tab.meta.freqUnit = RAO4DMeta::FreqUnit::RadPerSec;
    tab.meta.freqType = RAO4DMeta::FreqType::IncidentOmega;
    tab.meta.normType = RAO4DMeta::NormType::TransA_RotkA;

    tab.modeIds = SeakeepingCfg.modes;
    tab.Fns = SeakeepingCfg.Fn;

    std::vector<double> betas, oms;
    for (const auto& w : SeakeepingCfg.waves)
    {
        if (auto rw = std::dynamic_pointer_cast<RegularWave>(w))
        {
            betas.push_back(rw->direction());
            oms.push_back(rw->getFreq());
        }
        else if (auto ir = std::dynamic_pointer_cast<IrregularWave>(w))
        {
            // Irregular: representative (direction, spectral-peak ω) so the
            // RAO axes / findExactIndex lookups in buildCaseContext resolve.
            betas.push_back(ir->direction());
            oms.push_back(ir->getFreq());
        }
    }

    std::sort(betas.begin(), betas.end());
    betas.erase(std::unique(betas.begin(), betas.end()), betas.end());

    std::sort(oms.begin(), oms.end());
    oms.erase(std::unique(oms.begin(), oms.end()), oms.end());

    betaAxisDeg = betas;
    omegaAxisIncident = oms;

    tab.dir = betas;
    tab.w = oms;
    if (raoEnabled_)
        tab.finalizeAndAllocate(true);
    return tab;
}

void LinearCumminsTDGF::initGreenTable()
{
    GreenTable::Params p = GreenTable::Params::AdaptiveDefault(300.0);
    p.bmin = 1e-6;
    p.bmax = 300.0;
    p.mmin = 0.0;
    p.mmax = 1.0;
    p.mSplit1 = 0.005;
    p.mSplit2 = 0.01;
    p.mSplit3 = 0.1;
    p.throwOnOutOfRange = true;
    p.autoTune();

    gGreenTable = GreenTable(p);
    gGreenTable.init(filePath + "green_table.bin", G, true);
}

void LinearCumminsTDGF::initTDGFProvider()
{
    const auto& cfg = SeakeepingCfg.GreenFunction;

    auto resolvePath = [&](const std::string& p) -> std::string
        {
            if (p.empty())
                return filePath + "green_table.bin";

            std::filesystem::path path(p);

            if (path.is_absolute())
                return path.string();

            return (std::filesystem::path(filePath) / path).string();
        };

    if (cfg.Method == "Exact")
    {
        tdgfProvider_.initExact();

        std::cout << "[TDGF] Using Exact Greenf::GreenFunctionCal\n";
    }
    else if (cfg.Method == "Interp")
    {
        GreenTable::Params p = GreenTable::Params::AdaptiveDefault(300.0);
        p.bmin = 1.0e-6;
        p.bmax = 300.0;
        p.mmin = 0.0;
        p.mmax = 1.0;
        p.mSplit1 = 0.005;
        p.mSplit2 = 0.01;
        p.mSplit3 = 0.1;
        p.throwOnOutOfRange = true;
        p.autoTune();

        gGreenTable = GreenTable(p);

        const std::string tablePath = resolvePath(
            cfg.TablePath.empty() ? "green_table.bin" : cfg.TablePath);

        gGreenTable.init(tablePath, G, true);
        tdgfProvider_.initInterp(&gGreenTable);

        std::cout << "[TDGF] Using Interp table: " << tablePath << "\n";
    }
    else if (cfg.Method == "ODETable")
    {
        const std::string tablePath = resolvePath(
            cfg.TablePath.empty() ? "tdgf_ode_table.bin" : cfg.TablePath);

        tdgfProvider_.initODETable(tablePath, cfg.RK4Step);

        std::cout << "[TDGF] Using ODETable: " << tablePath
            << ", RK4Step=" << cfg.RK4Step << "\n";
    }
    else
    {
        throw std::runtime_error(
            "LinearCumminsTDGF::initTDGFProvider: unknown GreenFunction.Method: "
            + cfg.Method);
    }
}

void LinearCumminsTDGF::buildEncounterBandForFn(double Fn)
{
    const double Ulocal = Fn * std::sqrt(G * ShipCfg.Geometry.Length);

    std::vector<double> wes;
    wes.reserve(SeakeepingCfg.waves.size());

    auto pushWe = [&](double wInc, double dir)
        {
            try
            {
                const auto enc = calcEncounterInfo(wInc, Ulocal, dir);
                if (enc.we > 1e-8)
                    wes.push_back(enc.we);
            }
            catch (...) { /* following-seas degeneracy: skip this sample */ }
        };

    // Recurse so a CrossWave contributes both its sub-waves' encounter bands.
    std::function<void(const std::shared_ptr<WaveBase>&)> addWave =
        [&](const std::shared_ptr<WaveBase>& w)
        {
            if (auto rw = std::dynamic_pointer_cast<RegularWave>(w))
            {
                pushWe(rw->getFreq(), rw->direction());
            }
            else if (auto ir = std::dynamic_pointer_cast<IrregularWave>(w))
            {
                // Bracket the encounter band over the spectrum's incident range
                // so the online Green grid is fine enough for every component.
                const auto& sp = ir->spectrum();
                const double dir = ir->direction();
                const int ns = 9;
                for (int i = 0; i <= ns; ++i)
                {
                    const double f = sp.omegaLo
                        + (sp.omegaHi - sp.omegaLo) * (double)i / ns;
                    pushWe(f, dir);
                }
            }
            else if (auto cw = std::dynamic_pointer_cast<CrossWave>(w))
            {
                addWave(cw->wave1());
                addWave(cw->wave2());
            }
        };

    for (const auto& w : SeakeepingCfg.waves)
        addWave(w);

    if (wes.empty())
        throw std::runtime_error("LinearCumminsTDGF: no valid encounter frequencies (Regular/Irregular/Cross).");

    for (int i = 0; i < wes.size(); ++i) {
        std::cout << "we[" << i << "]:\t" << wes[i] << std::endl;
    }

    std::sort(wes.begin(), wes.end());

    omegaMin = wes.front();
    omegaMax = wes.back();
}

void LinearCumminsTDGF::setupOnlineGreenGrid(double dt, int TG)
{
    if (dt <= 0.0 || TG <= 0)
        throw std::runtime_error("LinearCumminsTDGF: invalid runtime grid.");

    dtShared_ = dt;
    TGShared_ = TG;
    TSShared_ = TG + 1;
    bufCols = TG + 1;

    SeakeepingCfg.Time.dt = dt;
    SeakeepingCfg.Time.GreenStep = TG;
}

void LinearCumminsTDGF::setupOfflineKernelGridForFn(double Fn)
{
    // ChiTimeGrid 是权威源：启用时，dtShared_ 与 TGShared_ 由它派生。
    // 未启用时退回硬编码默认（保持老行为）。
    const auto& cg = SeakeepingCfg.ChiTimeGrid;
    const bool cgValid = (cg.enabled && cg.dt_fine > 0.0 && !cg.blocks.empty());

    double dtKernel;
    int    TGKernel;

    if (cgValid)
    {
        // dt_fine 直接作为 Green 表的等距步长；TG 取 chi 最远节点的 fine-grid 索引。
        // Green 表需要覆盖到 lag = chiIdxOnFineMax - 1，等价于 TG >= chiIdxOnFineMax。
        int chiIdxOnFineMax = 0;
        for (const auto& b : cg.blocks)
            chiIdxOnFineMax += b.stride * b.count;

        dtKernel = cg.dt_fine;
        TGKernel = chiIdxOnFineMax;     // Green 表覆盖 [0, TG·dt_fine]

        std::cout << "[LinearCumminsTDGF] offline kernel grid driven by ChiTimeGrid:"
                  << " dt_fine=" << dtKernel
                  << ", TG=" << TGKernel
                  << " (tMax=" << TGKernel * dtKernel << "s)\n";
    }
    else
    {
        // 硬编码默认（与原行为完全一致）
        dtKernel = 0.015 / std::sqrt(G / ShipCfg.Geometry.Length);
        TGKernel = static_cast<int>(12 / 0.015);   // = 800

        std::cout << "[LinearCumminsTDGF] offline kernel grid: using built-in defaults"
                  << " (ChiTimeGrid disabled / invalid).\n";
    }

    setupOnlineGreenGrid(dtKernel, TGKernel);

    std::cout << "\n[LinearCumminsTDGF] offline kernel grid:"
        << " dt=" << dtShared_ << "\n"
        << " TG=" << TGShared_ << "\n"
        << " tMem=" << (dtShared_ * TGShared_)
        << "\n\n";

    // 解析非均匀 chi 时间网格（若未启用则 chiIdxOnFine_ 为空，下游走等距路径）。
    // dt_fine == dtShared_ 已自动满足（我们刚把 dtShared_ 设为 dt_fine）。
    buildChiTimeGridFromConfig();
}

// 把 SeakeepingCfg.ChiTimeGrid 展开成内部数组。失败/未启用 -> 数组留空。
void LinearCumminsTDGF::buildChiTimeGridFromConfig()
{
    chiIdxOnFine_.clear();
    chiTimes_.clear();
    chiDtFine_ = 0.0;
    chiMemoryCutoffLag_ = 0;

    const auto& cfg = SeakeepingCfg.ChiTimeGrid;
    if (!cfg.enabled || cfg.dt_fine <= 0.0 || cfg.blocks.empty())
        return;

    // 自洽性：dt_fine 必须与离线 K build 的 dtShared_ 一致（否则 Green 表索引和 chi 节点对不上）
    // 容许小的浮点误差；偏差较大就降级到等距路径并 warn。
    if (std::abs(cfg.dt_fine - dtShared_) > 1e-9 * std::max(1.0, dtShared_))
    {
        std::cerr << "[LinearCumminsTDGF] ChiTimeGrid.dt_fine(" << cfg.dt_fine
                  << ") != dtShared_(" << dtShared_ << "); "
                     "non-uniform path disabled, fallback to uniform.\n";
        return;
    }

    // 展开 blocks
    chiIdxOnFine_.reserve(1 + [&](){
        int n = 0;
        for (auto& b : cfg.blocks) n += b.count;
        return n;
    }());
    chiIdxOnFine_.push_back(0);
    int cur = 0;
    for (const auto& b : cfg.blocks)
    {
        for (int k = 0; k < b.count; ++k)
        {
            cur += b.stride;
            chiIdxOnFine_.push_back(cur);
        }
    }

    chiTimes_.resize(chiIdxOnFine_.size());
    for (std::size_t i = 0; i < chiIdxOnFine_.size(); ++i)
        chiTimes_[i] = static_cast<double>(chiIdxOnFine_[i]) * cfg.dt_fine;

    chiDtFine_          = cfg.dt_fine;
    chiMemoryCutoffLag_ = cfg.memory_cutoff_lag;

    // 必要性：chi 最大节点的 fine 索引必须 ≤ TGShared_，否则 Green 表查表越界
    if (chiIdxOnFine_.back() > TGShared_)
    {
        std::cerr << "[LinearCumminsTDGF] ChiTimeGrid total span ("
                  << chiIdxOnFine_.back() << " · dt_fine = "
                  << chiTimes_.back() << "s) exceeds Green table range ("
                  << TGShared_ << " · dt = "
                  << TGShared_ * dtShared_ << "s); "
                     "non-uniform path disabled, fallback to uniform.\n";
        chiIdxOnFine_.clear();
        chiTimes_.clear();
        chiDtFine_ = 0.0;
        chiMemoryCutoffLag_ = 0;
        return;
    }

    std::cout << "[LinearCumminsTDGF] ChiTimeGrid expanded: "
              << chiIdxOnFine_.size() << " nodes, "
              << "tMax=" << chiTimes_.back() << "s, "
              << "cutoff_lag=" << chiMemoryCutoffLag_
              << " (=" << chiMemoryCutoffLag_ * chiDtFine_ << "s)\n";
}

void LinearCumminsTDGF::allocSharedGreenBuffers()
{
    const int NE = SeakeepingCfg.Panel.NE;
    const int nWL = element->n_WL;

    Green = std::make_shared<std::vector<GMat>>(TGShared_, GMat::Zero(NE, NE));
    Green_dnP = std::make_shared<std::vector<GMat>>(TGShared_, GMat::Zero(NE, NE));
    Green_dnS = std::make_shared<std::vector<GMat>>(TGShared_, GMat::Zero(NE, NE));

    Gw = std::make_shared<std::vector<GMat>>(TGShared_, GMat::Zero(NE, nWL));
    Gw_dnP = std::make_shared<std::vector<GMat>>(TGShared_, GMat::Zero(NE, nWL));
    Gw_dnS = std::make_shared<std::vector<GMat>>(TGShared_, GMat::Zero(NE, nWL));
    Gw_dx = std::make_shared<std::vector<GMat>>(TGShared_, GMat::Zero(NE, nWL));
    Gw_dl = std::make_shared<std::vector<GMat>>(TGShared_, GMat::Zero(NE, nWL));
    Gw_dt = std::make_shared<std::vector<GMat>>(TGShared_, GMat::Zero(NE, nWL));
}

void LinearCumminsTDGF::writeGreenTable(int i_element)
{
    std::string outFile = filePath + "green_table.csv";
    std::ofstream ofs(outFile);
    if (!ofs)
    {
        std::cerr << "Failed to open file for writing: " << outFile << std::endl;
        return;
    }

    ofs << "TimeStep,Time,";
    for (int j = 0; j < SeakeepingCfg.Panel.NE; ++j)
    {
        ofs << "G(" << i_element << ")(" << j << ")";
        if (j < SeakeepingCfg.Panel.NE - 1)
            ofs << ",";
    }

    for (int i = 0; i < TGShared_; ++i)
    {
        ofs << i << "," << dtShared_ * i;
        for (int j = 0; j < SeakeepingCfg.Panel.NE; ++j)
        {
            ofs << "," << (*Green)[static_cast<std::size_t>(i)](i_element, j);
        }
        ofs << "\n";
    }
    std::cout << "Green table written to: " << outFile << std::endl;
}


CaseContext LinearCumminsTDGF::buildCaseContext(int i_case, double Fn, int iFn) const
{
    const double L = ShipCfg.Geometry.Length;
    const double displacement = ShipCfg.Geometry.Displacement;

    CaseContext ctx;
    ctx.i_case = i_case;
    ctx.iFn = iFn;
    ctx.DOF = SeakeepingCfg.DOF;

    ctx.wave = SeakeepingCfg.waves.at(i_case);
    ctx.outputWaveHistory = ctx.wave->outputHistory;

    // Build the incident-elevation WaveField for a single (regular/irregular)
    // wave. Cross sub-waves are decomposed into these.
    auto fieldFor = [](const std::shared_ptr<WaveBase>& w) -> WaveField
        {
            if (auto ir = std::dynamic_pointer_cast<IrregularWave>(w))
                return WaveField::makeIrregular(ir->spectrum(), ir->direction());
            if (auto rw = std::dynamic_pointer_cast<RegularWave>(w))
                return WaveField::makeRegular(rw->getAmp(), rw->direction(),
                    rw->getFreq(), rw->initialPhase());
            throw std::runtime_error(
                "LinearCumminsTDGF: cross sub-wave must be Regular or Irregular.");
        };

    if (auto cw = std::dynamic_pointer_cast<CrossWave>(ctx.wave))
    {
        // Step 4: two crossing systems, each its own direction/kernel/field.
        // Grid/region sizing follows sub-wave 1 (ctx.wave := wave1).
        ctx.subWaves = { cw->wave1(), cw->wave2() };
        ctx.wave = cw->wave1();
        ctx.crossStartX = cw->startX();
        ctx.crossStartY = cw->startY();
    }
    else
    {
        ctx.subWaves = { ctx.wave };
    }

    ctx.reg = std::dynamic_pointer_cast<RegularWave>(ctx.wave);
    auto irr = std::dynamic_pointer_cast<IrregularWave>(ctx.wave);
    if (!ctx.reg && !irr)
        throw std::runtime_error(
            "LinearCumminsTDGF: unsupported wave type (expected Regular/Irregular/Cross).");

    ctx.subFields.clear();
    for (const auto& sw : ctx.subWaves)
        ctx.subFields.push_back(fieldFor(sw));

    if (i_case != 0)
    {
        const double direction = ctx.wave->direction();
        const double old_direction = SeakeepingCfg.waves.at(i_case - 1)->direction();
        ctx.new_direction = std::abs(direction - old_direction) >= 0.1;
    }

    ctx.Amp = ctx.wave->getAmp();   // regular: H; irregular: 2√m0 (display scale)
    ctx.W = ctx.wave->getFreq();  // regular: ω; irregular: spectral-peak ω
    ctx.Fn = Fn;
    ctx.U = Fn * std::sqrt(G * ShipCfg.Geometry.Length);
    ctx.UsquareG = ctx.U * ctx.U / G;

    ctx.dirRad = ctx.wave->direction();
    ctx.we = calcEncounterInfo(ctx.W, ctx.U, ctx.dirRad).we;

    // ctx.waveField mirrors sub-wave 0 (grid sizing + regular-only sentinel).
    // The actual per-sub convolution uses ctx.subFields[ctx.activeSub].
    ctx.waveField = ctx.subFields.front();

    if (raoEnabled_)
    {
        ctx.iDir = findExactIndex(betaAxisDeg, ctx.dirRad, 1e-6);
        ctx.iw = findExactIndex(omegaAxisIncident, ctx.W, 1e-9);
    }
    else
    {
        ctx.iDir = -1; // RAO disabled (irregular/cross): no ω/β axis index
        ctx.iw = -1;
    }

    if (ctx.we <= 1e-10)
        throw std::runtime_error("LinearCumminsTDGF: encounter frequency <= 0.");

    //const double T = 2.0 * PI / ctx.we;
    //const double T = 2.0 * PI / omegaMax;

    //ctx.dt = T / 60.0;

    //ctx.TG = std::max(8, static_cast<int>(std::ceil(
    //    SeakeepingCfg.Time.GreenCircle * T / ctx.dt)));
    //ctx.tMot = std::max(1, static_cast<int>(std::ceil(
    //    SeakeepingCfg.Time.PreCircle * T / ctx.dt)));

    //const int runSteps = std::max(8, static_cast<int>(std::ceil(
    //    SeakeepingCfg.Time.TimeCircle * T / ctx.dt)));
    //ctx.TS = ctx.tMot + runSteps + 2;

    //ctx.dt_const = ctx.dt;
    //ctx.TG_const = ctx.TG * ctx.dt;
    //ctx.tMot_const = ctx.tMot * ctx.dt;

    //const double T = 2.0 * PI / omegaMax;
    //ctx.dt_const = T / 100.0;

    // Physics-driven (dt, tMot) per case — see FKImpulseKernelConfig and deriveKernelGrid().
    // Region-adaptive defaults: finer dt + shorter τ window in head seas; longer cap in
    // following (especially F2) so the impulse tail is not truncated while |ω_e| is small.
    const KernelGrid grid =
        deriveKernelGrid(ctx, L, SeakeepingCfg.FKImpulseKernel);
    ctx.dt_const = grid.dt;
    ctx.tMot = grid.tMot;
    ctx.tMot_const = grid.memory;
    ctx.TG_const = 2 * ctx.tMot_const;
    ctx.TG = 2 * ctx.tMot;

    const double T = 2.0 * PI / ctx.we;
    ctx.dt = std::min(T / 100.0, ctx.dt_const);
    //ctx.dt = ctx.dt_const;
    ctx.TS = std::max(
        static_cast<int>(std::ceil(10.0 * ctx.tMot_const / ctx.dt)),
        static_cast<int>(std::ceil(50.0 * T / ctx.dt))
    );

    const double a = ctx.Amp;

    ctx.NE = SeakeepingCfg.Panel.NE;
    ctx.n_WL = element->n_WL;

    ctx.a_scale[0] = 1.0;
    ctx.a_scale[1] = 1.0;
    ctx.a_scale[2] = rho * G * displacement * a / L;
    ctx.a_scale[3] = 1.0;
    ctx.a_scale[4] = rho * G * displacement * a;
    ctx.a_scale[5] = 1.0;

    // Tag = "_Fn{Fn}_" + condition descriptor + "_LinearCumminsKernel".
    // Wave::conditionTag names both systems of a crossing sea ("X_{s1}__{s2}")
    // and the spectrum of an irregular sea, so the condition is identifiable
    // from the file name alone. A single regular wave also carries its
    // encounter frequency (_We), which is unambiguous there.
    std::ostringstream ss;
    ss << "_Fn" << shortNum(Fn)
       << "_" << Wave::conditionTag(SeakeepingCfg.waves.at(i_case));

    if (ctx.reg && ctx.subWaves.size() == 1)
        ss << "_We" << shortNum(ctx.we);

    ss << "_LinearCumminsKernel";
    ctx.tag = ss.str();

    //ctx.ExcitingForce.resize(ctx.TS - ctx.tMot, ctx.DOF);
    ctx.ExcitingForce.resize(ctx.TS, ctx.DOF);
    ctx.FKForceHist.resize(ctx.TS, ctx.DOF);
    ctx.DiffForceHist.resize(ctx.TS, ctx.DOF);
    ctx.ExcitingForce.setZero();
    ctx.FKForceHist.setZero();
    ctx.DiffForceHist.setZero();

    std::cout << "T:\t" << 10.0 * T << "\n"
        << "TG:\t" << ctx.TG << "\n"
        << "TS:\t" << ctx.TS << "\n"
        << "dt_const:\t" << ctx.dt_const << "\n"
        << "dt:\t" << ctx.dt << "\n";

    //ctx.F3_amplitude = F3_amp[i_case];
    //ctx.F5_amplitude = F5_amp[i_case];

    return ctx;
}

std::string LinearCumminsTDGF::keyDouble(double x)
{
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(4) << x;
    std::string s = ss.str();
    for (char& c : s)
    {
        if (c == '.') c = 'p';
        if (c == '-') c = 'm';
    }
    return s;
}

std::string LinearCumminsTDGF::makeKernelCachePath(double Fn) const
{
    // 模态码连成一串数字（modes 通常是单位数，如 {2,4} -> "24"）。
    const std::string modesTag = [&]()
        {
            std::ostringstream ss;
            for (std::size_t i = 0; i < SeakeepingCfg.modes.size(); ++i)
                ss << SeakeepingCfg.modes[i];
            return ss.str();
        }();

    std::filesystem::path dir = std::filesystem::path(filePath) / "kernel_cache";

    // <ShipName>_Fn<Fn:.3f>_L<Lpp:.4f>_dof<modes>.bin
    // 例：S175_Fn0.275_L3.0337_dof24.bin
    // 不再附 dt / TG / NE：这些都做成"文件内容自描述"，下游 resample 自适应。
    std::ostringstream name;
    name << ShipCfg.Name
        << "_Fn" << std::fixed << std::setprecision(3) << Fn
        << "_L"  << std::fixed << std::setprecision(4) << ShipCfg.Geometry.Length
        << "_dof" << modesTag
        << ".bin";

    return (dir / name.str()).string();
}

bool LinearCumminsTDGF::tryLoadKernelCache(double Fn, RadiationKernelData& kernel) const
{
    const std::string file = makeKernelCachePath(Fn);
    if (!RadiationKernelCache::load(file, kernel))
        return false;

    // ---- 关键身份校验（必须严格匹配，不一致直接拒绝命中、强制重建）----
    if (std::abs(kernel.Fn - Fn) > 1e-9) return false;
    if (kernel.DOF != SeakeepingCfg.DOF) return false;
    if (kernel.modes != SeakeepingCfg.modes) return false;

    // dt / TG / NE 仍然要求一致：避免陈旧 cache 被误命中导致 K 截断不足
    // （文件名虽然不再含 dt/TG，但内容里有，所以重建会自动覆盖旧文件）
    if (std::abs(kernel.dt - dtShared_) > 1e-12) return false;
    if (kernel.TG != TGShared_) return false;
    if (static_cast<int>(kernel.Klag.size()) != TGShared_ + 1) return false;

    // 矩阵尺寸校验
    if (kernel.A_inf.rows() != SeakeepingCfg.DOF || kernel.A_inf.cols() != SeakeepingCfg.DOF) return false;
    if (kernel.B.rows() != SeakeepingCfg.DOF || kernel.B.cols() != SeakeepingCfg.DOF) return false;
    if (kernel.C_prime.rows() != SeakeepingCfg.DOF || kernel.C_prime.cols() != SeakeepingCfg.DOF) return false;

    // Klag_times 长度必须与 Klag 一致（v7 老格式由 cache 内部自动补齐 i*dt）
    if (kernel.Klag_times.size() != kernel.Klag.size())
    {
        std::cerr << "[LinearCumminsTDGF] kernel cache rejected: "
                     "Klag_times/Klag size mismatch in " << file << "\n";
        return false;
    }

    // 文件指纹与当前案例一致性（仅 warn，不拒绝——方便跨案例搬运同一份 K）
    if (!kernel.ShipName.empty() && kernel.ShipName != ShipCfg.Name)
    {
        std::cerr << "[LinearCumminsTDGF] kernel cache WARNING: ShipName mismatch "
                  << "(file=" << kernel.ShipName << " cfg=" << ShipCfg.Name << ")\n";
    }
    if (kernel.Lpp > 0.0 &&
        std::abs(kernel.Lpp - ShipCfg.Geometry.Length) > 1e-6)
    {
        std::cerr << "[LinearCumminsTDGF] kernel cache WARNING: Lpp mismatch "
                  << "(file=" << kernel.Lpp
                  << " cfg=" << ShipCfg.Geometry.Length << ")\n";
    }

    std::cout << "[LinearCumminsTDGF] kernel cache loaded: " << file
              << " (N=" << kernel.Klag.size()
              << ", tMax=" << kernel.Klag_times.back() << "s)\n";
    return true;
}

void LinearCumminsTDGF::buildKernelCacheForFn(double Fn, RadiationKernelData& kernel)
{
    buildKernelDirectTimeDomain(Fn, kernel);

    kernel.Fn = Fn;
    kernel.U = Fn * std::sqrt(G * ShipCfg.Geometry.Length);
    kernel.dt = dtShared_;
    kernel.TG = static_cast<int>(kernel.Klag.size()) - 1;
    kernel.NE = SeakeepingCfg.Panel.NE;
    kernel.DOF = SeakeepingCfg.DOF;
    kernel.modes = SeakeepingCfg.modes;
    kernel.ShipName = ShipCfg.Name;
    kernel.Lpp = ShipCfg.Geometry.Length;

    // 等距路径下，确保 Klag_times = i*dt 与 Klag 同步存出。
    // buildKernelDirectTimeDomain 退出时应该已经填好；这里做兜底。
    if (kernel.Klag_times.size() != kernel.Klag.size())
    {
        kernel.Klag_times.resize(kernel.Klag.size());
        for (std::size_t i = 0; i < kernel.Klag.size(); ++i)
            kernel.Klag_times[i] = static_cast<double>(i) * dtShared_;
    }

    const std::string file = makeKernelCachePath(Fn);
    if (!RadiationKernelCache::save(file, kernel))
        throw std::runtime_error("LinearCumminsTDGF: failed to save radiation kernel cache.");

    std::cout << "[LinearCumminsTDGF] kernel cache saved: " << file << "\n";
    writeKernelDebugCsv(Fn, kernel);
}

void LinearCumminsTDGF::writeKernelDebugCsv(double Fn, const RadiationKernelData& kernel) const
{
    std::filesystem::create_directories(std::filesystem::path(filePath) / "kernel_cache");

    std::ostringstream ss;
    ss << filePath << "kernel_cache/tdkernel_ref34_debug_Fn" << keyDouble(Fn) << ".csv";

    std::ofstream out(ss.str());
    if (!out.is_open()) return;

    const int D = SeakeepingCfg.DOF;

    const double L = ShipCfg.Geometry.Length;
    const double V = ShipCfg.Geometry.Displacement;   // 这里默认 Displacement 存的是排水体积
    const double forceScale = rho * G * V;
    const double momentScale = rho * G * V * L;

    out << "# A_inf (dimensional)\n";
    for (int i = 0; i < D; ++i)
    {
        for (int j = 0; j < D; ++j)
        {
            if (j) out << ",";
            out << kernel.A_inf(i, j);
        }
        out << "\n";
    }

    out << "# B (dimensional)\n";
    for (int i = 0; i < D; ++i)
    {
        for (int j = 0; j < D; ++j)
        {
            if (j) out << ",";
            out << kernel.B(i, j);
        }
        out << "\n";
    }

    out << "# C_prime (dimensional)\n";
    for (int i = 0; i < D; ++i)
    {
        for (int j = 0; j < D; ++j)
        {
            if (j) out << ",";
            out << kernel.C_prime(i, j);
        }
        out << "\n";
    }

    out << "t_nd";
    for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j)
            out << ",Knd_" << SeakeepingCfg.modes[i] << SeakeepingCfg.modes[j];
    out << "\n";

    // 非均匀路径下 Klag_times[m] 是真实节点时间；等距路径下兜底 = m*kernel.dt。
    const double tScale = std::sqrt(G / L);
    const double non_dt = kernel.dt * tScale;
    const bool haveTimes =
        kernel.Klag_times.size() == kernel.Klag.size() && !kernel.Klag_times.empty();
    const int Nrows = static_cast<int>(kernel.Klag.size());

    for (int m = 0; m < Nrows; ++m)
    {
        const double t_nd = haveTimes
            ? (kernel.Klag_times[static_cast<std::size_t>(m)] * tScale)
            : (m * non_dt);
        out << t_nd;

        for (int i = 0; i < D; ++i)
        {
            const int rowMode = SeakeepingCfg.modes[i];
            const double rowScale = (rowMode >= 0 && rowMode <= 2) ? forceScale : momentScale;

            for (int j = 0; j < D; ++j)
            {
                const int colMode = SeakeepingCfg.modes[j];
                const double colScale = (colMode >= 0 && colMode <= 2) ? L : 1.0;

                const double Knd =
                    kernel.Klag[static_cast<std::size_t>(m)](i, j) * colScale / rowScale;

                out << "," << Knd;
            }
        }
        out << "\n";
    }
}

void LinearCumminsTDGF::writeVectorHistoryDebugCsv(
    const std::string& quantityName,
    double Fn,
    int srcMode,
    const std::vector<Eigen::VectorXd>& history) const
{
    std::filesystem::create_directories(std::filesystem::path(filePath) / "kernel_cache");

    std::ostringstream ss;
    ss << filePath
        << "kernel_cache/" << quantityName << "_debug"
        << "_Fn" << keyDouble(Fn)
        << "_mode" << srcMode
        << "_dt" << keyDouble(dtShared_)
        << "_TG" << TGShared_
        << ".csv";

    std::ofstream out(ss.str());
    if (!out.is_open())
    {
        std::cerr << "[LinearCumminsTDGF] cannot open "
            << quantityName << " debug csv: "
            << ss.str() << "\n";
        return;
    }

    const double L = ShipCfg.Geometry.Length;
    const double tScale = std::sqrt(G / L);   // t_nd = t * sqrt(g/L)
    const int NE = SeakeepingCfg.Panel.NE;

    out << std::setprecision(17);

    out << "step,t,t_nd," << quantityName << "_norm,"
        << quantityName << "_jump_norm";
    for (int p = 0; p < NE; ++p)
        out << ",element_" << p;
    out << "\n";

    // 非均匀路径下用 chi 节点真实时间；等距路径下回退到 n*dtShared_。
    const bool useNUTimes =
        !chiTimes_.empty() && chiTimes_.size() == history.size();

    for (int n = 0; n < static_cast<int>(history.size()); ++n)
    {
        const double t = useNUTimes
            ? chiTimes_[static_cast<std::size_t>(n)]
            : (n * dtShared_);
        const double t_nd = t * tScale;
        const Eigen::VectorXd& v = history[static_cast<std::size_t>(n)];

        double jump_norm = 0.0;
        if (n >= 1)
            jump_norm = (history[static_cast<std::size_t>(n)]
                - history[static_cast<std::size_t>(n - 1)]).norm();

        out << n << "," << t << "," << t_nd << ","
            << v.norm() << "," << jump_norm;

        for (int p = 0; p < v.size(); ++p)
            out << "," << v(p);

        out << "\n";
    }

    std::cout << "[LinearCumminsTDGF] " << quantityName
        << " debug csv saved: " << ss.str() << "\n";
}

void LinearCumminsTDGF::writeChiDebugCsv(
    double Fn,
    int srcMode,
    const std::vector<Eigen::VectorXd>& chi_history) const
{
    writeVectorHistoryDebugCsv("chi", Fn, srcMode, chi_history);
}

void LinearCumminsTDGF::writeDchiDebugCsv(
    double Fn,
    int srcMode,
    const std::vector<Eigen::VectorXd>& dchi_history) const
{
    writeVectorHistoryDebugCsv("dchi", Fn, srcMode, dchi_history);
}

std::vector<Eigen::VectorXd> LinearCumminsTDGF::smoothStartupHistory(
    const std::vector<Eigen::VectorXd>& raw_history,
    int fitSteps,
    int blendSteps,
    int polyDeg,
    double jaggedRatioThreshold) const
{
    // ---------------------------------------------------------------------
    // Conservative adaptive smoothing for chi history.
    //
    // Purpose:
    //   1. Do not assume the oscillation only occurs at startup.
    //   2. Detect local high-frequency abnormal oscillations automatically.
    //   3. Smooth chi mildly before differentiating it.
    //   4. Avoid aggressive global smoothing that may destroy physical memory effects.
    //
    // Method:
    //   - For each panel p, analyse chi_p(t).
    //   - Use robust statistics of first and second differences.
    //   - Only if local curvature/jump is abnormal, apply a limited 5-point
    //     quadratic Savitzky-Golay correction.
    //
    // Notes:
    //   - The old parameters are kept for compatibility with the existing header.
    //   - fitSteps / blendSteps / polyDeg are no longer used as fixed fitting controls.
    //   - jaggedRatioThreshold is used only as a sensitivity reference.
    // ---------------------------------------------------------------------

    (void)fitSteps;
    (void)blendSteps;
    (void)polyDeg;

    if (raw_history.empty())
        return raw_history;

    const int N = static_cast<int>(raw_history.size());
    const int NE = static_cast<int>(raw_history.front().size());

    std::vector<Eigen::VectorXd> smooth_history = raw_history;

    if (N < 7 || NE <= 0)
        return smooth_history;

    const double eps = 1.0e-14;

    auto medianValue = [](std::vector<double> v) -> double
        {
            v.erase(
                std::remove_if(v.begin(), v.end(),
                    [](double x)
                    {
                        return !std::isfinite(x);
                    }),
                v.end());

            if (v.empty())
                return 0.0;

            std::sort(v.begin(), v.end());

            const std::size_t n = v.size();
            if (n % 2 == 1)
                return v[n / 2];

            return 0.5 * (v[n / 2 - 1] + v[n / 2]);
        };

    auto robustSigma = [&](const std::vector<double>& v, double med) -> double
        {
            std::vector<double> dev;
            dev.reserve(v.size());

            for (double x : v)
            {
                if (std::isfinite(x))
                    dev.push_back(std::abs(x - med));
            }

            const double mad = medianValue(dev);

            // 1.4826 * MAD is a robust estimate of standard deviation.
            return 1.4826 * mad;
        };

    auto clamp01 = [](double x) -> double
        {
            return std::max(0.0, std::min(1.0, x));
        };

    // Conservative settings.
    // Larger kDetect -> fewer points are smoothed.
    // Larger maxAlpha -> stronger smoothing.
    const double sensitivity = std::max(1.0, jaggedRatioThreshold);

    const double kDetect =
        std::max(4.0, 3.0 * sensitivity);      // normally about 4.05 when threshold=1.35

    const double maxAlpha = 0.30;              // never replace more than 30%
    const double minAbnormalRatio = 1.20;      // avoid reacting to tiny numerical noise

    int changedCount = 0;
    double alphaSum = 0.0;

    // Temporary scalar arrays for one panel.
    std::vector<double> y(static_cast<std::size_t>(N), 0.0);
    std::vector<double> dy(static_cast<std::size_t>(N - 1), 0.0);
    std::vector<double> ddy(static_cast<std::size_t>(N - 2), 0.0);

    for (int p = 0; p < NE; ++p)
    {
        bool hasBadValue = false;

        for (int n = 0; n < N; ++n)
        {
            const double v = raw_history[static_cast<std::size_t>(n)](p);

            if (!std::isfinite(v))
            {
                hasBadValue = true;
                y[static_cast<std::size_t>(n)] = 0.0;
            }
            else
            {
                y[static_cast<std::size_t>(n)] = v;
            }
        }

        // If non-finite values appear, first repair them by neighbour interpolation.
        if (hasBadValue)
        {
            for (int n = 0; n < N; ++n)
            {
                const double v = raw_history[static_cast<std::size_t>(n)](p);

                if (std::isfinite(v))
                    continue;

                int l = n - 1;
                int r = n + 1;

                while (l >= 0 &&
                    !std::isfinite(raw_history[static_cast<std::size_t>(l)](p)))
                    --l;

                while (r < N &&
                    !std::isfinite(raw_history[static_cast<std::size_t>(r)](p)))
                    ++r;

                if (l >= 0 && r < N)
                    y[static_cast<std::size_t>(n)] =
                    0.5 * (raw_history[static_cast<std::size_t>(l)](p) +
                        raw_history[static_cast<std::size_t>(r)](p));
                else if (l >= 0)
                    y[static_cast<std::size_t>(n)] =
                    raw_history[static_cast<std::size_t>(l)](p);
                else if (r < N)
                    y[static_cast<std::size_t>(n)] =
                    raw_history[static_cast<std::size_t>(r)](p);
                else
                    y[static_cast<std::size_t>(n)] = 0.0;
            }
        }

        // First difference and second difference.
        for (int n = 0; n < N - 1; ++n)
        {
            dy[static_cast<std::size_t>(n)] =
                std::abs(y[static_cast<std::size_t>(n + 1)] -
                    y[static_cast<std::size_t>(n)]);
        }

        for (int n = 1; n < N - 1; ++n)
        {
            ddy[static_cast<std::size_t>(n - 1)] =
                std::abs(y[static_cast<std::size_t>(n + 1)]
                    - 2.0 * y[static_cast<std::size_t>(n)]
                    + y[static_cast<std::size_t>(n - 1)]);
        }

        const double medDy = medianValue(dy);
        const double sigDy = robustSigma(dy, medDy);

        const double medDdy = medianValue(ddy);
        const double sigDdy = robustSigma(ddy, medDdy);

        const double scaleY = std::max(
            eps,
            medianValue(std::vector<double>(y.begin(), y.end())));

        const double thrDy =
            medDy + kDetect * std::max(sigDy, 1.0e-10 * std::max(1.0, std::abs(scaleY)));

        const double thrDdy =
            medDdy + kDetect * std::max(sigDdy, 1.0e-10 * std::max(1.0, std::abs(scaleY)));

        // If the whole sequence is very smooth, skip this panel.
        if (thrDy <= eps && thrDdy <= eps)
            continue;

        // Keep both ends unchanged.
        smooth_history[0](p) = y[0];
        smooth_history[static_cast<std::size_t>(N - 1)](p) =
            y[static_cast<std::size_t>(N - 1)];

        for (int n = 1; n < N - 1; ++n)
        {
            const double jumpLeft =
                std::abs(y[static_cast<std::size_t>(n)] -
                    y[static_cast<std::size_t>(n - 1)]);

            const double jumpRight =
                std::abs(y[static_cast<std::size_t>(n + 1)] -
                    y[static_cast<std::size_t>(n)]);

            const double localJump = std::max(jumpLeft, jumpRight);

            const double localCurv =
                std::abs(y[static_cast<std::size_t>(n + 1)]
                    - 2.0 * y[static_cast<std::size_t>(n)]
                    + y[static_cast<std::size_t>(n - 1)]);

            const double ratioJump =
                localJump / std::max(eps, thrDy);

            const double ratioCurv =
                localCurv / std::max(eps, thrDdy);

            const double abnormal =
                std::max(ratioJump, ratioCurv);

            if (abnormal < minAbnormalRatio)
            {
                smooth_history[static_cast<std::size_t>(n)](p) =
                    y[static_cast<std::size_t>(n)];
                continue;
            }

            // Convert abnormal level to a mild smoothing weight.
            // abnormal = 1.2 -> alpha close to 0
            // abnormal >= 3.0 -> alpha close to maxAlpha
            double alpha =
                maxAlpha * clamp01((abnormal - minAbnormalRatio) / (3.0 - minAbnormalRatio));

            // 5-point quadratic Savitzky-Golay smoothing when possible.
            // Coefficients: [-3, 12, 17, 12, -3] / 35.
            double sg = y[static_cast<std::size_t>(n)];

            if (n >= 2 && n <= N - 3)
            {
                sg =
                    (-3.0 * y[static_cast<std::size_t>(n - 2)]
                        + 12.0 * y[static_cast<std::size_t>(n - 1)]
                        + 17.0 * y[static_cast<std::size_t>(n)]
                        + 12.0 * y[static_cast<std::size_t>(n + 1)]
                        - 3.0 * y[static_cast<std::size_t>(n + 2)]) / 35.0;
            }
            else
            {
                // Near boundaries, use a very mild 3-point smoother.
                sg =
                    0.25 * y[static_cast<std::size_t>(n - 1)]
                    + 0.50 * y[static_cast<std::size_t>(n)]
                    + 0.25 * y[static_cast<std::size_t>(n + 1)];

                alpha *= 0.5;
            }

            double correction =
                alpha * (sg - y[static_cast<std::size_t>(n)]);

            // Limit absolute correction to avoid destroying physically meaningful peaks.
            const double localAmplitude =
                std::max({
                    std::abs(y[static_cast<std::size_t>(n - 1)]),
                    std::abs(y[static_cast<std::size_t>(n)]),
                    std::abs(y[static_cast<std::size_t>(n + 1)]),
                    eps
                    });

            const double maxCorrection =
                0.20 * localJump + 0.02 * localAmplitude;

            if (std::abs(correction) > maxCorrection)
            {
                correction =
                    std::copysign(maxCorrection, correction);
            }

            const double newValue =
                y[static_cast<std::size_t>(n)] + correction;

            smooth_history[static_cast<std::size_t>(n)](p) = newValue;

            if (std::abs(correction) > 0.0)
            {
                ++changedCount;
                alphaSum += alpha;
            }
        }
    }

    // Optional diagnostic.
    if (changedCount > 0)
    {
        std::cout << "[smoothStartupHistory] adaptive robust smoothing applied."
            << " changed_values=" << changedCount
            << " mean_alpha=" << alphaSum / static_cast<double>(changedCount)
            << " max_alpha=" << maxAlpha
            << " kDetect=" << kDetect
            << "\n";
    }
    else
    {
        std::cout << "[smoothStartupHistory] no abnormal high-frequency oscillation detected.\n";
    }

    return smooth_history;
}



std::vector<Eigen::VectorXd> LinearCumminsTDGF::differentiateHistory(
    const std::vector<Eigen::VectorXd>& history) const
{
    if (history.empty())
        return history;

    const int N = static_cast<int>(history.size());
    const int NE = static_cast<int>(history.front().size());

    std::vector<Eigen::VectorXd> d(
        static_cast<std::size_t>(N),
        Eigen::VectorXd::Zero(NE));

    const double dt = dtShared_;
    if (N < 2 || dt <= 0.0)
        return d;

    if (N == 2)
    {
        d[0] = (history[1] - history[0]) / dt;
        d[1] = d[0];
        return d;
    }

    // 二阶单边端点
    d[0] =
        (-3.0 * history[0]
            + 4.0 * history[1]
            - 1.0 * history[2]) / (2.0 * dt);

    d[static_cast<std::size_t>(N - 1)] =
        (3.0 * history[static_cast<std::size_t>(N - 1)]
            - 4.0 * history[static_cast<std::size_t>(N - 2)]
            + 1.0 * history[static_cast<std::size_t>(N - 3)]) / (2.0 * dt);

    // 中间点
    for (int n = 1; n <= N - 2; ++n)
    {
        if (n >= 2 && n <= N - 3)
        {
            // 四阶中心差分
            d[static_cast<std::size_t>(n)] =
                (-1.0 * history[static_cast<std::size_t>(n + 2)]
                    + 8.0 * history[static_cast<std::size_t>(n + 1)]
                    - 8.0 * history[static_cast<std::size_t>(n - 1)]
                    + 1.0 * history[static_cast<std::size_t>(n - 2)])
                / (12.0 * dt);
        }
        else
        {
            // 靠近端点用二阶中心差分
            d[static_cast<std::size_t>(n)] =
                (history[static_cast<std::size_t>(n + 1)]
                    - history[static_cast<std::size_t>(n - 1)])
                / (2.0 * dt);
        }
    }

    return d;
}


//std::vector<Eigen::MatrixXd> LinearCumminsTDGF::resampleKernelToDt(
//    const RadiationKernelData& kernel,
//    double dtOnline) const
//{
//    const int D = kernel.DOF;
//    const double tMem = kernel.TG * kernel.dt;
//
//    const int M = std::max(
//        1,
//        static_cast<int>(std::floor(tMem / dtOnline + 1e-12)));
//
//    std::vector<Eigen::MatrixXd> out(
//        static_cast<std::size_t>(M + 1),
//        Eigen::MatrixXd::Zero(D, D));
//
//    auto sampleAt = [&](double tau) -> Eigen::MatrixXd
//    {
//        if (tau <= 0.0)
//            return kernel.Klag[0];
//
//        if (tau > tMem)
//            return Eigen::MatrixXd::Zero(D, D);
//
//        if (std::abs(tau - tMem) <= 1e-12 * std::max(1.0, tMem))
//            return kernel.Klag[static_cast<std::size_t>(kernel.TG)];
//
//        const double x = tau / kernel.dt;
//        int i0 = static_cast<int>(std::floor(x));
//        double a = x - static_cast<double>(i0);
//
//        if (i0 < 0)
//            return kernel.Klag[0];
//
//        if (i0 >= kernel.TG)
//            return kernel.Klag[static_cast<std::size_t>(kernel.TG)];
//
//        const auto& K0m = kernel.Klag[static_cast<std::size_t>(i0)];
//        const auto& K1m = kernel.Klag[static_cast<std::size_t>(i0 + 1)];
//
//        return (1.0 - a) * K0m + a * K1m;
//    };
//
//    for (int m = 0; m <= M; ++m)
//    {
//        out[static_cast<std::size_t>(m)] = sampleAt(m * dtOnline);
//    }
//
//    return out;
//}

std::vector<Eigen::MatrixXd>
LinearCumminsTDGF::resampleKernelToDt(
    const RadiationKernelData& kernel,
    double dtNew) const
{
    if (dtNew <= 0.0)
        throw std::runtime_error("resampleKernelToDt: dtNew <= 0");

    if (kernel.Klag.empty())
        throw std::runtime_error("resampleKernelToDt: empty Klag");

    // Klag_times 是新版的"时间网格"——等距/非均匀都走它。
    // 老数据如果没填，外层（cache 加载、build 收尾、override）已经统一兜底。
    if (kernel.Klag_times.size() != kernel.Klag.size())
        throw std::runtime_error(
            "resampleKernelToDt: Klag_times size must equal Klag size");

    const int D = kernel.DOF;
    const auto& T = kernel.Klag_times;
    const int oldN = static_cast<int>(T.size()) - 1;
    const double tMax = T.back();

    if (oldN < 1)
        return kernel.Klag;

    const int newN = std::max(1, static_cast<int>(std::ceil(tMax / dtNew)));

    std::vector<Eigen::MatrixXd> Knew(
        static_cast<std::size_t>(newN + 1),
        Eigen::MatrixXd::Zero(D, D));

    // 在 T 上做线性插值；T 严格递增；t<=0 保留 K(0)，t>=tMax 保留末点。
    // 等距网格下，这与原来的 floor(t/dt) 实现数值完全等价。
    auto sampleK = [&](double t) -> Eigen::MatrixXd
        {
            if (t <= 0.0)
                return kernel.Klag.front();
            if (t >= tMax)
                return kernel.Klag.back();

            auto it = std::upper_bound(T.begin(), T.end(), t);
            const int j1 = static_cast<int>(it - T.begin());
            const int j0 = j1 - 1;
            const double t0 = T[static_cast<std::size_t>(j0)];
            const double t1 = T[static_cast<std::size_t>(j1)];
            const double w = (t - t0) / std::max(1e-30, (t1 - t0));

            const Eigen::MatrixXd& K0 = kernel.Klag[static_cast<std::size_t>(j0)];
            const Eigen::MatrixXd& K1 = kernel.Klag[static_cast<std::size_t>(j1)];
            return (1.0 - w) * K0 + w * K1;
        };

    // 关键：明确保留 K(0)
    Knew[0] = kernel.Klag.front();

    for (int m = 1; m <= newN; ++m)
    {
        const double t = m * dtNew;
        Knew[static_cast<std::size_t>(m)] = sampleK(t);
    }

    return Knew;
}



void LinearCumminsTDGF::solve()
{
    raoEnabled_ = true;
    for (const auto& w : SeakeepingCfg.waves)
        if (!std::dynamic_pointer_cast<RegularWave>(w)) { raoEnabled_ = false; break; }
    if (!raoEnabled_)
        std::cout << "[RAO4D] irregular/cross sea -> frequency-domain RAO disabled.\n";

    RAO4DTable tab = initRAO4D();

    std::vector<double> non_we;
    double non_omiga = sqrt(G / ShipCfg.Geometry.Length);

    for (int iFn = 0; iFn < static_cast<int>(SeakeepingCfg.Fn.size()); ++iFn)
    {
        bool new_Fn = true;
        fkImpulsePrepared_ = false;
        fkImpulsePrepDtConst_ = -1.0;
        fkImpulsePrepTMot_ = -1;

        const double Fn = SeakeepingCfg.Fn[iFn];
        U = Fn * std::sqrt(G * ShipCfg.Geometry.Length);
        UsquareG = U * U / G;

        setupOfflineKernelGridForFn(Fn);

        RadiationKernelData kernel;

        if (!tryLoadKernelCache(Fn, kernel))
        {
            std::cout << "[LinearCumminsTDGF] load kernel cache failed, start computing ..." << "\n";

            ensureGreenTablesMatchOnlineGrid();

            //writeGreenTable(25);

            buildKernelCacheForFn(Fn, kernel);
            std::cout << "build Kernel cache for Fn=\t" << Fn << "\tdone.\n";
        }

        // ---- 辐射核后处理：两套策略二选一，注释切换 ----
        // 策略 A（原）：tail-split，扣常数并入 C_prime。
        // SHIPSIM_TAILSPLIT=0 可整体关闭（水线项符号扫描用：tail-split 按尾巴
        // 质量自适应触发，各符号组合被修正的程度不同，会污染 A/B(ω) 对比）。
        {
            const char* tsEnv = std::getenv("SHIPSIM_TAILSPLIT");
            if (tsEnv && std::string(tsEnv) == "0")
                std::cout << "[LinearCumminsTDGF] tail-split postprocess DISABLED"
                             " (SHIPSIM_TAILSPLIT=0, raw-kernel scan mode)\n";
            else
                postprocessRadiationKernel(Fn, kernel);
        }
        // 策略 B（新）：不扣常数、不动 A_inf/B/C_prime，时间轴延长一倍并补尾部常数。
        //postprocessRadiationKernel2(Fn, kernel);

        // 强制重写 K debug CSV：cache 命中时 build 函数不会被调用，
        // 但 CSV 仍要按当前的 Klag_times 和后处理结果输出，否则坐标会陈旧。
        writeKernelDebugCsv(Fn, kernel);

        TimeToFrequency(filePath, Fn, PI, kernel);

        //applyRadiationKernelOverride(kernel);

        for (int i_case = 0; i_case < static_cast<int>(SeakeepingCfg.waves.size()); ++i_case)
        {
            CaseContext ctx = buildCaseContext(i_case, Fn, iFn);
            ctx.new_Fn = new_Fn;
            non_we.push_back(ctx.we / non_omiga);

            std::cout << "\n solving wave case: we = " << ctx.we << "\n";

            U = ctx.U;
            UsquareG = ctx.UsquareG;

            //setupOnlineGreenGrid(ctx.dt, ctx.TG);
            setupOnlineGreenGrid(ctx.dt_const, ctx.TG);

            eforce.resize(ctx.NE, 2 * ctx.tMot + 1, ctx.TG);

            computeExciting(ctx);

            //solveKernelCase_RK4(ctx, kernel, tab);

            //ctx.ExcitingForce *= 0.5;

            solveKernelCase(ctx, kernel, tab);

            new_Fn = false;
        }
    }

    if (raoEnabled_)
        tab.writeCSV(filePath + "RAO_linear_cummins_kernel.csv", non_we, 12, true);
}


void LinearCumminsTDGF::computeGreenTables()
{
    std::cout << "\nstart computing GreenFunction ...\n\n";
#pragma omp parallel
    {
        shipsim::EigenSingleThreadGuard _g;
        Gsinteg green(element, U);
#pragma omp for schedule(static)
        for (int tN = 0; tN < TGShared_; ++tN)
        {
            const double tn = (tN + 1) * dtShared_;
            GreenFunction(tN, tn, green);

            std::cout << "GreenFunction step:\t" << tN + 1 << "/" << TGShared_ << "\n";
        }
    }

    //zeroWL(TGShared_);

    std::cout << "[LinearCumminsTDGF] Green tables done. dt=" << dtShared_
        << " TG=" << TGShared_ << "\n";
}

void LinearCumminsTDGF::ensureGreenTablesMatchOnlineGrid()
{
    const bool stale =
        !Green || !Green_dnP || !Gw ||
        static_cast<int>(Green->size()) != TGShared_ ||
        std::abs(greenBuiltForDt_ - dtShared_) > 1.0e-14;
    if (!stale)
        return;

    allocSharedGreenBuffers();
    computeGreenTables();
    greenBuiltForDt_ = dtShared_;
}


//void LinearCumminsTDGF::GreenFunction(int tN, double tn, Gsinteg& green)
//{
//    const int NE = SeakeepingCfg.Panel.NE;
//    const int nWL = element->n_WL;
//
//    GreenData Gdata;
//    double nV = 0.0;
//
//    int p;
//
//    //constexpr bool usePanelGauss = true;
//    constexpr int panelGaussOrder = 2;   // 2x2 = 4点；3表示3x3=9点；4表示4x4=16点
//
//    const double half_draft = -ShipCfg.Geometry.Draft / 2.0;
//
//    const int half_NE = NE / 2;
//    //源点(算一半，另一半对称，前提是船左右对称，前后一半分别在同一边）
//    for (int j = 0; j < half_NE; ++j)
//    {
//        //场点
//        for (int i = 0; i < NE; ++i)
//        {
//            if(element->znr(i)>half_draft || element->znr(j) > half_draft)
//                    Gdata = green.GreenCalPanelGauss(j, i, tn, green.GF1, gGreenTable, panelGaussOrder);
//            else
//                    Gdata = green.GreenCal(j, i, tn, green.GF1, gGreenTable);
//
//            //if (SeakeepingCfg.Solver == "Potential")
//            //{
//            //    nV = -(element->A31[i] * Gdata.xdG
//            //        + element->A32[i] * Gdata.ydG
//            //        + element->A33[i] * Gdata.zdG);
//            //}
//            //else
//            //{
//            //    nV = -(element->A31[j] * Gdata.xdG
//            //        + element->A32[j] * Gdata.ydG
//            //        + element->A33[j] * Gdata.zdG);
//            //}
//
//            Green->at(tN)(j, i) = static_cast<GScalar>(Gdata.sG);
//
//            nV = element->A31[i] * Gdata.xdG
//                        + element->A32[i] * Gdata.ydG
//                        - element->A33[i] * Gdata.zdG;
//            Green_dnP->at(tN)(j, i) = static_cast<GScalar>(nV);
//
//            nV = -(element->A31[j] * Gdata.xdG
//                + element->A32[j] * Gdata.ydG
//                + element->A33[j] * Gdata.zdG);
//            Green_dnS->at(tN)(j, i) = static_cast<GScalar>(nV);
//
//            Green->at(tN)(j + half_NE, i) = Green->at(tN)(j, i);
//            Green_dnP->at(tN)(j + half_NE, i) = Green_dnP->at(tN)(j, i);
//            Green_dnS->at(tN)(j + half_NE, i) = Green_dnS->at(tN)(j, i);
//        }
//
//        //if (j == 0)
//        //{
//        //    for (int i = 0; i < 5; i++)
//        //    {
//        //        std::cout << i << ":\t" << Green_dnP->at(tN)(j, i) << ",\t";
//        //    }
//        //    std::cout << std::endl;
//        //}
//
//        for (int i = 0; i < nWL; ++i)
//        {
//            p = PotL_idx(i);
//
//            //水线重排过，要用原始的水线顺序的编号
//            int origIdx = wlOrder_.empty() ? i : wlOrder_[static_cast<std::size_t>(i)];
//            //int origIdx = i;
//
//			//std::cout << "origIdex:\t" << origIdx << "\ti:\t" << i << "\n";
//
//			//std::cout << "greencal_wL[" << i << "]:\t" << origIdx << "\n";
//            Gdata = green.GreenCal_WL_Gauss12(j, origIdx, tn, green.GF1, gGreenTable);
//
//            //if (SeakeepingCfg.Solver == "Potential")
//            //    nV = -(element->A31[i] * Gdata.xdG
//            //        + element->A32[i] * Gdata.ydG
//            //        + element->A33[i] * Gdata.zdG);
//            //else
//            //{
//            //    nV = -(element->A31[j] * Gdata.xdG
//            //        + element->A32[j] * Gdata.ydG
//            //        + element->A33[j] * Gdata.zdG);
//            //}
//
//            Gw->at(tN)(j, i) = static_cast<GScalar>(Gdata.sG);
//
//            
//			//std::cout << "p:\t" << p << "\n";
//            nV = element->A31[p] * Gdata.xdG
//                + element->A32[p] * Gdata.ydG
//                - element->A33[p] * Gdata.zdG;
//            Gw_dnP->at(tN)(j, i) = static_cast<GScalar>(nV);
//
//            nV = -(element->A31[j] * Gdata.xdG
//                + element->A32[j] * Gdata.ydG
//                + element->A33[j] * Gdata.zdG);
//            Gw_dnS->at(tN)(j, i) = static_cast<GScalar>(nV);
//
//            Gw_dx->at(tN)(j, i) = static_cast<GScalar>(Gdata.xdG);
//
//            Gw_dl->at(tN)(j, i) = static_cast<GScalar>(wlTx_(i) * Gdata.xdG + wlTy_(i) * Gdata.ydG);
//			Gw_dt->at(tN)(j, i) = static_cast<GScalar>(Gdata.tdG);
//
//            Gw->at(tN)(j + half_NE, i) = Gw->at(tN)(j, i);
//            Gw_dnP->at(tN)(j + half_NE, i) = Gw_dnP->at(tN)(j, i);
//            Gw_dnS->at(tN)(j + half_NE, i) = Gw_dnS->at(tN)(j, i);
//            Gw_dx->at(tN)(j + half_NE, i) = Gw_dx->at(tN)(j, i);
//            Gw_dl->at(tN)(j + half_NE, i) = Gw_dl->at(tN)(j, i);
//            Gw_dt->at(tN)(j + half_NE, i) = Gw_dt->at(tN)(j, i);
//        }
//    }
////}


void LinearCumminsTDGF::GreenFunction(int tN, double tn, Gsinteg& green)
{
    const int NE = SeakeepingCfg.Panel.NE;
    const int nWL = element->n_WL;

    GreenData Gdata;
    double nV = 0.0;

    int p;

    constexpr int panelGaussOrder = 2;

    const double half_draft = -ShipCfg.Geometry.Draft / 2.0;
    const int half_NE = NE / 2;

    //ship element is symmetrical on left and right sides
    int j_mirror, i_mirror;

    // 源点只算一半，另一半按左右舷对称复制
    for (int j = 0; j < NE; ++j)
    {
        j_mirror = j + half_NE;
        // 面元场点
        for (int i = 0; i < NE; ++i)
        {
            i_mirror = (i + half_NE) % NE;

            //if (element->znr(i) > half_draft || element->znr(j) > half_draft)
            if (true)
            {
                Gdata = green.GreenCalPanelGauss(
                    j,
                    i,
                    tn,
                    tdgfProvider_,
                    panelGaussOrder);
            }
            else
            {
                Gdata = green.GreenCal(
                    j,
                    i,
                    tn,
                    tdgfProvider_);
            }

            Green->at(tN)(j, i) = static_cast<GScalar>(Gdata.sG);

            nV =
                element->A31[i] * Gdata.xdG
                + element->A32[i] * Gdata.ydG
                - element->A33[i] * Gdata.zdG;

            Green_dnP->at(tN)(j, i) = static_cast<GScalar>(nV);

            nV = -(
                element->A31[j] * Gdata.xdG
                + element->A32[j] * Gdata.ydG
                + element->A33[j] * Gdata.zdG);

            Green_dnS->at(tN)(j, i) = static_cast<GScalar>(nV);

            //Green->at(tN)(j_mirror, i_mirror) = Green->at(tN)(j, i);
            //Green_dnP->at(tN)(j_mirror, i_mirror) = Green_dnP->at(tN)(j, i);
            //Green_dnS->at(tN)(j_mirror, i_mirror) = Green_dnS->at(tN)(j, i);
        }
    }


    for (int j = 0; j < NE; ++j)
    {
        // 水线项
        for (int i = 0; i < nWL; ++i)
        {
            p = PotL_idx(i);

            // 水线重排过，要用原始水线顺序编号
            int origIdx = wlOrder_.empty()
                ? i
                : wlOrder_[static_cast<std::size_t>(i)];

            Gdata = green.GreenCal_WL_Gauss12(
                j,
                origIdx,
                tn,
                tdgfProvider_);

            Gw->at(tN)(j, i) = static_cast<GScalar>(Gdata.sG);

            nV =
                element->A31[p] * Gdata.xdG
                + element->A32[p] * Gdata.ydG
                - element->A33[p] * Gdata.zdG;

            Gw_dnP->at(tN)(j, i) = static_cast<GScalar>(nV);

            nV = -(
                element->A31[j] * Gdata.xdG
                + element->A32[j] * Gdata.ydG
                + element->A33[j] * Gdata.zdG);

            Gw_dnS->at(tN)(j, i) = static_cast<GScalar>(nV);

            Gw_dx->at(tN)(j, i) = static_cast<GScalar>(Gdata.xdG);

            Gw_dl->at(tN)(j, i) =
                static_cast<GScalar>(wlTx_(i) * Gdata.xdG + wlTy_(i) * Gdata.ydG);

            Gw_dt->at(tN)(j, i) = static_cast<GScalar>(Gdata.tdG);

            //Gw->at(tN)(j_mirror, i_mirror) = Gw->at(tN)(j, i);
            //Gw_dnP->at(tN)(j_mirror, i_mirror) = Gw_dnP->at(tN)(j, i);
            //Gw_dnS->at(tN)(j_mirror, i_mirror) = Gw_dnS->at(tN)(j, i);
            //Gw_dx->at(tN)(j_mirror, i_mirror) = Gw_dx->at(tN)(j, i);
            //Gw_dl->at(tN)(j_mirror, i_mirror) = Gw_dl->at(tN)(j, i);
            //Gw_dt->at(tN)(j_mirror, i_mirror) = Gw_dt->at(tN)(j, i);
        }
    }
}


void LinearCumminsTDGF::zeroWL(int tN)
{
    const int nWL = element->n_WL;
    int p1, p2;
    for (int i = 0; i < nWL; ++i)
    {
        p1 = PotL_idx(i);
        for (int j = 0; j < nWL; ++j)
        {
            p2 = PotL_idx(j);
            for (int n = 0; n < tN; ++n)
            {
                Green->at(n)(p1, p2) = 0.0;
                Green_dnP->at(n)(p1, p2) = 0.0;
                Green_dnS->at(n)(p1, p2) = 0.0;
            }
        }
    }
}


void LinearCumminsTDGF::initialFK(fkpData& fkpdata)
{
    fkpdata.NE = SeakeepingCfg.Panel.NE;
    fkpdata.U = U;
    fkpdata.A31 = element->A31;
    fkpdata.A32 = element->A32;
    fkpdata.A33 = element->A33;
    fkpdata.xcr = element->xcr;
    fkpdata.ycr = element->ycr;
    fkpdata.zcr = element->zcr;
}


void LinearCumminsTDGF::rebuildDiffractionForceFromIntegral(
    const CaseContext& ctx)
{
    const int nImp = std::min(
        static_cast<int>(eforce.dForce.rows()),
        2 * ctx.tMot + 1);

    if (nImp <= 2)
        return;

    const double dt = ctx.dt_const;

    if (dt <= 0.0)
        throw std::runtime_error(
            "rebuildDiffractionForceFromIntegral: invalid dt.");

    constexpr int kSgWinSmooth = 31;
    constexpr int kSgWinDeriv = 41;
    constexpr int kSgPoly = 3;

    for (int k = 0; k < 6; ++k)
    {
        Eigen::VectorXd I =
            eforce.forceIntHist.col(k).head(nImp);

        //Eigen::VectorXd Is = smooth5Binomial(I, 5);
        //Eigen::VectorXd dI = differentiateCentral(Is, dt);

        const auto sg = sg::smoothThenDeriv1(I, dt, kSgWinSmooth, kSgWinDeriv, kSgPoly);
        const Eigen::VectorXd& Is = sg.ySmooth;
        const Eigen::VectorXd& dI = sg.dydt;

        eforce.forceIntSmoothHist.col(k).head(nImp) = Is;
        eforce.forceDtHist.col(k).head(nImp) = dI;
    }

    for (int n = 0; n < nImp; ++n)
    {
        // 平动和 roll
        eforce.dForce(n, 0) = -rho * eforce.forceDtHist(n, 0);
        eforce.dForce(n, 1) = -rho * eforce.forceDtHist(n, 1);
        eforce.dForce(n, 2) = -rho * eforce.forceDtHist(n, 2);
        eforce.dForce(n, 3) = -rho * eforce.forceDtHist(n, 3);

        // pitch:
        // 原 ForceCal:
        // F4 = -rho * (dI4/dt + U * ArInt.col(2).dot(Pt))
        // I2 = -ArInt.col(2).dot(Pt)
        // 所以 ArInt.col(2).dot(Pt) = -I2
        eforce.dForce(n, 4) =
            -rho * (
                eforce.forceDtHist(n, 4)
                - ctx.U * eforce.forceIntSmoothHist(n, 2)
                );

        // yaw:
        // 原 ForceCal:
        // F5 = -rho * (dI5/dt - U * ArInt.col(1).dot(Pt))
        // I1 = -ArInt.col(1).dot(Pt)
        // 所以 ArInt.col(1).dot(Pt) = -I1
        eforce.dForce(n, 5) =
            -rho * (
                eforce.forceDtHist(n, 5)
                + ctx.U * eforce.forceIntSmoothHist(n, 1)
                );
    }

    // 超出脉冲表范围的行清零，避免 writeFKImplese 读到旧值
    for (int n = nImp; n < eforce.dForce.rows(); ++n)
        eforce.dForce.row(n).setZero();
}


void LinearCumminsTDGF::buildFKImpulseDirect(CaseContext& ctx, FKphi& fkphi)
{
    ensureGreenTablesMatchOnlineGrid();

    const WaveForceRegion region =
        wave_force_region::classify(ctx.U, ctx.dirRad, ctx.W);
    const int nSteps = 2 * ctx.tMot + 1;
    std::cout << "[FKImpulse] direct compute"
        << "  Fn=" << ctx.Fn
        << "  Dir=" << ctx.dirRad
        << "  LuBand=" << wave_force_region::tag(region)
        << "  Bucket=" << wave_force_region::impulseBucketTag(region)
        << "  ω_inc=" << ctx.W
        << "  ω_e=" << ctx.we
        << "\n           dt=" << ctx.dt_const
        << "  tMot=" << ctx.tMot
        << "  memory=" << ctx.tMot * ctx.dt_const << "s"
        << "  totalSteps=" << nSteps << "\n\n";

    const double dtImp = ctx.dt_const;
    double tau = -ctx.tMot * dtImp;

    eforce.fkForce.setZero();
    eforce.dForce.setZero();
    eforce.forceIntHist.setZero();
    eforce.forceIntSmoothHist.setZero();
    eforce.forceDtHist.setZero();

    for (int tN = 0; tN <= 2 * ctx.tMot; ++tN, tau += dtImp)
    {
        ctx.wave->Exciting(tau, fkphi);
        ExcitingCal(tN, tau, fkphi, ctx);

        if (tN % 20 == 0 || tN == 2 * ctx.tMot)
        {
            std::cout << "computing exciting impulse step:\t"
                << tN << "/" << 2 * ctx.tMot << "\n";
        }
    }
    std::cout << "\n";

    rebuildDiffractionForceFromIntegral(ctx);

    FKImpulseKernelIO::Data data;
    data.fkForce = eforce.fkForce;
    data.dForce = eforce.dForce;

    FKImpulseKernelIO::save(
        makeFKImpulseParams(filePath, ShipCfg, ctx), data);

    // 这个继续作为调试输出，不作为缓存
    writeFKImplese(filePath, ctx);
}



void LinearCumminsTDGF::prepareFKImpulse(CaseContext& ctx, FKphi& fkphi)
{
    // Direct : always recompute and overwrite.
    // Load   : library only — error if no matching file.
    // Auto   : try library first, recompute on miss.
    //
    // Library lookup is fuzzy on (dt, tMot): we match (ship, Fn, Dir, Lu bucket H or 1/2/3)
    // and resample onto the current (ctx.dt_const, ctx.tMot) grid.

    const std::string method = SeakeepingCfg.FKImpulseKernel.Method;

    if (method == "Direct")
    {
        buildFKImpulseDirect(ctx, fkphi);
        return;
    }

    const WaveForceRegion region =
        wave_force_region::classify(ctx.U, ctx.dirRad, ctx.W);

    // Port/starboard symmetry: only 0–180° kernels are stored. A heading on
    // the port side (β in (π,2π)) reuses its mirror βm = 2π−β ∈ (0,π).
    // cos β is unchanged by the mirror, so the encounter region is identical;
    // only the antisymmetric DOFs flip sign — sway(1), roll(3), yaw(5) — while
    // surge(0), heave(2), pitch(4) are unchanged.
    double lookupDir = ctx.dirRad;
    while (lookupDir < 0.0)        lookupDir += 2.0 * PI;
    while (lookupDir >= 2.0 * PI)  lookupDir -= 2.0 * PI;
    bool mirrorPS = false;
    if (lookupDir > PI + 1e-9)
    {
        lookupDir = 2.0 * PI - lookupDir;
        mirrorPS = true;
    }

    const std::string libFile = FKImpulseKernelIO::findByKey(
        filePath, ShipCfg.Name, ctx.Fn, lookupDir, region);

    bool loaded = false;
    if (!libFile.empty())
    {
        FKImpulseKernelIO::LoadInfo info;
        FKImpulseKernelIO::Data     raw;
        if (FKImpulseKernelIO::loadFromFile(libFile, info, raw))
        {
            FKImpulseKernelIO::Data fit;
            FKImpulseKernelIO::resample(raw, info.dt, info.tMot,
                ctx.dt_const, ctx.tMot, fit);
            eforce.fkForce = std::move(fit.fkForce);
            eforce.dForce = std::move(fit.dForce);

            if (mirrorPS)
            {
                // Antisymmetric DOFs about the centreline: sway, roll, yaw.
                for (int c : {1, 3, 5})
                {
                    if (c < eforce.fkForce.cols()) eforce.fkForce.col(c) *= -1.0;
                    if (c < eforce.dForce.cols())  eforce.dForce.col(c) *= -1.0;
                }
            }
            loaded = true;

            std::cout << "[FKImpulse] library hit: " << libFile
                << (mirrorPS ? "  [PS-mirrored: sway/roll/yaw negated]" : "")
                << "\n"
                << "           src (dt=" << info.dt << ", tMot=" << info.tMot
                << ", memory=" << info.dt * info.tMot << "s)"
                << " -> dst (dt=" << ctx.dt_const << ", tMot=" << ctx.tMot
                << ", memory=" << ctx.tMot_const << "s)\n";
        }
    }

    if (loaded) return;

    if (method == "Load")
    {
        std::ostringstream msg;
        msg << "FK impulse cache miss for Fn=" << ctx.Fn
            << " Dir=" << ctx.dirRad
            << " Bucket=" << wave_force_region::impulseBucketTag(region)
            << " (case path " << filePath << "/fkImpulse).";
        throw std::runtime_error(msg.str());
    }

    // Auto: compute and save.
    buildFKImpulseDirect(ctx, fkphi);
}


void LinearCumminsTDGF::computeExciting(CaseContext& ctx)
{
    fkpData fkpdata;
    FKphi   fkphi(SeakeepingCfg.Panel.NE);

    initialFK(fkpdata);
    ctx.wave->loadData(fkpdata);

    // Regression sentinel (regular wave only): incidentEtaAtShip must stay
    // bit-identical to RegularWave::Eta. Skipped for irregular/cross, whose
    // ctx.wave->Eta is a stub.
    if (ctx.waveField.components().size() == 1 && ctx.reg)
    {
        double maxAbsDiff = 0.0;
        double maxAbsEta = 1.0e-30;
        const int nChk = std::min(ctx.TS, 2000);
        for (int n = 0; n < nChk; ++n)
        {
            const double t = n * ctx.dt;
            const double etaNew = incidentEtaAtShip(ctx, t);
            const double etaRef = ctx.wave->Eta(t);
            maxAbsDiff = std::max(maxAbsDiff, std::abs(etaNew - etaRef));
            maxAbsEta = std::max(maxAbsEta, std::abs(etaRef));
        }
        std::cout << "[Step2 WaveField check] max|etaNew-Eta|=" << maxAbsDiff
            << "  (rel=" << (maxAbsDiff / maxAbsEta) << ")\n";
    }

    std::filesystem::create_directories(std::filesystem::path(filePath) / "ExcitingForce");

    std::string excitingFile =
        filePath + "/ExcitingForce/ExcitingForce" + ctx.tag + ".csv";
    std::ofstream efile(excitingFile);
    if (!efile.is_open())
        throw std::runtime_error("cannot create exciting file.");

    const double tScale = std::sqrt(G / ShipCfg.Geometry.Length);

    // Optional debug dump: incident wave-elevation time history at the fixed
    // body reference point — exactly the η used in the wave-force integration
    // (encounter frame, body origin): η = Σ_s Σ_j a cos(we·t + ε + φ_spatial),
    // φ_spatial = −k(x0 cosθ + y0 sinθ). Columns: t, eta_total, eta_sub0[,sub1].
    // No region filter / no smoothing — the raw wave the ship sees.
    if (ctx.outputWaveHistory && ctx.TS > 0 && ctx.dt > 0.0
        && !ctx.subFields.empty())
    {
        std::filesystem::create_directories(
            std::filesystem::path(filePath) / "wave_history");
        std::ofstream wfile(
            filePath + "/wave_history/wave_history" + ctx.tag + ".csv");
        if (!wfile.is_open())
            throw std::runtime_error("cannot create wave_history file.");

        const int nSubW = static_cast<int>(ctx.subFields.size());
        wfile << "t,eta_total";
        for (int s = 0; s < nSubW; ++s) wfile << ",eta_sub" << s;
        wfile << "\n";

        for (int tN = 0; tN < ctx.TS; ++tN)
        {
            const double tw = tN * ctx.dt;
            double total = 0.0;
            std::vector<double> per(nSubW, 0.0);
            for (int s = 0; s < nSubW; ++s)
            {
                double e = 0.0;
                for (const auto& w : ctx.subFields[static_cast<std::size_t>(s)].components())
                {
                    double we_j = 0.0;
                    try { we_j = calcEncounterInfo(w.omega, ctx.U, w.theta).we; }
                    catch (...) { continue; }
                    const double k_j = w.omega * w.omega / G;
                    const double phiSp =
                        -k_j * (ctx.crossStartX * std::cos(w.theta)
                            + ctx.crossStartY * std::sin(w.theta));
                    e += w.a * std::cos(we_j * tw + w.eps + phiSp);
                }
                per[static_cast<std::size_t>(s)] = e;
                total += e;
            }
            wfile << tw << "," << total;
            for (int s = 0; s < nSubW; ++s) wfile << "," << per[static_cast<std::size_t>(s)];
            wfile << "\n";
        }
    }

    auto modes = SeakeepingCfg.modes;

    // Step 4: a crossing sea is the linear sum of independent sub-waves, each
    // with its own FK kernel (direction/region) and incident field. One pass
    // per sub-wave, accumulated. Regular/Irregular -> exactly one pass, so the
    // result is bit-identical to Steps 1-3.
    ctx.ExcitingForce.setZero();
    ctx.FKForceHist.setZero();
    ctx.DiffForceHist.setZero();
    const auto  sub0Wave_ = ctx.subWaves.front();
    const int   nSub_ = static_cast<int>(ctx.subWaves.size());
    for (int sIdx_ = 0; sIdx_ < nSub_; ++sIdx_)
    {
        ctx.activeSub = sIdx_;
        ctx.wave = ctx.subWaves[static_cast<std::size_t>(sIdx_)];
        ctx.dirRad = ctx.wave->direction();
        ctx.W = ctx.wave->getFreq();
        ctx.wave->loadData(fkpdata);
        const bool lastSub_ = (sIdx_ + 1 == nSub_);

        // Step 5: split this sub-wave by encounter region. Head/beam (cos β ≤ 0)
        // -> one pass, no filter (bit-identical to Steps 1-4). Following seas
        // -> one pass per non-empty F1/F2/F3 band, each with its own kernel and
        // only its components summed; the forces add (convolution linearity).
        struct RegionPass { int filter; double repW; };
        std::vector<RegionPass> regionPasses;
        {
            const double cb = std::cos(ctx.dirRad);
            const bool following = (ctx.U > 0.0 && cb > 1e-6);
            if (!following)
            {
                regionPasses.push_back({ -1, ctx.wave->getFreq() });
            }
            else
            {
                const double wc1 = G / (2.0 * ctx.U * cb);
                const double wc2 = G / (ctx.U * cb);
                double s1 = 0, s2 = 0, s3 = 0;
                int    n1 = 0, n2 = 0, n3 = 0;
                for (const auto& c :
                    ctx.subFields[static_cast<std::size_t>(sIdx_)].components())
                {
                    if (c.omega < wc1) { s1 += c.omega; ++n1; }
                    else if (c.omega < wc2) { s2 += c.omega; ++n2; }
                    else { s3 += c.omega; ++n3; }
                }
                if (n1 > 0) regionPasses.push_back(
                    { static_cast<int>(WaveForceRegion::F1), s1 / n1 });
                if (n2 > 0) regionPasses.push_back(
                    { static_cast<int>(WaveForceRegion::F2), s2 / n2 });
                if (n3 > 0) regionPasses.push_back(
                    { static_cast<int>(WaveForceRegion::F3), s3 / n3 });
                if (regionPasses.empty())
                    regionPasses.push_back({ -1, ctx.wave->getFreq() });
            }
        }

        for (std::size_t rp_ = 0; rp_ < regionPasses.size(); ++rp_)
        {
            ctx.activeRegion = regionPasses[rp_].filter;
            ctx.W = regionPasses[rp_].repW;
            const bool lastPass_ =
                lastSub_ && (rp_ + 1 == regionPasses.size());

            const WaveForceRegion impulseReg =
                wave_force_region::classify(ctx.U, ctx.dirRad, ctx.W);
            const bool impulseGridChanged =
                fkImpulsePrepTMot_ < 0
                || ctx.tMot != fkImpulsePrepTMot_
                || std::abs(ctx.dt_const - fkImpulsePrepDtConst_) > 1.0e-14;
            const bool needPrepareFK =
                ctx.new_Fn ||
                ctx.new_direction ||
                !fkImpulsePrepared_ ||
                impulseGridChanged ||
                std::abs(ctx.Fn - fkImpulsePrepFn_) > 1e-12 ||
                std::abs(ctx.dirRad - fkImpulsePrepDirRad_) > 1e-5 ||
                impulseReg != fkImpulsePrepRegion_;

            if (needPrepareFK)
            {
                prepareFKImpulse(ctx, fkphi);
                fkImpulsePrepFn_ = ctx.Fn;
                fkImpulsePrepDirRad_ = ctx.dirRad;
                fkImpulsePrepRegion_ = impulseReg;
                fkImpulsePrepDtConst_ = ctx.dt_const;
                fkImpulsePrepTMot_ = ctx.tMot;
                fkImpulsePrepared_ = true;
            }

            //if (ctx.new_direction)
            //{
            //    std::cout << "new_direction ..." << std::endl;

            //    //element->method = true;
            //    //element->RankineSource();
            //    //std::string Solver_mid = Solver;
            //    //Solver = "Source";

            //    const double dtImp = ctx.dt_const;
            //    double tau = -ctx.tMot * dtImp;

            //    for (int tN = 0; tN <= 2 * ctx.tMot; ++tN, tau += dtImp)
            //    {      
            //        ctx.wave->Exciting(tau, fkphi);
            //        ExcitingCal(tN, tau, fkphi, ctx);
            //    }

            //    // 在 writeFKImplese 前统一重建 dForce
            //    rebuildDiffractionForceFromIntegral(ctx);

            //    writeFKImplese(filePath, ctx);

            //    //element->method = false;
            //    //element->RankineSource();
            //    //Solver = Solver_mid;
            //}

            //double dt = ctx.dt;
            //double tn = dt;
            //double T = 0.2 * ctx.TS * dt;

            //auto smoothFunction = [T](double t)->double
            //{
            //    if (t < T)
            //        return 0.5 * (1 - cos(PI * t / T));
            //    else
            //        return 1;
            //};

            //for (int tN = 0; tN < ctx.TS; ++tN, tn += dt)
            //{
            //    eforce.force.setZero();

            //    switch (SeakeepingCfg.FKMethod)
            //    {
            //    case FKForceMethod::PrescribedAmp:
            //        eforce.force(2) =
            //            ctx.F3_amplitude * std::cos(ctx.we * tn) * ctx.a_scale[2];

            //        eforce.force(4) =
            //            ctx.F5_amplitude * std::cos(ctx.we * tn) * ctx.a_scale[4];
            //        break;

            //    case FKForceMethod::TDGFImpulse:
            //        ExcitingCal(tN, tn, fkphi, ctx);
            //        break;

            //    case FKForceMethod::DirectPressure:
            //        //DForceCal(tN, tn, fkphi, ctx);
            //        DForceCal(tN, tn, fkphi, ctx);
            //        computeDirectPressureFK(ctx, tn);
            //        break;
            //    }

            //    efile << tn * tScale;
            //    for (int i = 0; i < 6; ++i)
            //        efile << "," << eforce.force(i) / ctx.a_scale[i];

            //    //if (tN >= ctx.tMot)
            //    if(true)
            //    {
            //        for (int i = 0; i < ctx.DOF; ++i)
            //        {
            //            //ctx.ExcitingForce(tN - ctx.tMot, i) = eforce.force(SeakeepingCfg.modes[i]) * smoothFunction((tN - ctx.tMot) * dt);
            //            //efile << "," << ctx.ExcitingForce(tN - ctx.tMot, i) / ctx.a_scale[modes[i]];

            //            ctx.ExcitingForce(tN, i) = eforce.force(SeakeepingCfg.modes[i]) * smoothFunction((tN) * dt);
            //            efile << "," << ctx.ExcitingForce(tN, i) / ctx.a_scale[modes[i]];
            //        }
            //    }
            //    else
            //        efile << "," << 0 << "," << 0;

            //    efile << "\n";
            //}

            const double dtWave = ctx.dt;
            const double Tfade = 0.1 * ctx.TS * dtWave;

            auto smoothFunction = [Tfade](double t) -> double
                {
                    if (t < Tfade)
                        return 0.5 * (1.0 - std::cos(PI * t / Tfade));
                    return 1.0;
                };

            for (int tN = 0; tN < ctx.TS; ++tN)
            {
                const double tNow = tN * dtWave;

                eforce.force.setZero();
                Eigen::RowVectorXd fkPart = Eigen::RowVectorXd::Zero(6);
                Eigen::RowVectorXd dfPart = Eigen::RowVectorXd::Zero(6);

                switch (SeakeepingCfg.FKMethod)
                {
                case FKForceMethod::PrescribedAmp:
                    eforce.force(2) =
                        ctx.F3_amplitude * std::cos(ctx.we * tNow) * ctx.a_scale[2];

                    eforce.force(4) =
                        ctx.F5_amplitude * std::cos(ctx.we * tNow) * ctx.a_scale[4];
                    break;

                case FKForceMethod::TDGFImpulse:
                    eforce.force.setZero();
                    fkDfCal(tN, tNow, fkphi, ctx, /*dfOnly=*/true);
                    dfPart = eforce.force;
                    eforce.force.setZero();
                    fkDfCal(tN, tNow, fkphi, ctx, /*dfOnly=*/false, /*fkOnly=*/true);
                    fkPart = eforce.force;
                    eforce.force = fkPart + dfPart;
                    break;

                case FKForceMethod::DirectPressure:
                    eforce.force.setZero();
                    fkDfCal(tN, tNow, fkphi, ctx, /*dfOnly=*/true);
                    dfPart = eforce.force;
                    eforce.force.setZero();
                    directFkOverComponents(ctx, tNow);
                    fkPart = eforce.force;
                    eforce.force = fkPart + dfPart;
                    break;
                }

                const double fade = smoothFunction(tNow);
                for (int i = 0; i < ctx.DOF; ++i)
                {
                    const int mode = SeakeepingCfg.modes[i];
                    ctx.ExcitingForce(tN, i) += eforce.force(mode) * fade;
                    ctx.FKForceHist(tN, i) += fkPart(mode) * fade;
                    ctx.DiffForceHist(tN, i) += dfPart(mode) * fade;
                }

                if (lastPass_)
                {
                    efile << tNow;
                    for (int i = 0; i < 6; ++i)
                        efile << "," << eforce.force(i) / ctx.a_scale[i];
                    for (int i = 0; i < ctx.DOF; ++i)
                    {
                        const int mode = SeakeepingCfg.modes[i];
                        efile << "," << ctx.ExcitingForce(tN, i) / ctx.a_scale[mode];
                    }
                    efile << "\n";
                }
            }

        } // end per-region pass
    } // end per-sub-wave loop

    // Restore ctx to sub-wave 0 for the downstream solve.
    ctx.activeSub = 0;
    ctx.activeRegion = -1;
    ctx.wave = sub0Wave_;
    ctx.dirRad = ctx.wave->direction();
    ctx.W = ctx.wave->getFreq();


    //element->method = true;
    //element->RankineSource();
    //std::string Solver_mid = Solver;
    //Solver = "Source";

    //element->saveRankine(filePath + "Rankine1.csv");


    /*for (int tN = 0; tN < ctx.TS; ++tN, tn += dt, tnM += dt)
    {
        eforce.force.setZero();

        switch (SeakeepingCfg.FKMethod)
        {
        case FKForceMethod::PrescribedAmp:
            eforce.force(2) =
                ctx.F3_amplitude * std::cos(ctx.we * tnM) * ctx.a_scale[2];

            eforce.force(4) =
                ctx.F5_amplitude * std::cos(ctx.we * tnM) * ctx.a_scale[4];
            break;

        case FKForceMethod::TDGFImpulse:
            if (ctx.new_direction)
                ctx.wave->Exciting(tnM, fkphi);

            ExcitingCal(tN, tn, fkphi, ctx);
            break;

        case FKForceMethod::DirectPressure:

            if (ctx.new_direction)
                ctx.wave->Exciting(tnM, fkphi);

            DForceCal(tN, tn, fkphi, ctx);
            computeDirectPressureFK(ctx, tnM);

            break;
        }

        efile << tnM * tScale;
        for (int i = 0; i < 6; ++i)
            efile << "," << eforce.force(i) / ctx.a_scale[i];


        if (tN >= ctx.tMot)
        {
            for (int i = 0; i < ctx.DOF; ++i)
            {
                ctx.ExcitingForce(tN - ctx.tMot, i) = eforce.force(SeakeepingCfg.modes[i]) * smoothFunction((tN - ctx.tMot) * dt);

                efile << "," << ctx.ExcitingForce(tN - ctx.tMot, i) / ctx.a_scale[modes[i]];
            }
        }
        else
            efile << "," << 0 << "," << 0;

        efile << "\n";
    }*/

    //element->method = false;
    //element->RankineSource();
    //element->saveRankine(filePath + "Rankine2.csv");
    //Solver = Solver_mid;


    //saveExcitationKernelCache(ctx);
    // 一阶波浪力：对 ExcitingForce 用 |ω_e| 做 sin/cos+常数项最小二乘（末段若干周期），
    // 便于与 Lu et al. (2024) 图 11/12 等力幅曲线对照（与运动 RAO 拟合不是同一量）。
    writeFirstOrderWaveForceAmpPhaseTable(
        filePath, ctx, SeakeepingCfg, ShipCfg.Geometry.Length);
    std::cout << "compute exciting done. " << "\n";
}

void LinearCumminsTDGF::ExcitingCal(int tN, double tn, FKphi& fkphi, const CaseContext& ctx)
{
    double temp = 0.0;
    auto& wave = ctx.wave;

    //int t0 = std::max(0, tN - 2 * ctx.tMot);


    //if (t0 == 0 && ctx.new_direction)
    if (true)
    {
        if (!(element->method))
        {
            faiConvolution(
                tN, ctx,
                fkphi.df,
                eforce.veVn_d, eforce.veVn_g,
                eforce.vePhi_d, eforce.vePhi_g,
                eforce.sPot);

            ForceCal(tN, eforce.sPot, eforce.eFdt, eforce.force, ctx.dt_const);
        }
        else
        {
            SourceConvolution(tN, fkphi.df, eforce.veSg_d, eforce.veSg_g, eforce.sPot);
            ForceCal(tN, eforce.sPot, eforce.eFdt, eforce.force, ctx.dt_const);
        }
        eforce.dForce.row(tN) = eforce.force;

        for (int j = 0; j < 6; ++j)
            eforce.fkForce(tN, j) = -element->ArInt.col(j).dot(fkphi.fk) * rho;
    }

    //const double dt = SeakeepingCfg.Time.dt;
    //const double dt = ctx.dt_const;

    //for (int j = 0; j < 6; ++j)
    //{
    //    double ft1 = 0.0;
    //  
    //    for (int tn1 = t0; tn1 < tN; ++tn1)
    //    {
    //        ft1 += (eforce.dForce(tN - tn1, j) + eforce.fkForce(tN - tn1, j))
    //            * wave->Eta(tn1 * dt) * dt;
    //    }
    //    eforce.force(j) = ft1;
    //}
}

/*
void LinearCumminsTDGF::DForceCal(int tN, double tn, FKphi& fkphi, const CaseContext& ctx)
{
    //double temp = 0.0;
    auto& wave = ctx.wave;

    //int t0 = std::max(0, tN - 2 * ctx.tMot);


    //if (t0 == 0 && ctx.new_direction)
    //{
    //    if (!(element->method))
    //    {
    //        faiConvolution(
    //            tN,
    //            fkphi.df,
    //            eforce.veVn_d, eforce.veVn_g,
    //            eforce.vePhi_d, eforce.vePhi_g,
    //            eforce.sPot);

    //        ForceCal(tN, eforce.sPot, eforce.eFdt, eforce.force);
    //    }
    //    else
    //    {
    //        SourceConvolution(tN, fkphi.df, eforce.veSg_d, eforce.veSg_g, eforce.sPot);
    //        ForceCal(tN, eforce.sPot, eforce.eFdt, eforce.force);
    //    }
    //    eforce.dForce.row(tN) = eforce.force;
    //}

    //const double dt = SeakeepingCfg.Time.dt;

    const double dt = ctx.dt;
    const double tMot_const = ctx.tMot_const;
    const double t_total = tN * dt;

    double t = t_total - tMot_const;

    auto dForce = [](double t, int j)-> double
    {

    }

    for (int j = 0; j < 6; ++j)
    {
        double ft1 = 0.0;

        //for (int tn1 = t0; tn1 < tN; ++tn1)
        while (t < t_total + tMot_const)
        {
            ft1 += eforce.dForce(t-t_total, j) * wave->Eta(t) * dt;
        }
        eforce.force(j) = ft1;
    }
}
*/



//void LinearCumminsTDGF::DForceCal(int tN, double tn, FKphi& fkphi, const CaseContext& ctx)
//{
//    auto& wave = ctx.wave;
//
//    const double dtWave = ctx.dt;
//    const double dtImp = ctx.dt_const;
//
//    const int nRow = static_cast<int>(eforce.dForce.rows());
//    const int nCol = static_cast<int>(eforce.dForce.cols());
//
//    const int iCenter = ctx.tMot;
//    const int iFirst = 0;
//    const int iLast = std::min(nRow - 1, 2 * iCenter);
//
//    const double tauMin = (iFirst - iCenter) * dtImp;
//    const double tauMax = (iLast - iCenter) * dtImp;
//
//    auto dForce = [&](double tau, int j) -> double
//    {
//        if (j < 0 || j >= nCol)
//            throw std::out_of_range("dForce: DOF index out of range.");
//
//        if (nRow <= 1 || dtImp <= 0.0)
//            throw std::runtime_error("dForce: invalid impulse table or dtImp.");
//
//        const double eps = 1.0e-12;
//
//        if (tau < tauMin - eps || tau > tauMax + eps)
//            throw std::out_of_range("dForce: time is outside impulse table range.");
//
//        if (tau <= tauMin) return eforce.dForce(iFirst, j);
//        if (tau >= tauMax) return eforce.dForce(iLast, j);
//
//        const double x = tau / dtImp + static_cast<double>(iCenter);
//
//        const int i0 = static_cast<int>(std::floor(x));
//        const int i1 = i0 + 1;
//
//        const double a = x - static_cast<double>(i0);
//
//        return (1.0 - a) * eforce.dForce(i0, j)
//            + a * eforce.dForce(i1, j);
//    };
//
//    // 当前有量纲时刻。建议这里用传进来的 tn，而不是 tN * ctx.dt。
//    const double tNow = tn;
//
//    for (int j = 0; j < 6; ++j)
//    {
//        double ft1 = 0.0;
//
//        for (double tau = tauMin; tau <= tauMax + 0.5 * dtWave; tau += dtWave)
//        //for (double tau = tauMin; tau <= 0; tau += dtWave)
//        {
//            if (tau > tauMax) tau = tauMax;
//
//            const double tWave = tNow + tau;
//
//            ft1 += dForce(tau, j) * wave->Eta(tWave) * dtWave;
//
//            if (tau >= tauMax)
//                break;
//        }
//
//        eforce.force(j) = ft1;
//    }
//}


// DForceCal removed: merged into the DirectPressure path
// (fkDfCal(dfOnly) + directFkOverComponents), verified bit-identical.


double LinearCumminsTDGF::incidentEtaAtShip(const CaseContext& ctx, double t) const
{
    // Encounter-frame incident elevation seen by the moving ship:
    //   η(t) = Σ_j a_j · cos(we_j · t + ε_j),
    //   we_j = calcEncounterInfo(ω_j, U, θ_j).we
    // For a single RegularWave component this is bit-identical to
    // RegularWave::Eta (verified in Step 1: rel == 0). Generalises directly to
    // irregular / cross seas — same kernel groups, just more components summed.
    const WaveField& fld =
        (!ctx.subFields.empty()
            && ctx.activeSub >= 0
            && ctx.activeSub < static_cast<int>(ctx.subFields.size()))
        ? ctx.subFields[static_cast<std::size_t>(ctx.activeSub)]
        : ctx.waveField;

    double eta = 0.0;
    for (const auto& w : fld.components())
    {
        // Step 5: in following seas only sum the components belonging to the
        // region currently being convolved (head/beam: activeRegion == -1 ->
        // no filter, bit-identical to Steps 1-4).
        if (ctx.activeRegion >= 0)
        {
            const int r = static_cast<int>(
                wave_force_region::classify(ctx.U, w.theta, w.omega));
            if (r != ctx.activeRegion) continue;
        }
        const double we_j = calcEncounterInfo(w.omega, ctx.U, w.theta).we;
        // Crossing-sea start-position phase: −k(x0 cosθ + y0 sinθ), k = ω²/g.
        // crossStart* is 0 unless this is a CrossWave -> single waves unchanged.
        const double k_j = w.omega * w.omega / G;
        const double phiSpatial =
            -k_j * (ctx.crossStartX * std::cos(w.theta)
                + ctx.crossStartY * std::sin(w.theta));
        eta += w.a * std::cos(we_j * t + w.eps + phiSpatial);
    }
    return eta;
}

void LinearCumminsTDGF::fkDfCal(
    int tN,
    double tn,
    FKphi& fkphi,
    const CaseContext& ctx,
    bool dfOnly,
    bool fkOnly)
{
    if (dfOnly && fkOnly)
        throw std::invalid_argument("fkDfCal: dfOnly and fkOnly are mutually exclusive.");

    auto& wave = ctx.wave;

    const double dtWave = ctx.dt;
    const double dtImp = ctx.dt_const;

    const int nRow = static_cast<int>(eforce.dForce.rows());
    const int nCol = static_cast<int>(eforce.dForce.cols());

    const int iCenter = ctx.tMot;

    const int iFirst = 0;
    const int iLast = std::min(nRow - 1, 2 * iCenter);

    const double tauMin = (iFirst - iCenter) * dtImp;
    const double tauMax = (iLast - iCenter) * dtImp;

    auto fkDfForceAt = [&](double tau, int j) -> double
        {
            if (j < 0 || j >= nCol)
                throw std::out_of_range("fkDfForceAt: DOF index out of range.");

            if (nRow <= 1 || dtImp <= 0.0)
                throw std::runtime_error("fkDfForceAt: invalid impulse table.");

            const double eps =
                1.0e-12 * std::max(1.0, std::abs(tauMax - tauMin));

            if (tau < tauMin - eps || tau > tauMax + eps)
                return 0.0;

            auto kern = [&](int i) -> double
                {
                    if (dfOnly)
                        return eforce.dForce(i, j);
                    if (fkOnly)
                        return eforce.fkForce(i, j);
                    return eforce.dForce(i, j) + eforce.fkForce(i, j);
                };

            if (tau <= tauMin) return kern(iFirst);
            if (tau >= tauMax) return kern(iLast);

            const double x = tau / dtImp + static_cast<double>(iCenter);

            const int i0 = static_cast<int>(std::floor(x));
            const int i1 = i0 + 1;

            const double a = x - static_cast<double>(i0);

            return (1.0 - a) * kern(i0) + a * kern(i1);
        };

    const double tNow = tn;

    for (int j = 0; j < 6; ++j)
    {
        double ft = 0.0;

        // 先按双边核写法，与你现在的 impulse table 对齐
        for (double tau = tauMin; tau <= tauMax + 0.5 * dtWave; tau += dtWave)
        {
            if (tau > tauMax)
                tau = tauMax;

            const double tWave = tNow - tau;

            ft += fkDfForceAt(tau, j) * incidentEtaAtShip(ctx, tWave) * dtWave;

            if (tau >= tauMax)
                break;
        }

        eforce.force(j) = ft;
    }
}


// computeDirectPressureFK removed: the scalar single-component FK was merged
// into directFkOverComponents (component-summed), verified bit-identical for
// regular and crossing seas. Start-position spatial phase and per-sub-wave
// phase0 now live in directFkOverComponents.

void LinearCumminsTDGF::directFkOverComponents(const CaseContext& ctx, double tNow)
{
    const WaveField& fld =
        (!ctx.subFields.empty()
            && ctx.activeSub >= 0
            && ctx.activeSub < static_cast<int>(ctx.subFields.size()))
        ? ctx.subFields[static_cast<std::size_t>(ctx.activeSub)]
        : ctx.waveField;

    for (const auto& w : fld.components())
    {
        // Same region filter as incidentEtaAtShip (head/beam: activeRegion=-1).
        if (ctx.activeRegion >= 0)
        {
            const int r = static_cast<int>(
                wave_force_region::classify(ctx.U, w.theta, w.omega));
            if (r != ctx.activeRegion) continue;
        }

        double we = 0.0;
        try { we = calcEncounterInfo(w.omega, ctx.U, w.theta).we; }
        catch (...) { continue; }   // following-seas degeneracy: skip

        // Pure seakeeping = straight transit: Φ0 at the body origin is the
        // encounter-frame phase we_j·t + ε_j (consistent with incidentEtaAtShip
        // Σ a_j cos(we_j t + ε_j)). β = θ (ship axis is the reference).
        // Crossing-sea start-position phase: −k(x0 cosθ + y0 sinθ), k = ω²/g.
        // crossStart* is 0 unless this is a CrossWave -> single waves unchanged.
        const double k = w.omega * w.omega / G;
        const double phiSpatial =
            -k * (ctx.crossStartX * std::cos(w.theta)
                + ctx.crossStartY * std::sin(w.theta));

        DirectPressureFKContext c;
        c.amp = w.a;
        c.omega = w.omega;
        c.beta = w.theta;
        c.phiOrigin = we * tNow + w.eps + phiSpatial;

        eforce.force += DirectPressureFK::computeForce6(*element, c, tNow);
    }
}



void LinearCumminsTDGF::writeFKImplese(const std::string filePath, const CaseContext& ctx)
{
    const WaveForceRegion region =
        wave_force_region::classify(ctx.U, ctx.dirRad, ctx.W);

    const std::filesystem::path dir = std::filesystem::path(filePath) / "fkImpulse";
    std::filesystem::create_directories(dir);

    std::ostringstream fname;
    fname << "fkImpulse_nd_";
    if (SeakeepingCfg.FKImpulseKernel.IncludeShipName && !ShipCfg.Name.empty())
        fname << ShipCfg.Name << "_";
    fname << "Fn" << FKImpulseKernelIO::keyDouble(ctx.Fn)
        << "_U" << FKImpulseKernelIO::keyDouble(ctx.U)
        << "_Dir" << FKImpulseKernelIO::keyDouble(ctx.dirRad)
        << "_Bkt" << wave_force_region::impulseBucketTag(region)
        << ".csv";

    const std::string outFile = (dir / fname.str()).string();
    std::ofstream out(outFile);
    if (!out.is_open())
    {
        std::cerr << "writeFKImplese: cannot open " << outFile << "\n";
        return;
    }

    const double L = ShipCfg.Geometry.Length;
    const double displacement = ShipCfg.Geometry.Displacement;

    const double t_scale = std::sqrt(G / L);

    // Split FK / diffraction columns use the project-specific scales below; only the
    // *_plus_* sum columns use the literature-style denominators (K3sum, K5sum).
    const double scaleK3 = rho * G * displacement / L * std::sqrt(G / L); // K30, K37
    const double scaleK5 = rho * G * displacement * std::sqrt(G / L);     // K50, K57

    const double scaleK3sum = rho * G * L * L * t_scale;
    const double scaleK5sum = rho * G * L * L * L * t_scale;

    //const double scaleK3 = rho * G * displacement * 2.0 * ctx.Amp / L;
    //const double scaleK5 = rho * G * displacement * 2.0 * ctx.Amp;

    //const double scaleK3sum = scaleK3;
    //const double scaleK5sum = scaleK5;

    out << std::setprecision(17);
    out << "t_nd,K30,K37,K50,K57,K30_plus_K37,K50_plus_K57,"
        << "fai72, fai74, fai72_smooth, fai74_smooth, fai72_dt, fai74_dt\n";

    double dt = ctx.dt_const;
    for (int i = 0; i < eforce.fkForce.rows(); ++i)
    {
        const double tnd = (i - ctx.tMot) * dt * t_scale;

        const double K30 = eforce.fkForce(i, 2) / scaleK3;
        const double K37 = eforce.dForce(i, 2) / scaleK3;

        const double K50 = eforce.fkForce(i, 4) / scaleK5;
        const double K57 = eforce.dForce(i, 4) / scaleK5;

        const double K3sum = (eforce.fkForce(i, 2) + eforce.dForce(i, 2)) / scaleK3sum;
        const double K5sum = (eforce.fkForce(i, 4) + eforce.dForce(i, 4)) / scaleK5sum;

        out << tnd << ","
            << K30 << "," << K37 << ","
            << K50 << "," << K57 << ","
            << K3sum << "," << K5sum << ","
            << eforce.forceIntHist(i, 2) << "," << eforce.forceIntHist(i, 4) << ","
            << eforce.forceIntSmoothHist(i, 2) << "," << eforce.forceIntSmoothHist(i, 4) << ","
            << eforce.forceDtHist(i, 2) << "," << eforce.forceDtHist(i, 4) << "," << std::endl;
    }
    std::cout << "write FK impulse (nondim) done: " << outFile << std::endl;
}



std::string LinearCumminsTDGF::makeExcitationKernelCachePath(const CaseContext& ctx) const
{
    std::filesystem::path dir = std::filesystem::path(filePath) / "kernel_cache" / "excitation";

    auto keyDouble = [](double x)
        {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(4) << x;
            std::string s = ss.str();
            for (char& c : s)
            {
                if (c == '.') c = 'p';
                if (c == '-') c = 'm';
            }
            return s;
        };

    std::ostringstream modes;
    for (std::size_t i = 0; i < SeakeepingCfg.modes.size(); ++i)
    {
        if (i) modes << "_";
        modes << SeakeepingCfg.modes[i];
    }

    std::ostringstream name;
    name << "tdexc_"
        << ShipCfg.Name
        << "_Fn" << keyDouble(ctx.Fn)
        << "_Dir" << keyDouble(ctx.dirRad)
        << "_W" << keyDouble(ctx.W)
        << "_We" << keyDouble(ctx.we)
        << "_modes_" << modes.str()
        << ".bin";

    return (dir / name.str()).string();
}

void LinearCumminsTDGF::saveExcitationKernelCache(const CaseContext& ctx) const
{
    TDGFExcitationKernelData data;
    data.Fn = ctx.Fn;
    data.U = ctx.U;
    data.betaRel = ctx.dirRad;
    data.omegaIncident = ctx.W;
    data.omegaEncounter = ctx.we;
    data.dt = ctx.dt;
    data.DOF = ctx.DOF;
    data.modes = SeakeepingCfg.modes;
    data.region = wave_force_region::classify(ctx.U, ctx.dirRad, ctx.W);
    data.Qlag = ctx.ExcitingForce;

    const std::string file = makeExcitationKernelCachePath(ctx);
    if (!TDGFExcitationKernelIO::save(file, data))
        throw std::runtime_error("LinearCumminsTDGF: failed to save excitation kernel cache.");
}


Eigen::VectorXd LinearCumminsTDGF::buildPotentialRhs(
    const Eigen::VectorXd& boundary_condition) const
{
    return element->Rz * boundary_condition;
}

void LinearCumminsTDGF::faiConvolution(
    const int tN, const CaseContext& ctx,
    const Eigen::VectorXd& Vn,
    Eigen::MatrixXd& VnHist_d,
    SgMatG& VnHist_g,
    Eigen::MatrixXd& PhiHist_d,
    SgMatG& PhiHist_g,
    Eigen::VectorXd& phiBody)
{
    // potential/direct method:
    // 未知量:  phi_D(., tN)
    // 已知量:  Vn = \partial phi_D / \partial n |_{tN}
    //
    // 1)当前时刻瞬时项
    // 2) 面元历史卷积项
    // 3) 水线 bracket 项（~ U^2/g）
    // 4) 水线差分项（~ 2U/g），其中 n = N-1 会引入当前未知 phi_N，
    //    因此需要把这一部分移到左端，用 fixed-point 迭代解。

    const auto& dGz = Green_dnP;   // potential 分支：场点法向导数核
    const auto& Gz = Green;

    const double dt = ctx.dt_const;

    const int NE = SeakeepingCfg.Panel.NE;


    //const int t0 = (tN > SeakeepingCfg.Time.GreenStep)
    //    ? (tN - SeakeepingCfg.Time.GreenStep) : 0;

    const int t0 = std::max(0, tN - ctx.TG);

    //const double dt = SeakeepingCfg.Time.dt;

    // ---------- 当前时刻瞬时项 ----------
    // 这里默认 element->Rz / element->lu 对应势函数法当前时刻的 1/r 离散系统
    Eigen::VectorXd rhs = buildPotentialRhs(Vn);

    // tN = 0 时没有历史项，直接求当前 phi 即可
    if (tN == 0)
    {
        phiBody = solveSourceLinearSystem(rhs);

        VnHist_d.col(col(tN)) = Vn;
        VnHist_g.col(col(tN)) = Vn.cast<GScalar>();
        PhiHist_d.col(col(tN)) = phiBody;
        PhiHist_g.col(col(tN)) = phiBody.cast<GScalar>();
        return;
    }

    // ---------- 面元历史卷积 + 水线 bracket 历史项 ----------
    GVec vn_g(NE);
    GVec phi_g(NE);

    for (int tn = t0; tn < tN; ++tn)
    {
        const int lag = tN - tn - 1;

        vn_g.noalias() = VnHist_g.col(col(tn));
        phi_g.noalias() = PhiHist_g.col(col(tn));

        // 面元历史项:
        // + dt * \int G * (dphi/dn)
        // - dt * \int dG/dn * phi
        rhs.noalias() += dt * (Gz->at(lag) * vn_g).cast<double>();
        rhs.noalias() -= dt * (dGz->at(lag) * phi_g).cast<double>();

        // (U^2/g) dt * ∮ n1^2 G * Vn dl


        rhs.noalias() -= (UsquareG * dt) * applyWaterlineN1SqG(
            lag, VnHist_d.col(col(tn)));

        // 水线 bracket 项:
        // - (U^2/g) dt * \oint phi * [ ... ] dl


        rhs.noalias() += (UsquareG * dt) * applyWaterlineBracket1(
            lag, PhiHist_d.col(col(tn)));
    }

    // ---------- 水线差分项 ----------
    // 图中最后一项:
    // - (2U/g) * sum_n \oint n1 [phi_{n+1} - phi_n] G dl
    //
    // 其中 n = tN - 1 这一项包含当前未知 phi_N:
    // - (2U/g) * W0 * phi_N + (2U/g) * W0 * phi_{N-1}
    //
    // 把含 phi_N 的第一部分移到左边，于是右端先只保留显式部分。

    if (std::abs(U) > 1e-14)
        //if(false)
    {
        // 显式历史差分项: n = t0, ..., tN-2
        for (int tn = t0; tn <= tN - 2; ++tn)
        {
            const int lag = tN - tn - 1;

            rhs.noalias() += (2.0 * U / G) * applyWaterlineG(
                lag,
                PhiHist_d.col(col(tn + 1)) - PhiHist_d.col(col(tn)));

            //rhs.noalias() += (2.0 * U / G * dtShared_) * applyWaterlineG(
            //lag,
            //PhiHist_d.col(col(tn)));
        }

        // n = tN-1 拆项后留下的显式部分: +(2U/g) * W0 * phi_{N-1}
        rhs.noalias() -= (2.0 * U / G) * applyWaterlineG(
            0, PhiHist_d.col(col(tN - 1)));

        //rhs.noalias() += (2.0 * U / G * dtShared_) * applyWaterlineG(
        //    0, PhiHist_d.col(col(tN - 1)));

        // ---------- 解 (A + M) phi_N = rhs ----------
        // A 由 element->lu 已分解；M = (2U/g) * W0
        // 固定点: A x^{m+1} = rhs - M x^{m}
        Eigen::VectorXd x = solveSourceLinearSystem(rhs);

        for (int iter = 0; iter < 1000; ++iter)
        {
            const Eigen::VectorXd corr =
                (2.0 * U / G) * applyWaterlineG(0, x);

            const Eigen::VectorXd xNew =
                solveSourceLinearSystem(rhs + corr);

            const double denom = std::max(1.0, xNew.norm());
            if ((xNew - x).norm() / denom < 1e-10)
            {
                phiBody = xNew;

                VnHist_d.col(col(tN)) = Vn;
                VnHist_g.col(col(tN)) = Vn.cast<GScalar>();
                PhiHist_d.col(col(tN)) = phiBody;
                PhiHist_g.col(col(tN)) = phiBody.cast<GScalar>();
                return;
            }
            x = 0.6 * x + 0.4 * xNew;
        }

        // 迭代没提前收敛时，返回最后一次结果
        phiBody = x;
    }
    else
    {
        // U = 0 时没有最后这一类水线差分项
        phiBody = solveSourceLinearSystem(rhs);
    }

    // ---------- 4. 存历史 ----------
    VnHist_d.col(col(tN)) = Vn;
    VnHist_g.col(col(tN)) = Vn.cast<GScalar>();
    PhiHist_d.col(col(tN)) = phiBody;
    PhiHist_g.col(col(tN)) = phiBody.cast<GScalar>();
}


void LinearCumminsTDGF::SourceConvolution(
    int& tN, const Eigen::VectorXd& Vn,
    Eigen::MatrixXd& Sg_d, SgMatG& Sg_g,
    Eigen::VectorXd& sPot)
{
    auto& dGz_wlT = (Solver == "Source") ? Gw_dnS : Gw_dnP;
    auto& dGz = (Solver == "Source") ? Green_dnS : Green_dnP;

    const int NE = SeakeepingCfg.Panel.NE;
    const int t0 = (tN > SeakeepingCfg.Time.GreenStep) ? (tN - SeakeepingCfg.Time.GreenStep) : 0;
    const double dt = SeakeepingCfg.Time.dt;

    ak = 4 * PI * Vn;

    const int nth = omp_get_max_threads();
    std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>> partial(nth);
    for (auto& v : partial) v = Eigen::VectorXd::Zero(NE);

#pragma omp parallel
    {
        shipsim::EigenSingleThreadGuard _g;
        const int tid = omp_get_thread_num();
        Eigen::VectorXd& acc = partial[tid];

        GVec sg_g(NE);
        GVec wl_weight_g(element->n_WL);
        GVec tmpNE1(NE), tmpNE2(NE);

#pragma omp for schedule(static)
        for (int tn = t0; tn < tN; ++tn)
        {
            const int lag = tN - tn - 1;

            const auto& dGz_wlT_m = dGz_wlT->at(lag);
            const auto& dGz_mat = dGz->at(lag);

            sg_g.noalias() = Sg_g.col(col(tn));
            wl_weight_g = sg_g(PotL_idx).cwiseProduct(N0sq);

            tmpNE1.noalias() = (dGz_wlT_m * wl_weight_g);
            tmpNE2.noalias() = (dGz_mat * sg_g);

            acc.noalias() += dt * UsquareG * tmpNE1.cast<double>();
            acc.noalias() -= dt * tmpNE2.cast<double>();
        }
    }

    for (auto& v : partial) ak += v;

    Eigen::VectorXd x = element->lu.solve(ak);

    Sg_d.col(col(tN)) = x;
    Sg_g.col(col(tN)) = x.cast<GScalar>();

    SourcePotential(t0, tN, sPot, Sg_d, Sg_g);
}

void LinearCumminsTDGF::SourcePotential(
    int t0, int tN, Eigen::VectorXd& Pt,
    const Eigen::MatrixXd& Sg_d, const SgMatG& Sg_g)
{
    Pt = element->Rz * Sg_d.col(col(tN));

    const double dt = SeakeepingCfg.Time.dt;

    GVec sg_g(SeakeepingCfg.Panel.NE);
    GVec wl_weight_g(element->n_WL);

    for (int tn = t0; tn < tN; ++tn)
    {
        const int lag = tN - tn - 1;

        const auto& Gz_mat = Green->at(lag);
        const auto& Gz_wlT_m = Gw->at(lag);

        sg_g.noalias() = Sg_g.col(col(tn));
        wl_weight_g = sg_g(PotL_idx).cwiseProduct(N0sq);

        Pt.noalias() += dt * (Gz_mat * sg_g).cast<double>();
        Pt.noalias() -= dt * UsquareG * (Gz_wlT_m * wl_weight_g).cast<double>();
    }

    Pt /= (4 * PI);
}

//void LinearCumminsTDGF::RadiationCal(int tN, double /*tn*/, const Eigen::VectorXd& rVnLocal)
//{
// 
//    if (SeakeepingCfg.Solver == "Potential")
//    {
//        SourceConvolution(tN, rVnLocal, rforce.vSg_d, rforce.vSg_g, rforce.sPot);
//        ForceCal(tN, rforce.sPot, rforce.rFdt, rforce.force);
//    }
//    else
//    {
//        SourceConvolution(tN, rVnLocal, rforce.vSg_d, rforce.vSg_g, rforce.sPot);
//        ForceCal(tN, rforce.sPot, rforce.rFdt, rforce.force);
//    }
//}

//void LinearCumminsTDGF::ForceCal(int& /*tN*/, Eigen::VectorXd& Pt,
//    Eigen::RowVectorXd& Fdt, Eigen::RowVectorXd& F, const double dt)
//{
//    for (int k = 0; k < 6; ++k)
//    {
//        const double force = -element->ArInt.col(k).dot(Pt);
//
//        //const double derivativeofforce = (force - Fdt(k)) / SeakeepingCfg.Time.dt;
//
//        const double derivativeofforce = (force - Fdt(k)) / dt;
//
//        Fdt(k) = force;
//
//        switch (k)
//        {
//        case 0:
//        case 1:
//        case 2:
//        case 3:
//            F(k) = -derivativeofforce * rho;
//            break;
//        case 4:
//            F(k) = -(derivativeofforce + U * element->ArInt.col(2).dot(Pt)) * rho;
//            break;
//        case 5:
//            F(k) = -(derivativeofforce - U * element->ArInt.col(1).dot(Pt)) * rho;
//            break;
//        }
//    }
//}


void LinearCumminsTDGF::ForceCal(
    int& tN,
    Eigen::VectorXd& Pt,
    Eigen::RowVectorXd& Fdt,
    Eigen::RowVectorXd& F,
    const double dt)
{
    (void)Fdt;
    (void)dt;

    if (tN < 0 || tN >= eforce.forceIntHist.rows())
        throw std::out_of_range("ForceCal: tN out of forceIntHist range.");

    for (int k = 0; k < 6; ++k)
    {
        // I_k(t) = -∫ phi n_k dS
        eforce.forceIntHist(tN, k) =
            -element->ArInt.col(k).dot(Pt);
    }

    // 这里先不算绕射力，后面统一平滑求导后再覆盖 eforce.dForce
    //F.setZero();
}



// -----------------------------------------------------------------------------
// 论文 3.3 / 3.4：刚体船舶运动
// -----------------------------------------------------------------------------

Eigen::VectorXd LinearCumminsTDGF::buildRankineRhs(const Eigen::VectorXd& boundary_condition) const
{
    return element->Rz * boundary_condition;
}

Eigen::VectorXd LinearCumminsTDGF::solveRankineBIE(
    const Eigen::VectorXd& boundary_condition) const
{
    const Eigen::VectorXd rhs = buildRankineRhs(boundary_condition);
    return element->lu.solve(rhs);
}

Eigen::VectorXd LinearCumminsTDGF::solveSourceLinearSystem(const Eigen::VectorXd& rhs) const
{
    return element->lu.solve(rhs);
}

Eigen::VectorXd LinearCumminsTDGF::buildBoundaryCondition_a(int mode) const
{
    return -element->Nvec.col(mode);
}

Eigen::VectorXd LinearCumminsTDGF::buildBoundaryCondition_b(int mode) const
{
    const int NE = SeakeepingCfg.Panel.NE;
    Eigen::VectorXd b = Eigen::VectorXd::Zero(NE);

    // 对刚体 6DOF：论文 3.22 的 b_r
    // surge, sway, heave, roll: 0
    // pitch:  +U n3
    // yaw:    -U n2
    if (mode == 4)
        b = -U * element->Nvec.col(2);
    else if (mode == 5)
        b = U * element->Nvec.col(1);

    return b;
}

double LinearCumminsTDGF::weightedSurfaceIntegral(
    const Eigen::VectorXd& panelWeight,
    const Eigen::VectorXd& phi) const
{
    return (element->Area.array() * panelWeight.array()).matrix().dot(phi);
}


double LinearCumminsTDGF::computeModeIntegral(int mode, const Eigen::VectorXd& phi) const
{
    return -element->ArInt.col(mode).dot(phi);
}


Eigen::VectorXd LinearCumminsTDGF::applyWaterlineBracket1(
    int lagIdx,
    const Eigen::VectorXd& panelValues) const
{
    const int NE = SeakeepingCfg.Panel.NE;
    const int nWL = element->n_WL;

    Eigen::VectorXd out = Eigen::VectorXd::Zero(NE);

    const auto& Gx = Gw_dx->at(lagIdx);
    const auto& Gl = Gw_dl->at(lagIdx);
    const auto& Gwi = Gw->at(lagIdx);

    for (int l = 0; l < nWL; ++l)
    {
        const int l2 = (l + 1) % nWL;

        const int p = PotL_idx(l);
        const int p2 = PotL_idx(l2);
        const double chiL = panelValues(p);

        const double n1 = wlN1_(l);
        const double n2 = wlN2_(l);

        const double n1_2 = wlN1_(l2);
        const double n2_2 = wlN2_(l2);
        const double dn1n2dl = wlDn1n2Dl_(l);


        // 论文 3.15 中的导数是对源点坐标(ξ,l)求导；当前表中存的是对场点坐标(x,y)求导。
        // 因 G = G(x-ξ, y-η, ...)，故 ∂G/∂ξ = -∂G/∂x，源点切向导数也同理取反。

        //out.noalias() += chiL * wlSegLen_(l) * (
        //    -n1 * Gx.col(l).cast<double>());
        //out.noalias() += chiL * (n1 * n2 *
        //    (Gwi.col(l2).cast<double>()  - Gwi.col(l).cast<double>())
        //    + Gwi.col(l).cast<double>() * (n1_2 * n2_2 - n1 * n2));


        out.noalias() += chiL * (
            -n1 * Gx.col(l).cast<double>()
            - (n1 * n2) * Gl.col(l).cast<double>()
            + wlDn1n2Dl_(l) * Gwi.col(l).cast<double>());
    }

    return out;
}


Eigen::VectorXd LinearCumminsTDGF::applyWaterlineN1SqG(
    int lagIdx,
    const Eigen::VectorXd& panelValues) const
{
    const int NE = SeakeepingCfg.Panel.NE;
    const int nWL = element->n_WL;

    Eigen::VectorXd out = Eigen::VectorXd::Zero(NE);
    const auto& Gwi = Gw->at(lagIdx);

    for (int l = 0; l < nWL; ++l)
    {
        const int p = PotL_idx(l);

        // 对应 ∮ n1^2 * G dl
        //const double val =
        //    wlSegLen_(l) * wlN1_(l) * wlN1_(l) * panelValues(p);

        const double val =
            wlN1_(l) * wlN1_(l) * panelValues(p);

        out.noalias() += val * Gwi.col(l).cast<double>();
    }

    return out;
}


Eigen::VectorXd LinearCumminsTDGF::applyWaterlineG(
    int lagIdx,
    const Eigen::VectorXd& panelValues) const
{
    const int NE = SeakeepingCfg.Panel.NE;
    const int nWL = element->n_WL;

    Eigen::VectorXd out = Eigen::VectorXd::Zero(NE);
    const auto& Gwi = Gw->at(lagIdx);
    const auto& Gwt = Gw_dt->at(lagIdx);

    for (int l = 0; l < nWL; ++l)
    {
        const int p = PotL_idx(l);
        //const double val = wlSegLen_(l) * wlN1_(l) * panelValues(p);
        const double val = wlN1_(l) * panelValues(p);
        out.noalias() += val * Gwi.col(l).cast<double>();
    }

    return out;
}


namespace
{
    // 水线积分项的符号/幅值定位开关（radiation χ 递推专用）。
    // 论文 (3.5) 与 (3.11c)/(3.12c)/(3.17c) 中这两项的符号相互矛盾（把
    // φ=ψ1δ+ψ2H+χ 代入是纯替换，不应改变符号），文献内部无法裁决正确符号。
    // 通过环境变量做符号扫描（共 4 组合），配合 K→B(ω) 与频域阻尼对比来判定：
    //   SHIPSIM_WL1_SCALE : U²/g 括号项（applyWaterlineBracket1）
    //   SHIPSIM_WL2_SCALE : 2U/g 前向差分项（applyWaterlineG，含移到左端的 χ_N 部分）
    // 取值：1（默认，当前行为）/ -1（翻号）/ 0（关闭该项，用于消融定位）。
    double envScaleOrOne(const char* name)
    {
        if (const char* v = std::getenv(name))
        {
            try { return std::stod(v); } catch (...) {}
        }
        return 1.0;
    }

    double waterlineTerm1Scale()
    {
        static const double s = envScaleOrOne("SHIPSIM_WL1_SCALE");
        return s;
    }

    double waterlineTerm2Scale()
    {
        static const double s = envScaleOrOne("SHIPSIM_WL2_SCALE");
        return s;
    }
}

Eigen::VectorXd LinearCumminsTDGF::solveChiAtTime(
    int tN,
    const Eigen::VectorXd& bc_a,
    const Eigen::VectorXd& bc_b,
    const Eigen::VectorXd& psi1,
    const Eigen::VectorXd& psi2,
    const std::vector<Eigen::VectorXd>& chi_history) const
{
    //const auto& dGz = (Solver == "Source") ? Green_dnS : Green_dnP;
    const auto& dGz = Green_dnP;

    const int NE = SeakeepingCfg.Panel.NE;

    if (tN == 0)
        return Eigen::VectorXd::Zero(NE);

    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(NE);

    // Eq. (3.17c): instantaneous term G*a_r and -dG/dn * psi1_r.
    rhs.noalias() += (Green->at(tN - 1) * bc_a.cast<GScalar>()).cast<double>();
    rhs.noalias() -= (dGz->at(tN - 1) * psi1.cast<GScalar>()).cast<double>();

    // History convolution terms in Eq. (3.17c).
    for (int n = 0; n <= tN - 1; ++n)
    {
        const int lag = tN - n - 1; // lag=0 corresponds to tau = dt

        rhs.noalias() += dtShared_ * (Green->at(lag) * bc_b.cast<GScalar>()).cast<double>();
        rhs.noalias() -= dtShared_ * (dGz->at(lag) * psi2.cast<GScalar>()).cast<double>();
        rhs.noalias() -= dtShared_ * (dGz->at(lag) * chi_history[static_cast<std::size_t>(n)].cast<GScalar>()).cast<double>();

        rhs.noalias() += (waterlineTerm1Scale() * UsquareG * dtShared_) * applyWaterlineBracket1(
            lag, chi_history[static_cast<std::size_t>(n)]);
    }

    //Eq. (3.17c) last waterline-history term contains (chi_{n+1} - chi_n).
    //The n = N-1 contribution introduces the current unknown chi_N. We keep the
    //discretization explicit for n <= N-2 and move the chi_N part to the left-hand side.
    // 

    const double wl2 = waterlineTerm2Scale();
    if (std::abs(U) > 1e-14 && wl2 != 0.0)
    {
        for (int n = 0; n <= tN - 2; ++n)
        {
            const int lag = tN - n - 1;
            rhs.noalias() += (wl2 * 2.0 * U / G) * applyWaterlineG(
                lag,
                chi_history[static_cast<std::size_t>(n + 1)] -
                chi_history[static_cast<std::size_t>(n)]);

            //rhs.noalias() += (2.0 * U / G * dtShared_) * applyWaterlineG(
            //    lag,
            //    chi_history[static_cast<std::size_t>(n)]);

            //std::cout<<"tn:\t"<< n<< "\n";
        }

        rhs.noalias() -= (wl2 * 2.0 * U / G) * applyWaterlineG(
            0, chi_history[static_cast<std::size_t>(tN - 1)]);

        //rhs.noalias() += (2.0 * U / G * dtShared_) * applyWaterlineG(
        //    0, chi_history[static_cast<std::size_t>(tN - 1)]);

         //Solve (A + M) chi_N = rhs using the stored LU of A only.
         //Fixed-point iteration: A x^{m+1} = rhs - M x^{m}.

        Eigen::VectorXd x = solveSourceLinearSystem(rhs);
        for (int iter = 0; iter < 1000; ++iter)
        {
            const Eigen::VectorXd corr = (wl2 * 2.0 * U / G) * applyWaterlineG(0, x);
            const Eigen::VectorXd xNew = solveSourceLinearSystem(rhs + corr);

            const double denom = std::max(1.0, xNew.norm());
            if ((xNew - x).norm() / denom < 1e-10)
                return xNew;

            x = 0.6 * x + 0.4 * xNew;
        }
        return x;
    }
    return solveSourceLinearSystem(rhs);
}


void LinearCumminsTDGF::buildKernelDirectTimeDomain(double Fn, RadiationKernelData& kernel)
{
    // 路径分流：非均匀网格已配置则走 _NonUniform，否则继续走原等距实现。
    if (useNonUniformChiGrid())
    {
        buildKernelDirectTimeDomain_NonUniform(Fn, kernel);
        return;
    }

    const int D = SeakeepingCfg.DOF;
    const int NE = SeakeepingCfg.Panel.NE;

    kernel.A_inf = Eigen::MatrixXd::Zero(D, D);
    kernel.B = Eigen::MatrixXd::Zero(D, D);
    kernel.C_prime = Eigen::MatrixXd::Zero(D, D);
    kernel.Klag.assign(static_cast<std::size_t>(TGShared_ + 1), Eigen::MatrixXd::Zero(D, D));
    kernel.K0 = Eigen::MatrixXd::Zero(D, D);

    std::cout << "[LinearCumminsTDGF] Building radiation kernel according to paper Sec. 3.3 / 3.4\n";
    std::cout << "  Fn = " << Fn << ", dt = " << dtShared_ << ", TG = " << TGShared_ << "\n";
    std::cout << "  waterline-term scales: WL1(U^2/g)=" << waterlineTerm1Scale()
              << "  WL2(2U/g)=" << waterlineTerm2Scale()
              << "  (env SHIPSIM_WL1_SCALE / SHIPSIM_WL2_SCALE; 1=default, -1=flip, 0=off)\n";

    // Mode iterations are independent (each writes a distinct column k of the kernel
    // matrices, and all class-state reads — Green tables, element LU, geometry — are const
    // and thread-safe). Inside the parallel region we force Eigen to single-thread to keep
    // the per-thread compute pinned and avoid nested-thread oversubscription.
#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < D; ++k)
    {
        shipsim::EigenSingleThreadGuard _g;

        const int srcMode = SeakeepingCfg.modes[k];

        std::cout << "  [Mode " << srcMode << "] Processing...\n";

        const Eigen::VectorXd bc_a = buildBoundaryCondition_a(srcMode);
        const Eigen::VectorXd bc_b = buildBoundaryCondition_b(srcMode);

        const Eigen::VectorXd psi1 = solveRankineBIE(bc_a);
        const Eigen::VectorXd psi2 = solveRankineBIE(bc_b);

        for (int i = 0; i < D; ++i)
        {
            const int respMode = SeakeepingCfg.modes[i];
            const Eigen::VectorXd b_resp = buildBoundaryCondition_b(respMode);

            // 3.21a / 3.22a / 3.22b
            kernel.A_inf(i, k) =
                rho * computeModeIntegral(respMode, psi1);

            kernel.B(i, k) =
                rho * (computeModeIntegral(respMode, psi2)
                    - weightedSurfaceIntegral(b_resp, psi1));

            kernel.C_prime(i, k) =
                -rho * weightedSurfaceIntegral(b_resp, psi2);
        }


        // -------------------------------------------------------------------------
        // Step 1: 先完整求 raw chi history
        // -------------------------------------------------------------------------
        std::vector<Eigen::VectorXd> chi_raw_history;
        chi_raw_history.reserve(static_cast<std::size_t>(TGShared_ + 1));
        chi_raw_history.push_back(Eigen::VectorXd::Zero(NE));

        for (int tN = 1; tN <= TGShared_; ++tN)
        {
            Eigen::VectorXd chiN = solveChiAtTime(
                tN, bc_a, bc_b, psi1, psi2, chi_raw_history);

            chi_raw_history.push_back(chiN);

            if ((tN % 50) == 0 || tN == TGShared_)
            {
                std::cout << "    raw chi solved up to step "
                    << tN << "/" << TGShared_
                    << " (t=" << tN * dtShared_ << "s)\n";
            }
        }

        // -------------------------------------------------------------------------
        // Step 2: 先做空间积分得到光滑的一维标量序列
        //         m_i(t) = ρ ∬ a_{resp_i} χ dS，
        //         再对该标量序列做时间微分。
        //   依据 3.22c: a_r 不含时，故
        //   ρ∬ a_r ∂χ/∂t dS = d/dt [ ρ∬ a_r χ dS ]，
        //   把 d/dt 移到空间积分外面：空间积分先把面元噪声平滑掉，
        //   再对干净的一维序列求导，避免对面元 χ 场逐点微分放大振荡。
        // -------------------------------------------------------------------------
        std::vector<Eigen::VectorXd> mHist(
            static_cast<std::size_t>(TGShared_ + 1),
            Eigen::VectorXd::Zero(D));

        for (int tN = 0; tN <= TGShared_; ++tN)
        {
            const Eigen::VectorXd& chi_t =
                chi_raw_history[static_cast<std::size_t>(tN)];
            for (int i = 0; i < D; ++i)
            {
                const int respMode = SeakeepingCfg.modes[i];
                mHist[static_cast<std::size_t>(tN)](i) =
                    rho * computeModeIntegral(respMode, chi_t);
            }
        }

        // d/dt [ ρ∬ a_resp χ dS ]：对一维标量序列复用同一差分模板
        const std::vector<Eigen::VectorXd> dmHist =
            differentiateHistory(mHist);

        // -------------------------------------------------------------------------
        // Step 3: 调试输出
        //   chi            : raw chi 面元历史
        //   m_integrated   : ρ∬ a_resp χ dS 空间积分后的标量序列 (列=响应模态序)
        //   dm_dt          : 对 m_integrated 时间求导 (即 k_dchi)
        // -------------------------------------------------------------------------
        writeChiDebugCsv(Fn, srcMode, chi_raw_history);

        {
            const double L = ShipCfg.Geometry.Length;
            const double tScale = std::sqrt(G / L);

            auto writeModeSeries =
                [&](const std::string& name,
                    const std::vector<Eigen::VectorXd>& hist)
                {
                    std::filesystem::create_directories(
                        std::filesystem::path(filePath) / "kernel_cache");

                    std::ostringstream ss;
                    ss << filePath << "kernel_cache/" << name << "_debug"
                        << "_Fn" << keyDouble(Fn)
                        << "_mode" << srcMode
                        << "_dt" << keyDouble(dtShared_)
                        << "_TG" << TGShared_ << ".csv";

                    std::ofstream out(ss.str());
                    if (!out.is_open()) return;

                    out << std::setprecision(17);
                    out << "step,t,t_nd";
                    for (int i = 0; i < D; ++i)
                        out << ",resp_mode" << SeakeepingCfg.modes[i];
                    out << "\n";

                    for (int n = 0; n < static_cast<int>(hist.size()); ++n)
                    {
                        const double t = n * dtShared_;
                        out << n << "," << t << "," << (t * tScale);
                        for (int i = 0; i < D; ++i)
                            out << "," << hist[static_cast<std::size_t>(n)](i);
                        out << "\n";
                    }
                    std::cout << "[LinearCumminsTDGF] " << name
                        << " debug csv saved: " << ss.str() << "\n";
                };

            writeModeSeries("m_integrated", mHist);
            writeModeSeries("dm_dt", dmHist);
        }

        // -------------------------------------------------------------------------
        // Step 4: 组装 Klag (3.22c: K = ρ∬ a_r ∂χ/∂t dS − ρ∬ b_r χ dS)
        //         第一项 k_dchi = d/dt[ρ∬ a_r χ dS] = dmHist
        // -------------------------------------------------------------------------
        for (int tN = 0; tN <= TGShared_; ++tN)
        {
            const Eigen::VectorXd& chi_used = chi_raw_history[static_cast<std::size_t>(tN)];

            for (int i = 0; i < D; ++i)
            {
                const int respMode = SeakeepingCfg.modes[i];
                const Eigen::VectorXd b_resp = buildBoundaryCondition_b(respMode);

                const double k_dchi = dmHist[static_cast<std::size_t>(tN)](i);
                const double k_bchi = -rho * weightedSurfaceIntegral(b_resp, chi_used);

                kernel.Klag[static_cast<std::size_t>(tN)](i, k) = k_dchi + k_bchi;
            }
        }

        std::cout << "  [Mode " << srcMode << "] Completed\n";
    }

    // Section 3.4 defines K(t) for the memory convolution. We do not invent K(0+)
    // by extrapolation; the online convolution uses only strictly positive lags.
    //kernel.Klag[0].setZero();
    //kernel.Klag[0] = kernel.Klag[1];
    //kernel.K0 = Eigen::MatrixXd::Zero(D, D);
    kernel.K0 = kernel.Klag[0];

    // 同步填写 Klag_times（当前 build 路径是等距 dtShared_，t[i]=i*dt；
    // 后续切到非均匀网格时，这里要换成实际节点时间）。
    kernel.Klag_times.resize(kernel.Klag.size());
    for (std::size_t i = 0; i < kernel.Klag.size(); ++i)
        kernel.Klag_times[i] = static_cast<double>(i) * dtShared_;

    std::cout << "[LinearCumminsTDGF] Radiation kernel built successfully\n"
        << "  A_inf range: [" << kernel.A_inf.minCoeff() << ", " << kernel.A_inf.maxCoeff() << "]\n"
        << "  B range: [" << kernel.B.minCoeff() << ", " << kernel.B.maxCoeff() << "]\n"
        << "  C' range: [" << kernel.C_prime.minCoeff() << ", " << kernel.C_prime.maxCoeff() << "]\n";

    writeKernelDiagnostics(Fn, kernel, "raw", dtShared_);
}

// =====================================================================
// 非均匀 chi 网格构建路径
// =====================================================================
//
// 假设：
//   - dt_fine == dtShared_（buildChiTimeGridFromConfig 已校验）
//   - 所有 chi 节点 idx*dt_fine 都落在 Green 表（按 dtShared_ 等距索引）的有效区间内
//   - chiIdxOnFine_ / chiTimes_ 严格递增，chiIdxOnFine_[0]==0
//
// chi 卷积近似（与 uniform 路径同阶 — left-rectangle）：
//   ∫₀^{t_curr} G(t_curr - τ) · f(τ) dτ
//     ≈ Σ_{j=0}^{iCurr-1} (chiTimes_[j+1] - chiTimes_[j]) · G[lag_j] · f(chiTimes_[j])
//   其中 lag_j = (idxCurr - chiIdxOnFine_[j]) - 1（Green 表索引，与 uniform 一致）。
//
// memory_cutoff_lag：若 (idxCurr - idxJ) > cutoff，该项直接跳过（认为 G 已衰减）。
// =====================================================================
void LinearCumminsTDGF::buildKernelDirectTimeDomain_NonUniform(double Fn, RadiationKernelData& kernel)
{
    const int D = SeakeepingCfg.DOF;
    const int NE = SeakeepingCfg.Panel.NE;
    const int N_chi = static_cast<int>(chiTimes_.size());

    kernel.A_inf   = Eigen::MatrixXd::Zero(D, D);
    kernel.B       = Eigen::MatrixXd::Zero(D, D);
    kernel.C_prime = Eigen::MatrixXd::Zero(D, D);
    kernel.Klag.assign(static_cast<std::size_t>(N_chi), Eigen::MatrixXd::Zero(D, D));
    kernel.K0      = Eigen::MatrixXd::Zero(D, D);
    kernel.Klag_times = chiTimes_;   // 输出时间向量直接接管 chi 网格

    std::cout << "[LinearCumminsTDGF] Building radiation kernel (NON-UNIFORM chi grid)\n"
              << "  Fn = " << Fn
              << ", dt_fine = " << chiDtFine_
              << ", chi nodes = " << N_chi
              << ", tMax = " << chiTimes_.back() << "s"
              << ", cutoff_lag = " << chiMemoryCutoffLag_ << "\n";

#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < D; ++k)
    {
        shipsim::EigenSingleThreadGuard _g;

        const int srcMode = SeakeepingCfg.modes[k];

        std::cout << "  [NU Mode " << srcMode << "] Processing...\n";

        const Eigen::VectorXd bc_a = buildBoundaryCondition_a(srcMode);
        const Eigen::VectorXd bc_b = buildBoundaryCondition_b(srcMode);

        const Eigen::VectorXd psi1 = solveRankineBIE(bc_a);
        const Eigen::VectorXd psi2 = solveRankineBIE(bc_b);

        // A_inf / B / C' —— 静态系数，与 uniform 路径完全一致。
        for (int i = 0; i < D; ++i)
        {
            const int respMode = SeakeepingCfg.modes[i];
            const Eigen::VectorXd b_resp = buildBoundaryCondition_b(respMode);

            kernel.A_inf(i, k) =
                rho * computeModeIntegral(respMode, psi1);

            kernel.B(i, k) =
                rho * (computeModeIntegral(respMode, psi2)
                    - weightedSurfaceIntegral(b_resp, psi1));

            kernel.C_prime(i, k) =
                -rho * weightedSurfaceIntegral(b_resp, psi2);
        }

        // ---- chi 历史：仅在 chi 节点上求解 ----
        std::vector<Eigen::VectorXd> chi_raw_history;
        chi_raw_history.reserve(static_cast<std::size_t>(N_chi));
        chi_raw_history.push_back(Eigen::VectorXd::Zero(NE));   // chi(0)=0

        for (int iCurr = 1; iCurr < N_chi; ++iCurr)
        {
            Eigen::VectorXd chiN = solveChiAtTime_NonUniform(
                iCurr, bc_a, bc_b, psi1, psi2, chi_raw_history);
            chi_raw_history.push_back(chiN);

            if ((iCurr % 25) == 0 || iCurr == N_chi - 1)
            {
                std::cout << "    [NU mode " << srcMode << "] chi solved up to node "
                          << iCurr << "/" << (N_chi - 1)
                          << " (t=" << chiTimes_[iCurr] << "s)\n";
            }
        }

        // ---- 空间积分得到 m(t) 标量序列 ----
        std::vector<Eigen::VectorXd> mHist(
            static_cast<std::size_t>(N_chi),
            Eigen::VectorXd::Zero(D));
        for (int t = 0; t < N_chi; ++t)
        {
            const Eigen::VectorXd& chi_t = chi_raw_history[static_cast<std::size_t>(t)];
            for (int i = 0; i < D; ++i)
            {
                const int respMode = SeakeepingCfg.modes[i];
                mHist[static_cast<std::size_t>(t)](i) =
                    rho * computeModeIntegral(respMode, chi_t);
            }
        }

        // ---- 对 m(t) 做非均匀 3 点中心差分得到 dm/dt ----
        const std::vector<Eigen::VectorXd> dmHist =
            differentiateHistoryNonUniform(mHist, chiTimes_);

        // ---- 装配 Klag：K = dm/dt − ρ ∫ b_r · χ dS ----
        for (int t = 0; t < N_chi; ++t)
        {
            const Eigen::VectorXd& chi_used = chi_raw_history[static_cast<std::size_t>(t)];

            for (int i = 0; i < D; ++i)
            {
                const int respMode = SeakeepingCfg.modes[i];
                const Eigen::VectorXd b_resp = buildBoundaryCondition_b(respMode);

                const double k_dchi = dmHist[static_cast<std::size_t>(t)](i);
                const double k_bchi = -rho * weightedSurfaceIntegral(b_resp, chi_used);

                kernel.Klag[static_cast<std::size_t>(t)](i, k) = k_dchi + k_bchi;
            }
        }

        std::cout << "  [NU Mode " << srcMode << "] Completed\n";
    }

    kernel.K0 = kernel.Klag[0];

    std::cout << "[LinearCumminsTDGF] Non-uniform radiation kernel built successfully\n"
              << "  A_inf range: [" << kernel.A_inf.minCoeff() << ", " << kernel.A_inf.maxCoeff() << "]\n"
              << "  B range: [" << kernel.B.minCoeff() << ", " << kernel.B.maxCoeff() << "]\n"
              << "  C' range: [" << kernel.C_prime.minCoeff() << ", " << kernel.C_prime.maxCoeff() << "]\n";

    writeKernelDiagnostics(Fn, kernel, "raw_nu", chiDtFine_);
}

// =====================================================================
// 非均匀 chi 求解器
// =====================================================================
//
// 与 solveChiAtTime（uniform 版）语义一致：
//   - 瞬时项：lag = idxCurr - 1（Green 表索引）
//   - 历史项：lag = (idxCurr - chiIdxOnFine_[j]) - 1，权重 = chiTimes_[j+1] - chiTimes_[j]
//   - U≠0 时的边界差分项 + fixed-point 迭代：lag 也是 (idxCurr - chiIdxOnFine_[iCurr-1]) - 1
//
// memory_cutoff_lag 用 fine-grid 单位的物理 lag 判断（idxCurr - chiIdxOnFine_[j]）。
// =====================================================================
Eigen::VectorXd LinearCumminsTDGF::solveChiAtTime_NonUniform(
    int iCurr,
    const Eigen::VectorXd& bc_a,
    const Eigen::VectorXd& bc_b,
    const Eigen::VectorXd& psi1,
    const Eigen::VectorXd& psi2,
    const std::vector<Eigen::VectorXd>& chi_history) const
{
    const auto& dGz = Green_dnP;

    const int NE = SeakeepingCfg.Panel.NE;
    if (iCurr <= 0)
        return Eigen::VectorXd::Zero(NE);

    const int idxCurr = chiIdxOnFine_[static_cast<std::size_t>(iCurr)];
    const int cutoff = (chiMemoryCutoffLag_ > 0)
                       ? chiMemoryCutoffLag_
                       : std::numeric_limits<int>::max();

    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(NE);

    // ---- 瞬时项 (Eq. 3.17c)：lag = idxCurr - 1 ----
    // 保留与 uniform 路径相同的索引语义；瞬时项不受 cutoff 影响。
    {
        const int lag0 = idxCurr - 1;
        rhs.noalias() += (Green->at(static_cast<std::size_t>(lag0))
                            * bc_a.cast<GScalar>()).cast<double>();
        rhs.noalias() -= (dGz->at(static_cast<std::size_t>(lag0))
                            * psi1.cast<GScalar>()).cast<double>();
    }

    // ---- 历史卷积（非均匀 left-rectangle）----
    for (int j = 0; j < iCurr; ++j)
    {
        const int idxJ = chiIdxOnFine_[static_cast<std::size_t>(j)];
        const int physLag = idxCurr - idxJ;          // ≥ 1 (fine units)
        if (physLag > cutoff) continue;              // memory truncation
        const int lag = physLag - 1;                 // Green 表索引

        const double dj = chiTimes_[static_cast<std::size_t>(j + 1)]
                          - chiTimes_[static_cast<std::size_t>(j)];
        const auto& chi_j = chi_history[static_cast<std::size_t>(j)];

        rhs.noalias() += dj * (Green->at(static_cast<std::size_t>(lag))
                                 * bc_b.cast<GScalar>()).cast<double>();
        rhs.noalias() -= dj * (dGz->at(static_cast<std::size_t>(lag))
                                 * psi2.cast<GScalar>()).cast<double>();
        rhs.noalias() -= dj * (dGz->at(static_cast<std::size_t>(lag))
                                 * chi_j.cast<GScalar>()).cast<double>();

        rhs.noalias() += (waterlineTerm1Scale() * UsquareG * dj) * applyWaterlineBracket1(lag, chi_j);
    }

    // ---- U ≠ 0 的边界差分项 ----
    const double wl2 = waterlineTerm2Scale();
    if (std::abs(U) > 1e-14 && wl2 != 0.0)
    {
        // (chi[j+1] - chi[j]) loop，覆盖到 iCurr-2（即不含最后一个区间）
        for (int j = 0; j < iCurr - 1; ++j)
        {
            const int idxJ = chiIdxOnFine_[static_cast<std::size_t>(j)];
            const int physLag = idxCurr - idxJ;
            if (physLag > cutoff) continue;
            const int lag = physLag - 1;

            const auto& chi_j  = chi_history[static_cast<std::size_t>(j)];
            const auto& chi_j1 = chi_history[static_cast<std::size_t>(j + 1)];

            rhs.noalias() += (wl2 * 2.0 * U / G) * applyWaterlineG(
                lag, chi_j1 - chi_j);
        }

        // 最后一个区间 [chiTimes_[iCurr-1], chiTimes_[iCurr]] 由 chi_curr 与 chi[iCurr-1] 分担：
        // 减去 chi[iCurr-1] 项；fixed-point 迭代里加上 chi_curr 项 → 等价 (chi_curr - chi[iCurr-1])
        const int lagBoundary =
            idxCurr - chiIdxOnFine_[static_cast<std::size_t>(iCurr - 1)] - 1;
        const auto& chi_prev = chi_history[static_cast<std::size_t>(iCurr - 1)];

        rhs.noalias() -= (wl2 * 2.0 * U / G) * applyWaterlineG(lagBoundary, chi_prev);

        // (A + M) chi_curr = rhs，  M = -(2U/g) · waterlineG(lagBoundary, ·)
        Eigen::VectorXd x = solveSourceLinearSystem(rhs);
        for (int iter = 0; iter < 1000; ++iter)
        {
            const Eigen::VectorXd corr = (wl2 * 2.0 * U / G) * applyWaterlineG(lagBoundary, x);
            const Eigen::VectorXd xNew = solveSourceLinearSystem(rhs + corr);

            const double denom = std::max(1.0, xNew.norm());
            if ((xNew - x).norm() / denom < 1e-10)
                return xNew;

            x = 0.6 * x + 0.4 * xNew;
        }
        return x;
    }

    return solveSourceLinearSystem(rhs);
}

// =====================================================================
// 非均匀差分：5 点 4 阶中心（内点）+ 3 点 2 阶（次端点）+ 单边（端点）
// 用 Lagrange 插值多项式求导得到非均匀权重；
// 等距情形下退化为标准 4 阶中心差分 (-1, 8, 0, -8, 1)/(12h)。
// =====================================================================
namespace {
    // 计算 5 点 Lagrange 插值在 xc 处的一阶导权重 w[0..4]。
    // 通用 (xc 可以是任一节点)：f'(xc) ≈ Σ w[k] · f(x[k])。
    inline void nu5PointDiffWeights(
        const double x[5], double xc, double w[5])
    {
        for (int k = 0; k < 5; ++k)
        {
            // denom = Π_{j≠k} (x[k] - x[j])
            double denom = 1.0;
            for (int j = 0; j < 5; ++j)
                if (j != k) denom *= (x[k] - x[j]);

            // num = d/dx [Π_{j≠k}(x - x[j])] at xc
            //     = Σ_{m≠k} Π_{j≠k,m} (xc - x[j])
            double num = 0.0;
            for (int m = 0; m < 5; ++m)
            {
                if (m == k) continue;
                double prod = 1.0;
                for (int j = 0; j < 5; ++j)
                {
                    if (j == k || j == m) continue;
                    prod *= (xc - x[j]);
                }
                num += prod;
            }
            w[k] = num / denom;
        }
    }
}

std::vector<Eigen::VectorXd> LinearCumminsTDGF::differentiateHistoryNonUniform(
    const std::vector<Eigen::VectorXd>& history,
    const std::vector<double>& times) const
{
    if (history.empty())
        return history;

    const int N = static_cast<int>(history.size());
    const int D = static_cast<int>(history.front().size());

    std::vector<Eigen::VectorXd> d(
        static_cast<std::size_t>(N),
        Eigen::VectorXd::Zero(D));

    if (N < 2) return d;
    if (static_cast<int>(times.size()) != N)
        throw std::runtime_error(
            "differentiateHistoryNonUniform: times/history size mismatch");

    // N == 2：只能一阶差分
    if (N == 2)
    {
        const double h = times[1] - times[0];
        if (h > 0.0)
        {
            d[0] = (history[1] - history[0]) / h;
            d[1] = d[0];
        }
        return d;
    }

    // 3 点非均匀中心差分权重
    auto diff3 = [&](int i)
    {
        const double h1 = times[static_cast<std::size_t>(i)]
                          - times[static_cast<std::size_t>(i - 1)];
        const double h2 = times[static_cast<std::size_t>(i + 1)]
                          - times[static_cast<std::size_t>(i)];
        if (h1 <= 0.0 || h2 <= 0.0)
        {
            d[static_cast<std::size_t>(i)] = Eigen::VectorXd::Zero(D);
            return;
        }
        const double w_m = -h2 / (h1 * (h1 + h2));
        const double w_0 = (h2 - h1) / (h1 * h2);
        const double w_p =  h1 / (h2 * (h1 + h2));
        d[static_cast<std::size_t>(i)] =
            w_m * history[static_cast<std::size_t>(i - 1)]
          + w_0 * history[static_cast<std::size_t>(i)]
          + w_p * history[static_cast<std::size_t>(i + 1)];
    };

    // 5 点非均匀中心差分（4 阶精度），用于内部点 2..N-3
    auto diff5 = [&](int i)
    {
        const double xs[5] = {
            times[static_cast<std::size_t>(i - 2)],
            times[static_cast<std::size_t>(i - 1)],
            times[static_cast<std::size_t>(i)],
            times[static_cast<std::size_t>(i + 1)],
            times[static_cast<std::size_t>(i + 2)]
        };
        double w[5];
        nu5PointDiffWeights(xs, xs[2], w);
        d[static_cast<std::size_t>(i)] =
            w[0] * history[static_cast<std::size_t>(i - 2)]
          + w[1] * history[static_cast<std::size_t>(i - 1)]
          + w[2] * history[static_cast<std::size_t>(i)]
          + w[3] * history[static_cast<std::size_t>(i + 1)]
          + w[4] * history[static_cast<std::size_t>(i + 2)];
    };

    // 首点：一阶前向（O(h) — 仅 1 个端点，对全场影响极小）
    {
        const double h = times[1] - times[0];
        if (h > 0.0) d[0] = (history[1] - history[0]) / h;
    }

    // 次首点 i=1：3 点中心
    diff3(1);

    // 内点：5 点中心（4 阶精度）
    for (int i = 2; i <= N - 3; ++i)
        diff5(i);

    // 次末点 i=N-2：3 点中心
    if (N - 2 >= 1) diff3(N - 2);

    // 末点：一阶后向
    {
        const double h = times[static_cast<std::size_t>(N - 1)]
                       - times[static_cast<std::size_t>(N - 2)];
        if (h > 0.0)
            d[static_cast<std::size_t>(N - 1)] =
                (history[static_cast<std::size_t>(N - 1)]
               - history[static_cast<std::size_t>(N - 2)]) / h;
    }

    return d;
}

void LinearCumminsTDGF::writeKernelDiagnostics(
    double Fn,
    const RadiationKernelData& kernel,
    const std::string& stage,
    double dtForT) const
{
    const std::filesystem::path diagDir =
        std::filesystem::path(filePath) / "kernel_cache" / "diagnostics";
    std::filesystem::create_directories(diagDir);

    const std::string tag = stage + "_Fn" + keyDouble(Fn);
    const int D = static_cast<int>(kernel.A_inf.rows());

    const double dt = (dtForT > 0.0) ? dtForT
        : (kernel.dt > 0.0 ? kernel.dt : 1.0);

    auto dumpMatrix = [&](const std::string& fname, const Eigen::MatrixXd& M)
        {
            std::ofstream out((diagDir / fname).string());
            out << std::setprecision(12);
            for (int i = 0; i < M.rows(); ++i)
            {
                for (int j = 0; j < M.cols(); ++j)
                {
                    out << M(i, j);
                    if (j + 1 < M.cols()) out << ",";
                }
                out << "\n";
            }
        };

    dumpMatrix("Ainf_" + tag + ".csv", kernel.A_inf);
    dumpMatrix("B_" + tag + ".csv", kernel.B);
    dumpMatrix("Cprime_" + tag + ".csv", kernel.C_prime);
    dumpMatrix("K0_" + tag + ".csv", kernel.K0);

    // Symmetry deviation: ||M - M^T|| / ||M||  (smaller = more symmetric).
    auto symDev = [](const Eigen::MatrixXd& M) -> double
        {
            const double n = M.norm();
            if (n == 0.0) return 0.0;
            return (M - M.transpose()).norm() / n;
        };

    std::ofstream sym((diagDir / ("symmetry_" + tag + ".txt")).string());
    sym << "A_inf rel-asym: " << symDev(kernel.A_inf) << "\n"
        << "B     rel-asym: " << symDev(kernel.B) << "\n"
        << "C'    rel-asym: " << symDev(kernel.C_prime) << "\n";

    // K(t) main diagonals + dominant off-diagonals over time.
    std::ofstream kt((diagDir / ("Klag_diag_" + tag + ".csv")).string());
    kt << "t";
    for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j)
            kt << ",K_" << SeakeepingCfg.modes[i] << SeakeepingCfg.modes[j];
    kt << "\n";

    // 真实时间从 kernel.Klag_times 取（非均匀路径必须这样）；老的等距路径下
    // Klag_times[lag] = lag * dt，与原行为一致。
    const bool haveTimes =
        kernel.Klag_times.size() == kernel.Klag.size();

    kt << std::setprecision(10);
    for (std::size_t lag = 0; lag < kernel.Klag.size(); ++lag)
    {
        const double t = haveTimes
            ? kernel.Klag_times[lag]
            : (static_cast<double>(lag) * dt);
        kt << t;
        for (int i = 0; i < D; ++i)
            for (int j = 0; j < D; ++j)
                kt << "," << kernel.Klag[lag](i, j);
        kt << "\n";
    }

    std::cout << "[KernelDiag-" << stage << "] wrote " << diagDir.string()
        << " (Fn=" << Fn << ", dt=" << dt << ")\n";
}


void LinearCumminsTDGF::applyRadiationKernelOverride(
    RadiationKernelData& kernel) const
{
    for (auto& K : kernel.Klag)
        K.setZero();
    kernel.K0.setZero();

    const auto& ov = SeakeepingCfg.KradOverride;
    if (!ov.enabled || ov.file.empty())
        return;

    const std::string path = filePath + ov.file;
    std::ifstream in(path);
    if (!in.is_open())
    {
        std::cerr << "[Krad override] cannot open: " << path
            << "  (override skipped)\n";
        return;
    }

    auto stripBom = [](std::string& s)
    {
        if (s.size() >= 3
            && static_cast<unsigned char>(s[0]) == 0xEF
            && static_cast<unsigned char>(s[1]) == 0xBB
            && static_cast<unsigned char>(s[2]) == 0xBF)
        {
            s.erase(0, 3);
        }
    };

    auto split = [](const std::string& s) -> std::vector<std::string>
    {
        std::vector<std::string> out;
        std::stringstream ss(s);
        std::string c;
        while (std::getline(ss, c, ','))
        {
            while (!c.empty() && (c.front() == ' ' || c.front() == '\t' || c.front() == '\r'))
                c.erase(c.begin());
            while (!c.empty() && (c.back() == ' ' || c.back() == '\t' || c.back() == '\r'))
                c.pop_back();
            out.push_back(c);
        }
        return out;
    };

    // ---- map ship-DOF labels (2=heave,4=pitch) to kernel matrix indices ----
    auto idxOfCodeMode = [&](int codeMode) -> int
    {
        for (int j = 0; j < static_cast<int>(SeakeepingCfg.modes.size()); ++j)
            if (SeakeepingCfg.modes[static_cast<std::size_t>(j)] == codeMode) return j;
        return -1;
    };
    const int iH = idxOfCodeMode(2);   // heave (label 3)
    const int iP = idxOfCodeMode(4);   // pitch (label 5)
    if (iH < 0 || iP < 0)
    {
        std::cerr << "[Krad override] heave(2)/pitch(4) not both in modes."
            " (override skipped)\n";
        return;
    }
    const int rc[2] = { iH, iP };

    // ---- parse optional "# A_inf / # B / # C_prime" 2x2 blocks ----
    // File layout:
    //   # A_inf (dimensional)
    //   a00, a01
    //   a10, a11
    //   # B (dimensional)
    //   ...
    //   # C_prime (dimensional)
    //   ...
    //   t,K33,K35,K53,K55
    //   <rows>
    // Each block (if present) overrides the heave/pitch 2x2 sub-block of the
    // matching kernel matrix, with the rest of the matrix zeroed (mirrors the
    // K-override convention: heave/pitch only).
    auto readMatrix2x2 = [&](Eigen::Matrix2d& M) -> bool
    {
        std::string l;
        int got = 0;
        while (got < 2 && std::getline(in, l))
        {
            stripBom(l);
            if (l.empty()) continue;
            const auto f = split(l);
            if (f.size() < 2) continue;
            try {
                M(got, 0) = std::stod(f[0]);
                M(got, 1) = std::stod(f[1]);
            }
            catch (...) { return false; }
            ++got;
        }
        return got == 2;
    };

    auto applyTo = [&](Eigen::MatrixXd& dst, const Eigen::Matrix2d& M)
    {
        if (dst.rows() == 0 || dst.cols() == 0)
            dst = Eigen::MatrixXd::Zero(SeakeepingCfg.DOF, SeakeepingCfg.DOF);
        else
            dst.setZero();
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b)
                dst(rc[a], rc[b]) = M(a, b);
    };

    std::string line;
    bool sawA = false, sawB = false, sawC = false;
    bool kHeaderFound = false;
    std::vector<std::string> hdr;

    while (std::getline(in, line))
    {
        stripBom(line);
        std::string trimmed = line;
        while (!trimmed.empty() && (trimmed.front() == ' ' || trimmed.front() == '\t' || trimmed.front() == '\r'))
            trimmed.erase(trimmed.begin());
        if (trimmed.empty()) continue;

        if (trimmed.front() == '#')
        {
            // identify which matrix block this is
            std::string tag = trimmed.substr(1);
            while (!tag.empty() && (tag.front() == ' ' || tag.front() == '\t'))
                tag.erase(tag.begin());
            Eigen::Matrix2d M;
            if (tag.rfind("A_inf", 0) == 0 || tag.rfind("A (", 0) == 0)
            {
                if (!readMatrix2x2(M))
                {
                    std::cerr << "[Krad override] failed to read A_inf 2x2 block.\n";
                    return;
                }
                applyTo(kernel.A_inf, M);
                sawA = true;
            }
            else if (tag.rfind("B", 0) == 0)
            {
                if (!readMatrix2x2(M))
                {
                    std::cerr << "[Krad override] failed to read B 2x2 block.\n";
                    return;
                }
                applyTo(kernel.B, M);
                sawB = true;
            }
            else if (tag.rfind("C_prime", 0) == 0 || tag.rfind("C'", 0) == 0 || tag.rfind("C ", 0) == 0)
            {
                if (!readMatrix2x2(M))
                {
                    std::cerr << "[Krad override] failed to read C_prime 2x2 block.\n";
                    return;
                }
                applyTo(kernel.C_prime, M);
                sawC = true;
            }
            // otherwise: ignore unknown comment block
            continue;
        }

        // first non-comment, non-empty line is the K column header
        hdr = split(trimmed);
        kHeaderFound = true;
        break;
    }

    if (!kHeaderFound)
    {
        std::cerr << "[Krad override] missing t,K33,K35,K53,K55 header. (override skipped)\n";
        return;
    }
    auto colOf = [&](const std::string& name) -> int
    {
        for (int j = 0; j < static_cast<int>(hdr.size()); ++j)
            if (hdr[static_cast<std::size_t>(j)] == name) return j;
        return -1;
    };

    const int cT = colOf("t");
    const int cK[2][2] = {
        { colOf("K33"), colOf("K35") },
        { colOf("K53"), colOf("K55") }
    };
    if (cT < 0 || cK[0][0] < 0 || cK[0][1] < 0 || cK[1][0] < 0 || cK[1][1] < 0)
    {
        std::cerr << "[Krad override] header must contain t,K33,K35,K53,K55."
            " (override skipped)\n";
        return;
    }

    // ---- read rows ----
    std::vector<double> tExt;
    std::vector<std::array<double, 4>> kExt;  // [K33,K35,K53,K55]
    while (std::getline(in, line))
    {
        if (line.empty()) continue;
        const std::vector<std::string> f = split(line);
        const int need = std::max({ cT, cK[0][0], cK[0][1], cK[1][0], cK[1][1] });
        if (static_cast<int>(f.size()) <= need) continue;

        tExt.push_back(std::stod(f[static_cast<std::size_t>(cT)]));
        kExt.push_back({
            std::stod(f[static_cast<std::size_t>(cK[0][0])]),
            std::stod(f[static_cast<std::size_t>(cK[0][1])]),
            std::stod(f[static_cast<std::size_t>(cK[1][0])]),
            std::stod(f[static_cast<std::size_t>(cK[1][1])]) });
    }
    if (tExt.size() < 2)
    {
        std::cerr << "[Krad override] not enough data rows. (override skipped)\n";
        return;
    }

    // ---- 完全替换 kernel 的时间网格：dt 和长度都由 CSV 决定 ----
    // 这样下游 resampleKernelToDt 会按 CSV 自身的总时长铺满 Konline，
    // 运动方程卷积会把 CSV 里存的 K 全部用完，不再受 TGShared_ 截断。
    const int D = SeakeepingCfg.DOF;
    const int N = static_cast<int>(tExt.size());
    const double tMin = tExt.front();
    const double tMax = tExt.back();
    const double dtCsv = (N > 1)
        ? (tMax - tMin) / static_cast<double>(N - 1)
        : ((kernel.dt > 0.0) ? kernel.dt : dtShared_);

    // 重建 Klag：CSV 多长，kernel 就多长
    kernel.dt = dtCsv;
    kernel.TG = N - 1;
    kernel.DOF = D;
    kernel.Klag.assign(static_cast<std::size_t>(N),
        Eigen::MatrixXd::Zero(D, D));
    // Klag_times 直接接管 CSV 的时间列（允许 CSV 非均匀；下游 resample 自适应）
    kernel.Klag_times = tExt;

    // 直接逐行拷贝 heave/pitch 2x2 子块，不再插值
    for (int n = 0; n < N; ++n)
    {
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b)
                kernel.Klag[static_cast<std::size_t>(n)](rc[a], rc[b]) =
                kExt[static_cast<std::size_t>(n)][a * 2 + b];
    }
    kernel.K0 = kernel.Klag[0];

    std::cout << "[Krad override] applied from " << path
        << "  (ext rows=" << tExt.size()
        << ", ext t=[" << tMin << "," << tMax << "]"
        << ", kernel grid N=" << N << ", dt=" << dtCsv
        << " (taken from CSV)"
        << ", heave idx=" << iH << ", pitch idx=" << iP
        << ", A_inf=" << (sawA ? "yes" : "no")
        << ", B=" << (sawB ? "yes" : "no")
        << ", C'=" << (sawC ? "yes" : "no") << ")\n";
}

void LinearCumminsTDGF::solveKernelCase(
    const CaseContext& ctx,
    const RadiationKernelData& kernel,
    RAO4DTable& tab)
{
    // Shared time-stepping core — see CumminsTimeStepper.h.
    // The coupled-window solver uses the same stepper, so the seakeeping
    // logic stays in lock-step between the two code paths.

    Eigen::MatrixXd Mphys, Bphys, C;
    SeakeepingDOF::buildSystemMatrices(ShipCfg, SeakeepingCfg, hs, Mphys, Bphys, C);

    const int D = SeakeepingCfg.DOF;
    const int N = static_cast<int>(ctx.ExcitingForce.rows());
    const double dt = ctx.dt;

    cummins::StepperConfig cfg;
    cfg.Mtotal = Mphys + kernel.A_inf;
    cfg.Bconst = Bphys + kernel.B + 0.5 * dt * kernel.K0;
    cfg.Ctotal = C + kernel.C_prime;
    cfg.Konline = resampleKernelToDt(kernel, dt);
    cfg.dt = dt;
    cfg.rollIndex = SeakeepingDOF::findModeIndex(SeakeepingCfg, MODE_ROLL);
    cfg.rollVisc = rollVisc_;

    cummins::CumminsTimeStepper stepper;
    stepper.init(std::move(cfg));

    Eigen::MatrixXd motionsLocal = Eigen::MatrixXd::Zero(N, D);

    for (int n = 0; n < N; ++n)
    {
        stepper.step(ctx.ExcitingForce.row(n).transpose());
        motionsLocal.row(n) = stepper.state().q.transpose();
    }

    SeakeepingDOF::scaleMotionsForOutput(
        SeakeepingCfg,
        CaseContextLite{ ctx.Amp, ctx.W, ctx.U, ctx.dt, ctx.we },
        motionsLocal);

    std::filesystem::create_directories(std::filesystem::path(filePath) / "motion");

    std::ofstream out(filePath + "motion/motion" + ctx.tag + ".csv");
    out << "t_nd";
    for (int k = 0; k < D; ++k) out << ",mode" << SeakeepingCfg.modes[k];
    out << "\n";

    for (int i = 0; i < motionsLocal.rows(); ++i)
    {
        out << (i + 1) * ctx.dt * std::sqrt(G / ShipCfg.Geometry.Length);
        for (int k = 0; k < D; ++k) out << "," << motionsLocal(i, k);
        out << "\n";
    }

    Fit fitRAO;
    auto names = SeakeepingDOF::modeNames(SeakeepingCfg);
    fitRAO.setdata(motionsLocal, ctx.dt, ctx.we, 1.0, names);
    fitRAO.run();

    std::vector<std::complex<double>> Hk;
    Hk.reserve(static_cast<std::size_t>(D));
    for (const auto& r : fitRAO.results())
    {
        Hk.emplace_back(
            r.amplitude * std::cos(r.phase),
            r.amplitude * std::sin(r.phase));
    }

    if (raoEnabled_)
        tab.setAt(ctx.iFn, ctx.iDir, ctx.iw, Hk);
}

void LinearCumminsTDGF::solveKernelCase_RK4(
    const CaseContext& ctx,
    const RadiationKernelData& kernel,
    RAO4DTable& tab)
{
    Eigen::MatrixXd Mphys, Bphys, C;
    SeakeepingDOF::buildSystemMatrices(ShipCfg, SeakeepingCfg, hs, Mphys, Bphys, C);

    const auto Konline = resampleKernelToDt(kernel, ctx.dt);

    const int D = SeakeepingCfg.DOF;
    const int N = static_cast<int>(ctx.ExcitingForce.rows());
    const int M = static_cast<int>(Konline.size()) - 1;
    const double dt = ctx.dt;

    if (D <= 0 || N <= 0 || dt <= 0.0)
        throw std::runtime_error("solveKernelCase_RK4: invalid D/N/dt.");

    const Eigen::MatrixXd Mtotal = Mphys + kernel.A_inf;
    const Eigen::MatrixXd Bconst = Bphys + kernel.B;
    const Eigen::MatrixXd Ctotal = C + kernel.C_prime;

    Eigen::FullPivLU<Eigen::MatrixXd> massSolver(Mtotal);
    if (!massSolver.isInvertible())
        throw std::runtime_error("solveKernelCase_RK4: Mtotal is singular.");

    motions = Eigen::MatrixXd::Zero(N, D);

    // 历史位移/速度，用于卷积记忆项
    /*std::vector<Eigen::VectorXd> qHist(static_cast<std::size_t>(N), Eigen::VectorXd::Zero(D));
    std::vector<Eigen::VectorXd> vHist(static_cast<std::size_t>(N), Eigen::VectorXd::Zero(D));*/

    Eigen::MatrixXd qHist = Eigen::MatrixXd::Zero(D, N);
    Eigen::MatrixXd vHist = Eigen::MatrixXd::Zero(D, N);

    Eigen::VectorXd q = Eigen::VectorXd::Zero(D);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(D);

    // ---------------------------
    // 激励力：直接从 ctx.ExcitingForce 取
    // frac = 0, 0.5, 1 时做线性插值，供 RK4 各 stage 使用
    // ---------------------------
    auto FextAt = [&](int n, double frac) -> Eigen::VectorXd
        {
            if (n <= 0)
                return ctx.ExcitingForce.row(0).transpose();

            if (n >= N - 1)
                return ctx.ExcitingForce.row(N - 1).transpose();

            const Eigen::VectorXd F0 = ctx.ExcitingForce.row(n).transpose();
            const Eigen::VectorXd F1 = ctx.ExcitingForce.row(n + 1).transpose();
            return (1.0 - frac) * F0 + frac * F1;
        };

    // ---------------------------
    // 横摇黏性阻尼
    // ---------------------------
    const int iR = SeakeepingDOF::findModeIndex(SeakeepingCfg, MODE_ROLL);

    auto rollForceAt = [&](const Eigen::VectorXd& vel) -> Eigen::VectorXd
        {
            Eigen::VectorXd fr = Eigen::VectorXd::Zero(D);
            if (iR >= 0)
                fr(iR) += rollVisc_.moment(vel(iR));
            return fr;
        };

    // ---------------------------
    // 记忆卷积项
    //
    // 在第 n -> n+1 步、某个 stage 速度 vStage 下，
    // 取：
    //   lag=1 用当前 stage 速度
    //   lag>=2 用已知历史速度
    //
    // mem = dt * sum_{lag>=1} K(lag*dt) * v(t_{n+1-lag})
    // ---------------------------
    //auto memoryAt = [&](int n, const Eigen::VectorXd& vStage) -> Eigen::VectorXd
    //{
    //    Eigen::VectorXd mem = Eigen::VectorXd::Zero(D);

    //    if (M <= 0)
    //        return mem;

    //    // lag = 1, 对应最近一段，用当前 stage 速度
    //    if (n >= 0 && M >= 1)
    //        mem.noalias() += dt * Konline[1] * vStage;

    //    // lag >= 2, 用已知历史速度
    //    const int maxLag = std::min(M, n + 1);
    //    for (int lag = 2; lag <= maxLag; ++lag)
    //    {
    //        const int idx = n - (lag - 2); // lag=2 -> n, lag=3 -> n-1, ...
    //        if (idx >= 0)
    //            mem.noalias() += dt * Konline[static_cast<std::size_t>(lag)]
    //            * vHist[static_cast<std::size_t>(idx)];
    //    }

    //    return mem;
    //};

    Eigen::VectorXd memHist = Eigen::VectorXd::Zero(D);
    auto memoryAt = [&](int n, const Eigen::VectorXd& vStage) -> Eigen::VectorXd {
        if (M <= 0) return Eigen::VectorXd::Zero(D);
        Eigen::VectorXd mem = memHist;           // 从预算结果出发
        if (n >= 0 && M >= 1)
            mem.noalias() += dt * Konline[1] * vStage;   // 只加 lag=1
        return mem;
        };


    // ---------------------------
    // 一阶状态方程：
    // y = [q; v]
    // ydot = [v; a]
    // a = M^{-1}[ F - Bv - Cq - mem + F_roll ]
    // ---------------------------
    auto rhsState = [&](int n,
        double frac,
        const Eigen::VectorXd& qStage,
        const Eigen::VectorXd& vStage) -> Eigen::VectorXd
        {
            Eigen::VectorXd rhs(2 * D);
            rhs.setZero();

            const Eigen::VectorXd Fext = FextAt(n, frac);
            const Eigen::VectorXd mem = memoryAt(n, vStage);
            const Eigen::VectorXd Froll = rollForceAt(vStage);

            const Eigen::VectorXd acc =
                massSolver.solve(Fext - Bconst * vStage - Ctotal * qStage - mem + Froll);
            //massSolver.solve(Fext - Bconst * vStage - Ctotal * qStage  + Froll);

            rhs.head(D) = vStage;
            rhs.tail(D) = acc;
            return rhs;
        };

    // 初值
    motions.row(0) = q.transpose();
    //qHist[0] = q;
    //vHist[0] = v;

    qHist.col(0) = q;
    vHist.col(0) = v;

    // n=0 时先用初始激励算一次加速度
    {
        const Eigen::VectorXd F0 = ctx.ExcitingForce.row(0).transpose();
        const Eigen::VectorXd a0 =
            massSolver.solve(F0 - Bconst * v - Ctotal * q + rollForceAt(v));
        (void)a0; // 这里不单独存，加速度由 RK4 rhs 自动使用
    }

    for (int n = 0; n < N - 1; ++n)
    {
        memHist = Eigen::VectorXd::Zero(D);
        {
            const int maxLag = std::min(M, n + 1);
            for (int lag = 2; lag <= maxLag; ++lag) {
                const int idx = n - (lag - 2);  // lag=2->n, lag=3->n-1, ...
                if (idx >= 0)
                    memHist.noalias() += dt * Konline[static_cast<size_t>(lag)]
                    * vHist.col(idx);
            }
        }


        Eigen::VectorXd y(2 * D);
        y.head(D) = q;
        y.tail(D) = v;

        const Eigen::VectorXd k1 = rhsState(n, 0.0, y.head(D), y.tail(D));

        const Eigen::VectorXd y2 = y + 0.5 * dt * k1;
        const Eigen::VectorXd k2 = rhsState(n, 0.5, y2.head(D), y2.tail(D));

        const Eigen::VectorXd y3 = y + 0.5 * dt * k2;
        const Eigen::VectorXd k3 = rhsState(n, 0.5, y3.head(D), y3.tail(D));

        const Eigen::VectorXd y4 = y + dt * k3;
        const Eigen::VectorXd k4 = rhsState(n, 1.0, y4.head(D), y4.tail(D));

        const Eigen::VectorXd yNext =
            y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

        q = yNext.head(D);
        v = yNext.tail(D);

        //qHist[static_cast<std::size_t>(n + 1)] = q;
        //vHist[static_cast<std::size_t>(n + 1)] = v;

        qHist.col(n + 1) = q;
        vHist.col(n + 1) = v;

        motions.row(n + 1) = q.transpose();
    }

    SeakeepingDOF::scaleMotionsForOutput(
        SeakeepingCfg,
        CaseContextLite{ ctx.Amp, ctx.W, ctx.U, ctx.dt, ctx.we },
        motions);

    std::filesystem::create_directories(std::filesystem::path(filePath) / "motion");

    std::ofstream out(filePath + "motion/motion" + ctx.tag + ".csv");
    out << "t_nd";
    for (int k = 0; k < D; ++k) out << ",mode" << SeakeepingCfg.modes[k];
    out << "\n";

    for (int i = 0; i < motions.rows(); ++i)
    {
        out << (i + 1) * ctx.dt * std::sqrt(G / ShipCfg.Geometry.Length);
        for (int k = 0; k < D; ++k) out << "," << motions(i, k);
        out << "\n";
    }

    Fit fitRAO;
    auto names = SeakeepingDOF::modeNames(SeakeepingCfg);
    fitRAO.setdata(motions, ctx.dt, ctx.we, 1.0, names);
    fitRAO.run();

    std::vector<std::complex<double>> Hk;
    Hk.reserve(static_cast<std::size_t>(D));
    for (const auto& r : fitRAO.results())
    {
        Hk.emplace_back(
            r.amplitude * std::cos(r.phase),
            r.amplitude * std::sin(r.phase));
    }

    if (raoEnabled_)
        tab.setAt(ctx.iFn, ctx.iDir, ctx.iw, Hk);
}


int LinearCumminsTDGF::findExactIndex(
    const std::vector<double>& xs, double x, double tol)
{
    auto it = std::lower_bound(xs.begin(), xs.end(), x);
    if (it != xs.end() && std::abs(*it - x) <= tol)
        return static_cast<int>(std::distance(xs.begin(), it));

    if (it != xs.begin())
    {
        auto it2 = it - 1;
        if (std::abs(*it2 - x) <= tol)
            return static_cast<int>(std::distance(xs.begin(), it2));
    }

    throw std::runtime_error("value not on axis.");
}



void LinearCumminsTDGF::runFreeRollDecay(const std::string& casePath)
{
    std::cout << "==== Free roll decay mode ====\n";

    // 1) 先建立横摇阻尼
    double Ieff = ShipCfg.Mass.Ixx; // 默认直接用船体横摇转动惯量

    if (SeakeepingCfg.RollDamping.enabled &&
        SeakeepingCfg.RollDamping.mode == RollDampingMode::FromDecayCsv)
    {
        auto samples = RollDampingBuilder::readDecayCsv(
            filePath + SeakeepingCfg.RollDamping.decay.csvPath,
            SeakeepingCfg.RollDamping.decay.angleInDeg
        );

        auto result = RollDampingBuilder::identifyFromDecay(
            samples,
            ShipCfg.Mass.Mass,
            ShipCfg.Mass.GM,
            SeakeepingCfg.RollDamping.decay.polyOrder,
            SeakeepingCfg.RollDamping.decay.minPeakGap,
            SeakeepingCfg.RollDamping.decay.IeffOverride,
            SeakeepingCfg.RollDamping.decay.refine,
            SeakeepingCfg.RollDamping.decay.refineSkipFirstPeaks
        );

        rollVisc_ = result.damping;
        Ieff = result.Ieff;

        std::cout << "Roll damping identified from CSV:\n"
            << "  omega_n   = " << result.poly.omega_n << "\n"
            << "  T_n       = " << result.poly.T_n << "\n"
            << "  Ieff      = " << result.Ieff << "\n"
            << "  B44_lin   = " << rollVisc_.B44_lin << "\n"
            << "  B44_quad  = " << rollVisc_.B44_quad << "\n"
            << "  B44_cube  = " << rollVisc_.B44_cube << "\n";
    }
    else
    {
        rollVisc_ = RollDampingBuilder::build(
            casePath,
            SeakeepingCfg.RollDamping,
            ShipCfg.Mass.Mass,
            ShipCfg.Mass.GM
        );

        std::cout << "Roll damping loaded directly:\n"
            << "  Ieff      = " << Ieff << "\n"
            << "  B44_lin   = " << rollVisc_.B44_lin << "\n"
            << "  B44_quad  = " << rollVisc_.B44_quad << "\n"
            << "  B44_cube  = " << rollVisc_.B44_cube << "\n";
    }

    // 2) 线性恢复矩系数
    // C44 = m * g * GM
    const double C44 = ShipCfg.Mass.Mass * G * ShipCfg.Mass.GM;

    // 3) 初始条件
    const double dt = SeakeepingCfg.FreeRollDecay.dt;
    const double duration = SeakeepingCfg.FreeRollDecay.duration;

    double phi = SeakeepingCfg.FreeRollDecay.phi0_deg * PI / 180.0;
    double phidot = SeakeepingCfg.FreeRollDecay.phidot0_deg_s * PI / 180.0;

    const int nStep = static_cast<int>(duration / dt) + 1;

    const double Tn = 2.0 * PI * std::sqrt(Ieff / C44);
    const double omega_n = std::sqrt(C44 / Ieff);

    std::cout << "Free roll decay setup:\n"
        << "  phi0 [deg]   = " << SeakeepingCfg.FreeRollDecay.phi0_deg << "\n"
        << "  phidot0 [deg/s] = " << SeakeepingCfg.FreeRollDecay.phidot0_deg_s << "\n"
        << "  dt [s]       = " << dt << "\n"
        << "  duration [s] = " << duration << "\n"
        << "  C44          = " << C44 << "\n"
        << "  Ieff         = " << Ieff << "\n"
        << "  omega_n      = " << omega_n << "\n"
        << "  Tn [s]       = " << Tn << "\n";

    // 4) 输出文件
    std::string outFile = filePath + "free_roll_decay.csv";
    std::ofstream out(outFile);
    if (!out.is_open())
        throw std::runtime_error("cannot create free_roll_decay.csv");

    out << "t,phi_rad,phi_deg,phidot_rad_s,phidot_deg_s,M_visc\n";

    // 5) 定义 1DOF 横摇方程
    auto rhs = [&](double phi_now, double phidot_now)
        {
            const double Mvisc = rollVisc_.moment(phidot_now);
            const double phiddot = (Mvisc - C44 * phi_now) / Ieff;
            return std::pair<double, double>{ phidot_now, phiddot };
        };

    // 6) RK4 时间积分
    double t = 0.0;
    for (int i = 0; i < nStep; ++i)
    {
        const double Mvisc = rollVisc_.moment(phidot);

        out << t << ","
            << phi << ","
            << phi * 180.0 / PI << ","
            << phidot << ","
            << phidot * 180.0 / PI << ","
            << Mvisc << "\n";

        auto [k1_phi, k1_phidot] = rhs(phi, phidot);
        auto [k2_phi, k2_phidot] = rhs(
            phi + 0.5 * dt * k1_phi,
            phidot + 0.5 * dt * k1_phidot
        );
        auto [k3_phi, k3_phidot] = rhs(
            phi + 0.5 * dt * k2_phi,
            phidot + 0.5 * dt * k2_phidot
        );
        auto [k4_phi, k4_phidot] = rhs(
            phi + dt * k3_phi,
            phidot + dt * k3_phidot
        );

        phi += dt * (k1_phi + 2.0 * k2_phi + 2.0 * k3_phi + k4_phi) / 6.0;
        phidot += dt * (k1_phidot + 2.0 * k2_phidot + 2.0 * k3_phidot + k4_phidot) / 6.0;

        t += dt;
    }
    out.close();
    std::cout << "Free roll decay finished. Output: " << outFile << std::endl;
}


void LinearCumminsTDGF::TimeToFrequency(const std::string& casePat,
    double Fn, double beta, RadiationKernelData& kernel)
{
    std::vector<double> nada_L = { 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.0,
        1.1, 1.2, 1.3, 1.4, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0, 10.0, 30.0 };
    double L = ShipCfg.Geometry.Length;
    double displacement = ShipCfg.Geometry.Displacement;
    std::vector<double> omiga = nada_L;
    std::vector<double> omiga_e = nada_L;
    double k;
    double U = Fn * sqrt(G * L);
    const int N = nada_L.size();

    for (int i = 0; i < N; ++i)
    {
        k = 2.0 * PI / (nada_L[i] * L);
        omiga[i] = sqrt(k * G);
        omiga_e[i] = omiga[i] - U * k * cos(beta);
    }

    Eigen::MatrixXd& A_inf = kernel.A_inf;
    Eigen::MatrixXd& B = kernel.B;
    Eigen::MatrixXd& C = kernel.C_prime;
    std::vector<Eigen::MatrixXd>& Klag = kernel.Klag;
    double dt = kernel.dt;
    int DOF = kernel.DOF;
    std::vector<int> modes = kernel.modes;

    // 非均匀路径：Klag_times 存了真实节点时间；下方傅里叶积分要用真实时间和真实区间长。
    // 等距路径：Klag_times[i] = i*dt（buildKernelDirectTimeDomain 末尾兜底了），结果与之前一致。
    const bool haveTimes =
        kernel.Klag_times.size() == kernel.Klag.size() && !kernel.Klag_times.empty();

    //外层是遭遇频率数，内层是DOF*DOF的矩阵
    std::vector<Eigen::MatrixXd> A_added(N);
    std::vector<Eigen::MatrixXd> B_added(N);

    //无量纲化
    std::array<std::array<double, 6>, 6> A_non;
    std::array<std::array<double, 6>, 6> B_non;
    // 默认所有因子置 1.0（= 不做无量纲化，直接输出有量纲值）。
    // 注意：std::array::fill({1.0}) 只把每行“首元素”设 1、其余为 0，
    // 会导致横摇(mode 3) 等未显式赋值的因子为 0 → 除零 → CSV 里横摇列全是 inf。
    // 这里逐行 fill(1.0) 让所有元素默认 1.0：横摇列直接给出有量纲 A44/B44，
    // 而下面对 heave(2)/pitch(4) 的显式无量纲化行为完全不变。
    for (auto& r : A_non) r.fill(1.0);
    for (auto& r : B_non) r.fill(1.0);

    A_non[2][2] = rho * displacement;
    B_non[2][2] = rho * displacement * sqrt(G / L);

    A_non[4][2] = A_non[2][2] * L;
    B_non[4][2] = B_non[2][2] * L;

    A_non[4][4] = rho * displacement * L * L;
    B_non[4][4] = rho * displacement * L * L * sqrt(G / L);

    A_non[2][4] = A_non[4][4] / L;
    B_non[2][4] = B_non[4][4] / L;

    double omega_non = sqrt(G / L);

    double sum_A, sum_B;
    double tau, cos_term, sin_term;

    std::string outFile = filePath + "added_mass_damping.csv";
    std::ofstream out(outFile);

    out << "non_omiga," << "non_omiga_e,";
    for (int m = 0; m < DOF; ++m)
        for (int n = 0; n < DOF; ++n)
            out << "A_added" + std::to_string(modes[m]) + std::to_string(modes[n])
            << "," << "B_added" + std::to_string(modes[m]) + std::to_string(modes[n]) << ",";
    out << std::endl;

    // Per-frequency computation is independent; parallelise the heavy loop, then
    // serialise the file write to keep CSV row order stable.
#pragma omp parallel for schedule(static)
    for (int i = 0; i < N; ++i)
    {
        shipsim::EigenSingleThreadGuard _g;

        A_added[i] = Eigen::MatrixXd::Zero(DOF, DOF);
        B_added[i] = Eigen::MatrixXd::Zero(DOF, DOF);

        for (int m = 0; m < DOF; ++m)
        {
            const int mode_m = modes[m];
            for (int n = 0; n < DOF; ++n)
            {
                const int mode_n = modes[n];

                double sumA = 0.0;
                double sumB = 0.0;
                // 非均匀左矩形积分：每个 lag 用真实时间 + 真实前向区间长。
                // 等距路径下 tauL = lag*dt、w = dt，与原行为一致。
                for (size_t lag = 1; lag < Klag.size(); ++lag)
                {
                    const double tauL = haveTimes
                        ? kernel.Klag_times[lag]
                        : (lag * dt);
                    const double w = haveTimes
                        ? (kernel.Klag_times[lag] - kernel.Klag_times[lag - 1])
                        : dt;
                    const double cosT = std::cos(omiga_e[i] * tauL);
                    const double sinT = std::sin(omiga_e[i] * tauL);
                    sumA += Klag[lag](m, n) * sinT * w;
                    sumB += Klag[lag](m, n) * cosT * w;
                }

                A_added[i](m, n) = (A_inf(m, n) - sumA / omiga_e[i]
                    - C(m, n) / (omiga_e[i] * omiga_e[i])) / A_non[mode_m][mode_n];

                B_added[i](m, n) = (B(m, n) + sumB) / B_non[mode_m][mode_n];
            }
        }
    }

    for (int i = 0; i < N; ++i)
    {
        out << omiga[i] / omega_non << "," << omiga_e[i] / omega_non << ",";
        for (int m = 0; m < DOF; ++m)
            for (int n = 0; n < DOF; ++n)
                out << A_added[i](m, n) << "," << B_added[i](m, n) << ",";
        out << std::endl;
    }
    std::cout << "compute added mass and damping successfully!" << std::endl;
}


void LinearCumminsTDGF::postprocessRadiationKernel(
    double Fn,
    RadiationKernelData& kernel) const
{
    KernelPost::KernelScaleInfo scale;
    scale.L = ShipCfg.Geometry.Length;
    scale.displacement = ShipCfg.Geometry.Displacement;
    scale.rho = rho;
    scale.g = G;

    KernelPost::TailSplitOptions splitOpt;

    splitOpt.enabled = true;

    // 方案 A：
    // Klag -= Kinf
    // C_prime += Kinf
    splitOpt.addTailToCPrime = true;

    // 你的编号：2 = heave，4 = pitch
    splitOpt.targetModes = { 2, 4 };
    splitOpt.onlyTargetBlock = true;

    splitOpt.tailFraction = 0.20;
    splitOpt.minTailSamples = 30;

    splitOpt.minMeanToPeak = 0.05;
    splitOpt.maxOscToMean = 0.70;

    // K(2,4)/K(4,2) 实际上趋于 0，不需要强制 tail-split；只 K(4,4) 需要。
    // 让 tail-split 自动根据 meanToPeak 决定。

    splitOpt.verbose = true;

    KernelPost::TailSplitWorkflowOptions workflowOpt;

    workflowOpt.enabled = true;

    // 调试阶段建议 true；正式批量计算可改 false。
    workflowOpt.writeBeforeHistory = true;
    workflowOpt.writeAfterHistory = true;

    workflowOpt.writeDimensional = true;
    workflowOpt.writeNondimensional = true;

    workflowOpt.outputDir =
        (std::filesystem::path(filePath) / "kernel_cache").string();

    workflowOpt.tag = "Fn" + keyDouble(Fn);

    KernelPost::applyTailSplitWorkflowInPlace(
        kernel,
        scale,
        splitOpt,
        workflowOpt);

    if (kRadiationKernelSgSmooth && !kernel.Klag.empty())
    {
        impulse_sg::Options sgOpt;
        std::vector<Eigen::MatrixXd> kBefore;
        if (kRadiationKernelSgCompareCsv)
            kBefore = kernel.Klag;
        impulse_sg::smoothMatrixLagHistoryInPlace(kernel.Klag, sgOpt);
        if (!kernel.Klag.empty())
            kernel.K0 = kernel.Klag[0];
        if (kRadiationKernelSgCompareCsv && !kBefore.empty())
        {
            const std::string csvDir =
                (std::filesystem::path(workflowOpt.outputDir) / "kernel_sg_compare").string();
            const std::string csvPath =
                (std::filesystem::path(csvDir) / ("Klag_compare_" + workflowOpt.tag + ".csv")).string();
            const std::vector<std::pair<int, int>> pairs = {
                {2, 2}, {2, 4}, {4, 4}, {4, 2} }; // 0-based: heave=2, pitch=4
            if (impulse_sg::writeMatrixLagCompareCsv(csvPath, kernel.dt, pairs, kBefore, kernel.Klag))
                std::cout << "[ImpulseKernelSG] wrote " << csvPath << "\n";
        }
    }

    // Final diagnostic dump after tail-split + SG smoothing.
    // Compare against the "raw_*" files to see what postprocess actually did.
    writeKernelDiagnostics(Fn, kernel, "post", kernel.dt);
}

// =====================================================================
// 备选后处理策略（与 postprocessRadiationKernel 二选一）
// ---------------------------------------------------------------------
// 动机：原 tail-split 把 K(t) 末端的“常数” Kinf 扣掉并塞进 C_prime，
//       这对 heave/pitch 有两个副作用：(a) 改了恢复力 -> 共振频率漂移；
//       (b) Kinf 若来自“尚未收敛的尾巴”会被高估 -> 误删低频阻尼 -> 共振峰偏大。
//
// 本策略不再扣常数、也不动 A_inf/B/C_prime，而是：
//   (1) 保留原始 [0, T] 段 K(t) 不变；
//   (2) 把时间轴向后延长 extendFactor 倍（题目要求 2 倍）；
//   (3) 延长段补上 K(t) 末端“趋近的常数” Kc。
//
// Kc 的取法（“直接复制最后一个点” vs “尾段识别一个常数”）：
//   这里用“尾段均值 + 最小二乘斜率”来估 Kc，比直接复制最后一个点更稳，原因：
//     * 最后一个样本来自 differentiateHistory 的单边端点差分，精度最低、最噪；
//     * 若尾巴还在振荡，单点只是抓到某个相位的瞬时值，会在拼接处造成台阶，
//       台阶会在 B(ω)/A(ω) 里引入虚假高频纹波；
//     * 尾段均值把振荡平均成它的“直流电平”，正是我们想要的渐近常数；
//     * 顺带算出斜率 Kslope：|Kslope| 偏大说明核根本没收敛，常数外推本身不可靠
//       （此时应优先加大离线核记忆 tMem，而不是靠补常数）。
// =====================================================================
void LinearCumminsTDGF::postprocessRadiationKernel2(
    double Fn,
    RadiationKernelData& kernel) const
{
    if (kernel.Klag.empty()) return;

    // -------------------- 可调参数（切换/调参只改这里）--------------------
    const double extendFactor = 2.0;    // 时间轴延长为原来的倍数（题目：2 倍）
    const double tailWinFrac  = 0.25;   // 用末尾这一比例的样本估计 Kc / 斜率
    const bool   decayToZero  = false;  // false: 尾部保持常数 Kc（题目要求）
                                        // true : 尾部把 Kc 平滑衰减到 0（更符合物理，建议对比）
    const bool   doSgSmooth   = kRadiationKernelSgSmooth; // 是否沿用原 SG 平滑（作用在原始段）
    // ---------------------------------------------------------------------

    const int    D  = (kernel.DOF > 0) ? kernel.DOF : SeakeepingCfg.DOF;
    const int    N  = static_cast<int>(kernel.Klag.size());      // 原样本数
    const double dt = (kernel.dt > 0.0) ? kernel.dt : dtShared_;

    // 时间轴兜底：等距时 t[i]=i*dt；非均匀时由 build 阶段填好。
    if (kernel.Klag_times.size() != kernel.Klag.size())
    {
        kernel.Klag_times.assign(static_cast<std::size_t>(N), 0.0);
        for (int i = 0; i < N; ++i)
            kernel.Klag_times[static_cast<std::size_t>(i)] = i * dt;
    }

    // ---- 可选：先对原始 K 做 SG 平滑（与原后处理保持一致，便于横向对比）----
    if (doSgSmooth && !kernel.Klag.empty())
    {
        impulse_sg::Options sgOpt;
        impulse_sg::smoothMatrixLagHistoryInPlace(kernel.Klag, sgOpt);
        if (!kernel.Klag.empty())
            kernel.K0 = kernel.Klag.front();
    }

    // ---- 估计“趋近常数” Kc：尾段窗口最小二乘 a + b*t（逐 (i,j) 元素）----
    int win = std::max(3, static_cast<int>(std::lround(tailWinFrac * N)));
    win = std::min(win, N);
    const int wi0 = N - win;

    // 标量量（横轴对所有 (i,j) 相同），矩阵量（每个元素各自累加）
    double Sx = 0.0, Sxx = 0.0;
    Eigen::MatrixXd Sy  = Eigen::MatrixXd::Zero(D, D);
    Eigen::MatrixXd Sxy = Eigen::MatrixXd::Zero(D, D);
    for (int k = 0; k < win; ++k)
    {
        const int    idx = wi0 + k;
        const double x   = kernel.Klag_times[static_cast<std::size_t>(idx)];
        const Eigen::MatrixXd& y = kernel.Klag[static_cast<std::size_t>(idx)];
        Sx  += x;
        Sxx += x * x;
        Sy  += y;
        Sxy += x * y;
    }
    const double n      = static_cast<double>(win);
    const double xbar   = Sx / n;
    const double denom  = Sxx - Sx * xbar;           // = Σ(x-xbar)^2
    const Eigen::MatrixXd ybar = Sy / n;             // 窗口均值（直流电平）
    Eigen::MatrixXd Kslope = Eigen::MatrixXd::Zero(D, D);
    if (std::abs(denom) > 1e-30)
        Kslope = (Sxy - Sx * ybar) / denom;          // 最小二乘斜率 b
    const Eigen::MatrixXd Kc = ybar;                 // 取直流电平作为渐近常数

    // ---- 组装延长后的时间轴与核 ----
    const int    Nnew = std::max(N + 1, static_cast<int>(std::lround(N * extendFactor)));
    const int    ext  = Nnew - N;                    // 新增样本数
    const double dtTail = (N >= 2)
        ? (kernel.Klag_times[static_cast<std::size_t>(N - 1)]
            - kernel.Klag_times[static_cast<std::size_t>(N - 2)])
        : dt;                                        // 延长段步长（兼容非均匀网格）

    const Eigen::MatrixXd Klast = kernel.Klag[static_cast<std::size_t>(N - 1)];

    // 拼接处过渡：用 raised-cosine 把“原末端值 Klast”平滑过渡到“常数 Kc”，
    // 消除台阶（台阶会在 B(ω) 里产生虚假纹波）。过渡长度取延长段的 ~1/5。
    const int blendLen = std::min(ext, std::max(5, ext / 5));

    std::vector<Eigen::MatrixXd> newKlag(static_cast<std::size_t>(Nnew));
    std::vector<double>          newTimes(static_cast<std::size_t>(Nnew));

    for (int m = 0; m < N; ++m)
    {
        newKlag[static_cast<std::size_t>(m)]  = kernel.Klag[static_cast<std::size_t>(m)];
        newTimes[static_cast<std::size_t>(m)] = kernel.Klag_times[static_cast<std::size_t>(m)];
    }

    const double tJoin = kernel.Klag_times[static_cast<std::size_t>(N - 1)];
    for (int k = 0; k < ext; ++k)
    {
        const int m = N + k;
        newTimes[static_cast<std::size_t>(m)] = tJoin + (k + 1) * dtTail;

        // 1) 先得到“平台值” plateau：blend 区间内 Klast->Kc，之后保持 Kc。
        Eigen::MatrixXd plateau;
        if (k < blendLen)
        {
            const double w = 0.5 * (1.0 - std::cos(PI * (k + 1) / blendLen)); // 0->1
            plateau = (1.0 - w) * Klast + w * Kc;
        }
        else
        {
            plateau = Kc;
        }

        // 2) 可选：把平台值平滑衰减到 0（更符合 retardation function 物理）。
        //    衰减窗从 blend 结束处开始，到延长段末端降到 0（raised-cosine）。
        if (decayToZero)
        {
            double g = 1.0;
            const int decStart = blendLen;
            const int decLen    = std::max(1, ext - decStart);
            if (k >= decStart)
            {
                const double s = static_cast<double>(k - decStart) / decLen; // 0->1
                g = 0.5 * (1.0 + std::cos(PI * std::min(1.0, s)));           // 1->0
            }
            plateau *= g;
        }

        newKlag[static_cast<std::size_t>(m)] = plateau;
    }

    // ---- 写回 kernel（A_inf / B / C_prime 一律不动）----
    kernel.Klag       = std::move(newKlag);
    kernel.Klag_times = std::move(newTimes);
    kernel.K0         = kernel.Klag.front();
    kernel.TG         = static_cast<int>(kernel.Klag.size()) - 1;

    // ---- 诊断输出 ----
    const double maxKc    = Kc.cwiseAbs().maxCoeff();
    const double maxSlope = Kslope.cwiseAbs().maxCoeff();
    std::cout << "[postprocessRadiationKernel2] Fn=" << Fn
              << "  N:" << N << "->" << kernel.TG + 1
              << "  tMem:" << tJoin << "->" << kernel.Klag_times.back() << "s"
              << "  max|Kc|=" << maxKc
              << "  max|Kslope|=" << maxSlope
              << (decayToZero ? "  [tail=decay->0]" : "  [tail=const]")
              << "\n";
    if (maxSlope * (kernel.Klag_times.back() - tJoin) > 0.10 * (maxKc + 1e-30))
    {
        std::cout << "  [warn] 尾段仍有明显斜率：核可能未收敛，补常数不可靠，"
                     "建议优先增大离线核记忆 tMem（setupOfflineKernelGridForFn）。\n";
    }

    writeKernelDiagnostics(Fn, kernel, "post2", dt);
}

