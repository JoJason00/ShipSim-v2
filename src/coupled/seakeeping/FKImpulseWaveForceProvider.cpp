#include "FKImpulseWaveForceProvider.h"
#include "../common/CoupledMath.h"
#include "../../const/Const.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace
{
    std::vector<double> splitCsvRow(const std::string& line)
    {
        std::vector<double> v;
        v.reserve(20);
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
        {
            if (!cell.empty())
                v.push_back(std::stod(cell));
        }
        return v;
    }

    struct GroupKey
    {
        WaveForceRegion region = WaveForceRegion::Head;
        double          lookupDir = 0.0;
        bool            mirrorPS = false;

        bool operator<(const GroupKey& o) const
        {
            if (region != o.region) return region < o.region;
            if (mirrorPS != o.mirrorPS) return mirrorPS < o.mirrorPS;
            return lookupDir < o.lookupDir - 1e-9;
        }
    };
}

FKImpulseWaveForceProvider::FKImpulseWaveForceProvider(
    const SeakeepingConfig& skCfg,
    const std::string& casePath,
    const std::string& shipName,
    const HydroRefreshConfig& refreshCfg,
    FKImpulsePart part)
    : skCfg_(skCfg),
      casePath_(casePath),
      shipName_(shipName),
      refreshCfg_(refreshCfg),
      part_(part)
{
    scanCaseDir();
    if (kernels_.empty())
        throw std::runtime_error(
            "FKImpulseWaveForceProvider: no fkImpulse_*.csv found in "
            + casePath_ + "fkImpulse/");

    const char* partTag = (part_ == FKImpulsePart::DfOnly) ? "df-only"
                        : (part_ == FKImpulsePart::FkOnly) ? "fk-only"
                        : "fk+df";
    std::cout << "[FKImpulse] loaded " << kernels_.size()
              << " kernel(s) from " << casePath_ << "fkImpulse/"
              << "  part=" << partTag << "\n";
}

void FKImpulseWaveForceProvider::scanCaseDir()
{
    const auto dir = std::filesystem::path(casePath_) / "fkImpulse";
    if (!std::filesystem::exists(dir)) return;

    for (const auto& de : std::filesystem::directory_iterator(dir))
    {
        if (!de.is_regular_file())                  continue;
        if (de.path().extension() != ".csv")        continue;
        if (de.path().filename().string().rfind("fkImpulse_", 0) != 0) continue;

        FileKernel k;
        if (loadFile(de.path().string(), k))
        {
            // 每个 library 文件一行的明细打印会刷屏，关掉。汇总行见构造函数末尾。
            // std::cout << "[FKImpulse] read library " << de.path().filename().string()
            //           << "  dt=" << k.dtSrc
            //           << "  tMot=" << k.tMot
            //           << "  full_rows=" << (2 * k.tMot + 1)
            //           << "  causal_rows=" << k.fk.rows() << "\n";
            kernels_.push_back(std::move(k));
        }
    }
}

bool FKImpulseWaveForceProvider::loadFile(const std::string& path, FileKernel& out)
{
    std::ifstream in(path);
    if (!in.is_open()) return false;

    int tMot = -1, rows = -1;
    double Fn = 0.0, dirRad = 0.0, dt = 0.0;
    WaveForceRegion region = WaveForceRegion::Head;

    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty()) continue;

        if (line[0] == '#')
        {
            const auto comma = line.find(',');
            if (comma == std::string::npos) continue;
            int keyStart = 1;
            while (keyStart < static_cast<int>(comma)
                   && std::isspace(static_cast<unsigned char>(line[keyStart])))
                ++keyStart;
            const std::string key = line.substr(keyStart, comma - keyStart);
            const std::string val = line.substr(comma + 1);
            if      (key == "Fn")     Fn     = std::stod(val);
            else if (key == "Dir")    dirRad = std::stod(val);
            else if (key == "Region") region = wave_force_region::fromTag(val);
            else if (key == "dt")     dt     = std::stod(val);
            else if (key == "tMot")   tMot   = std::stoi(val);
            else if (key == "rows")   rows   = std::stoi(val);
            continue;
        }

        if (line.rfind("i,tau", 0) == 0) break;
    }

    if (tMot <= 0 || rows != 2 * tMot + 1 || dt <= 0.0)
        return false;

    const int causal = tMot + 1;
    out.Fn     = Fn;
    out.dirRad = dirRad;
    out.region = region;
    out.dtSrc  = dt;
    out.tMot   = tMot;
    out.fk     = Eigen::MatrixXd::Zero(causal, 6);
    out.df     = Eigen::MatrixXd::Zero(causal, 6);

    while (std::getline(in, line))
    {
        if (line.empty()) continue;

        const auto v = splitCsvRow(line);
        if (v.size() < 14) continue;

        const int i = static_cast<int>(v[0]);
        if (i < tMot || i > 2 * tMot) continue;

        const int r = i - tMot;
        for (int j = 0; j < 6; ++j)
        {
            out.fk(r, j) = v[2 + j];
            out.df(r, j) = v[8 + j];
        }
    }
    return true;
}

const FKImpulseWaveForceProvider::FileKernel*
FKImpulseWaveForceProvider::findNearest(double Fn, double betaRel,
                                        WaveForceRegion region) const
{
    if (kernels_.empty()) return nullptr;

    const FileKernel* best = nullptr;
    double bestMetric = std::numeric_limits<double>::infinity();

    auto scan = [&](bool matchRegion)
    {
        for (const auto& k : kernels_)
        {
            if (matchRegion && k.region != region) continue;

            const double dFn   = std::abs(k.Fn - Fn);
            const double dBeta = std::abs(coupled::wrapPiSigned(k.dirRad - betaRel));
            const double m     = dFn + dBeta;
            if (m < bestMetric) { bestMetric = m; best = &k; }
        }
    };

    scan(true);
    if (!best) scan(false);
    return best;
}

int FKImpulseWaveForceProvider::buildGroupKernel(
    const FileKernel& src, double dtFast, bool mirrorPS,
    Eigen::MatrixXd& Kout) const
{
    const int D = skCfg_.DOF;
    const double tMax = src.dtSrc * src.tMot;
    const int nLag = static_cast<int>(std::ceil(tMax / dtFast)) + 1;

    Kout.setZero(nLag, D);

    auto impulseAt = [&](int row, int mode) -> double
    {
        switch (part_)
        {
        case FKImpulsePart::FkOnly:  return src.fk(row, mode);
        case FKImpulsePart::DfOnly:  return src.df(row, mode);
        case FKImpulsePart::FkAndDf:
        default:                     return src.fk(row, mode) + src.df(row, mode);
        }
    };

    for (int k = 0; k < nLag; ++k)
    {
        const double tau = k * dtFast;
        if (tau >= tMax)
        {
            for (int j = 0; j < D; ++j)
            {
                const int mode = skCfg_.modes[j];
                if (mode < 0 || mode >= 6) continue;
                double val = impulseAt(src.tMot, mode);
                if (mirrorPS && (mode == 1 || mode == 3 || mode == 5))
                    val = -val;
                Kout(k, j) = val;
            }
            continue;
        }

        const double x  = tau / src.dtSrc;
        const int    i0 = static_cast<int>(std::floor(x));
        const int    i1 = std::min(src.tMot, i0 + 1);
        const double a  = x - i0;

        for (int j = 0; j < D; ++j)
        {
            const int mode = skCfg_.modes[j];
            if (mode < 0 || mode >= 6) continue;
            double val = (1.0 - a) * impulseAt(i0, mode) + a * impulseAt(i1, mode);
            if (mirrorPS && (mode == 1 || mode == 3 || mode == 5))
                val = -val;
            Kout(k, j) = val;
        }
    }

    return nLag;
}

std::size_t FKImpulseWaveForceProvider::groupsSignature(
    const CoupledSlowState3DOF& slow0,
    const CoupledEncounterState& enc0)
{
    std::hash<double> hD;
    std::hash<int>    hI;
    std::size_t sig = hI(static_cast<int>(enc0.comps.size()));

    auto mix = [&](std::size_t v)
    {
        sig ^= v + 0x9e3779b97f4a7c15ULL + (sig << 6) + (sig >> 2);
    };

    for (const auto& w : enc0.comps)
    {
        const double betaRel = coupled::wrapPiSigned(w.theta - slow0.psi);
        double lookupDir = 0.0;
        bool mirrorPS = false;
        double betaRelDeg = 0.0;
        double betaQueryDeg = 0.0;
        coupled::waveHeadingLookup(betaRel, lookupDir, mirrorPS, betaRelDeg, betaQueryDeg);
        const auto region =
            wave_force_region::classify(slow0.U, betaRel, w.omega);

        mix(hD(w.omega));
        mix(hD(lookupDir));
        mix(hI(static_cast<int>(region)));
        mix(hI(mirrorPS ? 1 : 0));
    }
    return sig;
}

bool FKImpulseWaveForceProvider::needGroupRefresh(
    const CoupledSlowState3DOF& slow0,
    const CoupledEncounterState& enc0,
    double dtFast) const
{
    if (!groupsValid_)                              return true;
    if (std::abs(dtFast - dtCached_) > 1e-12)      return true;
    if (groupsSignature(slow0, enc0) != groupsSig_) return true;

    const double FnRef = std::max(std::abs(FnCached_), 1e-3);
    if (std::abs(slow0.Fn - FnCached_) > refreshCfg_.speedTolRatio * FnRef)
        return true;

    return false;
}

bool FKImpulseWaveForceProvider::needsMidWindowKernelRefresh(
    const CoupledSlowState3DOF& slow,
    const CoupledEncounterState& enc) const
{
    if (!groupsValid_) return false;

    const double dBetaDeg =
        std::abs(coupled::wrapPiSigned(enc.betaRel - betaRebuildCached_)) * 180.0 / PI;
    if (dBetaDeg > refreshCfg_.betaTolDeg)
        return true;

    const double dPsiDeg =
        std::abs(coupled::wrapPiSigned(slow.psi - psiRebuildCached_)) * 180.0 / PI;
    if (dPsiDeg > refreshCfg_.betaTolDeg)
        return true;

    const double FnRef = std::max(std::abs(FnCached_), 1e-3);
    if (std::abs(slow.Fn - FnCached_) > refreshCfg_.speedTolRatio * FnRef)
        return true;

    return false;
}

void FKImpulseWaveForceProvider::rebuildGroups(
    const CoupledSlowState3DOF& slow0,
    const CoupledEncounterState& enc0,
    double dtFast)
{
    // Crossfade setup: hold the currently-active groups as `oldGroups_` and
    // schedule a linear blend over the next N fast steps. The η history is
    // shared by both old and new convolutions (we never clear etaHist_),
    // so the crossfade is consistent in the past direction.
    // N is chosen to span ≈1 wave period at the highest expected mode
    // frequency, which is enough to suppress the step's broadband content
    // around any plausible ω_n. Falls back to a small constant if dtFast
    // is unset (very first call).
    if (!groups_.empty() && groupsValid_)
    {
        oldGroups_ = std::move(groups_);
        const double blendSec = 1.0; // ~1 s ≈ one wave period at ω≈6 rad/s
        const double dtUse = (dtFast > 1.0e-9 ? dtFast
                              : (dtCached_ > 1.0e-9 ? dtCached_ : 0.02));
        blendStepsTotal_ = std::max(4, static_cast<int>(std::lround(blendSec / dtUse)));
        blendStepsLeft_  = blendStepsTotal_;
    }
    else
    {
        oldGroups_.clear();
        blendStepsLeft_  = 0;
        blendStepsTotal_ = 0;
    }

    groups_.clear();
    memoryLen_ = 0;

    std::vector<WaveComponent> comps = enc0.comps;
    if (comps.empty())
    {
        const double thetaAbs = slow0.psi + enc0.betaRel;
        comps.push_back(WaveComponent{ enc0.waveAmp, enc0.omega, thetaAbs, 0.0 });
    }

    std::map<GroupKey, std::size_t> index;
    for (const auto& w : comps)
    {
        const double betaRel = coupled::wrapPiSigned(w.theta - slow0.psi);
        double lookupDir = 0.0;
        bool mirrorPS = false;
        double betaRelDeg = 0.0;
        double betaQueryDeg = 0.0;
        coupled::waveHeadingLookup(betaRel, lookupDir, mirrorPS, betaRelDeg, betaQueryDeg);

        const GroupKey key{
            wave_force_region::classify(slow0.U, betaRel, w.omega),
            lookupDir,
            mirrorPS
        };

        auto it = index.find(key);
        if (it == index.end())
        {
            index[key] = groups_.size();
            groups_.push_back(Group{});
            groups_.back().dt = dtFast;
        }
        groups_[index[key]].comps.push_back(w);
    }

    for (const auto& kv : index)
    {
        const GroupKey& key = kv.first;
        Group& g = groups_[kv.second];

        const FileKernel* lib = findNearest(slow0.Fn, key.lookupDir, key.region);
        if (!lib)
        {
            std::cerr << "[FKImpulse] no kernel for Fn=" << slow0.Fn
                      << " Dir=" << (key.lookupDir * 180.0 / PI) << "deg"
                      << " Reg=" << wave_force_region::tag(key.region) << "\n";
            continue;
        }

        const int nLag = buildGroupKernel(*lib, dtFast, key.mirrorPS, g.K);
        g.dt = dtFast;
        memoryLen_ = std::max(memoryLen_, nLag);

        // 每个 group 每次 rebuild 都会刷屏（仿真过程中频繁发生），关掉。
        // std::cout << "[FKImpulse] group Reg=" << wave_force_region::tag(key.region)
        //           << " lookupDir=" << (key.lookupDir * 180.0 / PI) << "deg"
        //           << (key.mirrorPS ? " [PS-mirror]" : "")
        //           << " comps=" << g.comps.size()
        //           << " lib Fn=" << lib->Fn
        //           << " lib Dir=" << (lib->dirRad * 180.0 / PI) << "deg"
        //           << " nLag=" << nLag << "\n";
    }

    groups_.erase(
        std::remove_if(groups_.begin(), groups_.end(),
                       [](const Group& g) { return g.K.rows() == 0; }),
        groups_.end());

    coupled::resizeHistoryPreservingEnd(etaHist_, memoryLen_);

    FnCached_           = slow0.Fn;
    betaRebuildCached_  = enc0.betaRel;
    psiRebuildCached_   = slow0.psi;
}

double FKImpulseWaveForceProvider::componentElevation(
    double x, double y, double t, const WaveComponent& w)
{
    const double k = w.omega * w.omega / G;
    const double phase =
        k * (x * std::cos(w.theta) + y * std::sin(w.theta))
        - w.omega * t + w.eps;
    return w.a * std::cos(phase);
}

double FKImpulseWaveForceProvider::etaAllComponents(
    const std::vector<WaveComponent>& comps,
    double x, double y, double t) const
{
    double eta = 0.0;
    for (const auto& w : comps)
        eta += componentElevation(x, y, t, w);
    return eta;
}

double FKImpulseWaveForceProvider::etaAtGroup(const Group& g, double tWave) const
{
    const CoupledSlowState3DOF s = slowAtTime(tWave);
    return etaAllComponents(g.comps, s.xe, s.ye, tWave);
}

CoupledSlowState3DOF FKImpulseWaveForceProvider::interpSlow(
    const CoupledSlowState3DOF& a,
    const CoupledSlowState3DOF& b,
    double alpha)
{
    CoupledSlowState3DOF s{};
    const double om = 1.0 - alpha;

    s.t     = om * a.t     + alpha * b.t;
    s.u     = om * a.u     + alpha * b.u;
    s.v     = om * a.v     + alpha * b.v;
    s.r     = om * a.r     + alpha * b.r;
    s.xe    = om * a.xe    + alpha * b.xe;
    s.ye    = om * a.ye    + alpha * b.ye;
    s.psi   = om * a.psi   + alpha * b.psi;
    s.delta = om * a.delta + alpha * b.delta;
    s.U     = om * a.U     + alpha * b.U;
    s.betaDrift  = om * a.betaDrift  + alpha * b.betaDrift;
    s.Fn         = om * a.Fn         + alpha * b.Fn;
    s.rDotRadPerS2 = om * a.rDotRadPerS2 + alpha * b.rDotRadPerS2;
    return s;
}

CoupledSlowState3DOF FKImpulseWaveForceProvider::extrapolateSlow(
    const CoupledSlowState3DOF& ref,
    double tRef,
    double tQuery)
{
    CoupledSlowState3DOF s = ref;
    s.t = tQuery;

    const double dt = tQuery - tRef;
    const double c  = std::cos(ref.psi);
    const double sn = std::sin(ref.psi);

    s.xe  = ref.xe + dt * (ref.u * c - ref.v * sn);
    s.ye  = ref.ye + dt * (ref.u * sn + ref.v * c);
    s.psi = ref.psi + dt * ref.r;
    return s;
}

CoupledSlowState3DOF FKImpulseWaveForceProvider::slowAtTime(double tQuery) const
{
    if (!trajValid_ || std::abs(trajT1_ - trajT0_) < 1e-14)
        return trajSlow0_;

    if (tQuery <= trajT0_)
        return extrapolateSlow(trajSlow0_, trajT0_, tQuery);
    if (tQuery >= trajT1_)
        return extrapolateSlow(trajSlow1_, trajT1_, tQuery);

    const double alpha = (tQuery - trajT0_) / (trajT1_ - trajT0_);
    return interpSlow(trajSlow0_, trajSlow1_, alpha);
}

void FKImpulseWaveForceProvider::onWindowBegin(
    const CoupledSlowState3DOF& slow0,
    const CoupledEncounterState& enc0,
    double dtFast,
    const CoupledFastWindowInfo& win)
{
    trajValid_ = false;
    if (win.slowEnd != nullptr)
    {
        trajSlow0_ = slow0;
        trajSlow1_ = *win.slowEnd;
        trajT0_    = win.t0;
        trajT1_    = win.t1;
        trajValid_ = true;
    }

    if (!needGroupRefresh(slow0, enc0, dtFast))
        return;

    rebuildGroups(slow0, enc0, dtFast);

    dtCached_    = dtFast;
    groupsSig_   = groupsSignature(slow0, enc0);
    groupsValid_ = !groups_.empty();

    // 每次窗口刷新都会打印，仿真过程中持续刷屏，关掉。
    // std::cout << "[FKImpulse] window refresh: Fn=" << slow0.Fn
    //           << " groups=" << groups_.size()
    //           << " memoryLen=" << memoryLen_
    //           << " dtFast=" << dtFast << "\n";
}

void FKImpulseWaveForceProvider::setEtaHistory(const std::vector<double>& hist)
{
    etaHist_ = hist;
    if (memoryLen_ > 0)
        coupled::resizeHistoryPreservingEnd(etaHist_, memoryLen_);
}

std::vector<double> FKImpulseWaveForceProvider::etaHistory() const
{
    return etaHist_;
}

Eigen::VectorXd FKImpulseWaveForceProvider::evaluate(
    const CoupledSlowState3DOF& slow,
    double t,
    const CoupledEncounterState& enc,
    Eigen::RowVectorXd* force6)
{
    const int D = skCfg_.DOF;
    Eigen::VectorXd F = Eigen::VectorXd::Zero(D);

    if (!groupsValid_ || groups_.empty())
        return F;

    // Turning / Z-track: β_rel changes quickly; refresh kernels with current ψ, not
    // only at the slow-window boundary (otherwise the wrong 0–180° file / mirror can
    // persist for up to dtSlow).
    if (needsMidWindowKernelRefresh(slow, enc))
        rebuildGroups(slow, enc, dtCached_);

    // Total η at CG in the earth-frame wave field (diagnostics / history hand-off).
    const CoupledSlowState3DOF sNow =
        trajValid_ ? slowAtTime(t) : slow;
    double etaNow = 0.0;
    if (!enc.comps.empty())
        etaNow = etaAllComponents(enc.comps, sNow.xe, sNow.ye, t);
    else
        etaNow = enc.waveAmp * std::cos(enc.phaseAtCG);
    coupled::pushRing(etaHist_, etaNow);

    // Causal convolution: η_g sampled at CG trajectory times t−τ (WaveField phase).
    // Crossfade across a kernel rebuild: F(t) = α·F_old(t) + (1−α)·F_new(t),
    // α = blendStepsLeft_/blendStepsTotal_ stepping from 1 → 0 over N fast steps.
    // Each group convolves against its OWN components-derived η (etaAtGroup),
    // so old and new groups are self-consistent in their past lookups; we
    // never touch the shared η history.
    auto convolveOnto = [&](const std::vector<Group>& gs, Eigen::VectorXd& Fout)
    {
        for (const Group& g : gs)
        {
            if (g.K.rows() == 0) continue;
            const int nLag = static_cast<int>(g.K.rows());
            const double dt = g.dt;
            for (int lag = 0; lag < nLag; ++lag)
            {
                const double tWave = t - static_cast<double>(lag) * dt;
                const double etaLag = etaAtGroup(g, tWave);
                Fout.noalias() += g.K.row(lag).transpose() * (etaLag * dt);
            }
        }
    };

    if (blendStepsLeft_ > 0 && !oldGroups_.empty() && blendStepsTotal_ > 0)
    {
        Eigen::VectorXd Fold = Eigen::VectorXd::Zero(D);
        Eigen::VectorXd Fnew = Eigen::VectorXd::Zero(D);
        convolveOnto(oldGroups_, Fold);
        convolveOnto(groups_,    Fnew);

        const double alpha =
            static_cast<double>(blendStepsLeft_) / static_cast<double>(blendStepsTotal_);
        F = alpha * Fold + (1.0 - alpha) * Fnew;

        --blendStepsLeft_;
        if (blendStepsLeft_ <= 0)
        {
            oldGroups_.clear();
            blendStepsLeft_  = 0;
            blendStepsTotal_ = 0;
        }
    }
    else
    {
        convolveOnto(groups_, F);
    }

    if (force6) coupled::scatterToSixDof(F, skCfg_, *force6);

    // ===== 临时诊断（仅 fkImpulse 主模式 part=FkAndDf 时启用）=====
    // 把每步 |F|, F0..F{D-1} 写到 <casePath>/coupled_turning/fkImpulse_force_dump.csv
    // 列格式刻意与 DirectFKPlusDfWaveForceProvider 的 dump 对齐前 3 列，
    // 这样可以直接拿两份 csv 在 python 里按 t 对比。
    // part = DfOnly 时不写（避免和 directFKdf 那条线重复 dump）。
    if (part_ == FKImpulsePart::FkAndDf)
    {
        static std::ofstream dump([this]() {
            const std::string dir = casePath_ + "coupled_turning/";
            std::filesystem::create_directories(dir);
            return dir + "fkImpulse_force_dump.csv";
        }());
        static bool wroteHeader = false;
        if (!wroteHeader)
        {
            dump << "t,|F|";
            for (int j = 0; j < F.size(); ++j) dump << ",F" << j;
            dump << "\n";
            wroteHeader = true;
        }
        dump << t << "," << F.norm();
        for (int j = 0; j < F.size(); ++j) dump << "," << F(j);
        dump << "\n";
    }
    // ===== 临时诊断结束 =====

    return F;
}
