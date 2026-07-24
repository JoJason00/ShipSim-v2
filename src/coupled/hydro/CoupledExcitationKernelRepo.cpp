#include "CoupledExcitationKernelRepo.h"
#include "../../const/Const.h"
#include <filesystem>
#include <algorithm>
#include <stdexcept>
#include <cmath>

CoupledExcitationKernelRepo::CoupledExcitationKernelRepo(
    const ShipConfig& ship,
    const SeakeepingConfig& skCfg,
    const std::string& casePath)
    : ship_(ship), skCfg_(skCfg), casePath_(casePath)
{
}

double CoupledExcitationKernelRepo::wrapAngle(double x)
{
    while (x <= -PI) x += 2.0 * PI;
    while (x > PI) x -= 2.0 * PI;
    return x;
}

void CoupledExcitationKernelRepo::scan() const
{
    if (scanned_) return;
    scanned_ = true;

    const auto dir = std::filesystem::path(casePath_) / "kernel_cache" / "excitation";
    if (!std::filesystem::exists(dir))
        throw std::runtime_error("CoupledExcitationKernelRepo: kernel_cache/excitation directory not found.");

    for (const auto& de : std::filesystem::directory_iterator(dir))
    {
        if (!de.is_regular_file()) continue;
        if (de.path().extension() != ".bin") continue;

        TDGFExcitationKernelData k;
        if (!TDGFExcitationKernelIO::load(de.path().string(), k))
            continue;
        if (k.DOF != skCfg_.DOF) continue;
        if (k.modes != skCfg_.modes) continue;

        entries_.push_back({ de.path().string(), std::move(k) });
    }

    if (entries_.empty())
        throw std::runtime_error("CoupledExcitationKernelRepo: no excitation kernels found.");
}

TDGFExcitationKernelData CoupledExcitationKernelRepo::resampleKernel(const TDGFExcitationKernelData& src, double dtFast)
{
    if (std::abs(src.dt - dtFast) < 1e-12)
        return src;

    TDGFExcitationKernelData out = src;
    const int nOld = static_cast<int>(src.Qlag.rows());
    if (nOld <= 1)
    {
        out.dt = dtFast;
        return out;
    }

    const double tMax = src.dt * (nOld - 1);
    const int nNew = std::max(2, static_cast<int>(std::ceil(tMax / dtFast)) + 1);
    out.dt = dtFast;
    out.Qlag.resize(nNew, src.Qlag.cols());

    for (int i = 0; i < nNew; ++i)
    {
        const double t = i * dtFast;
        if (t >= tMax)
        {
            out.Qlag.row(i) = src.Qlag.row(nOld - 1);
            continue;
        }

        const double x = t / src.dt;
        const int i0 = static_cast<int>(std::floor(x));
        const int i1 = std::min(nOld - 1, i0 + 1);
        const double a = x - i0;
        out.Qlag.row(i) = (1.0 - a) * src.Qlag.row(i0) + a * src.Qlag.row(i1);
    }

    return out;
}

namespace
{
    // Median of the gaps between sorted, de-duplicated nodes. Returns false in
    // `valid` when there are fewer than two distinct nodes.
    double medianNodeSpacing(std::vector<double> vals, double dedupeTol, bool& valid)
    {
        valid = false;
        std::sort(vals.begin(), vals.end());

        std::vector<double> gaps;
        gaps.reserve(vals.size());
        double prev = vals.empty() ? 0.0 : vals.front();
        bool havePrev = !vals.empty();
        for (std::size_t i = 1; i < vals.size(); ++i)
        {
            const double g = vals[i] - prev;
            if (g > dedupeTol)
            {
                gaps.push_back(g);
                prev = vals[i];
            }
        }
        (void)havePrev;

        if (gaps.empty())
            return 0.0;

        std::sort(gaps.begin(), gaps.end());
        valid = true;
        const std::size_t m = gaps.size() / 2;
        return (gaps.size() % 2 == 1)
                   ? gaps[m]
                   : 0.5 * (gaps[m - 1] + gaps[m]);
    }
}

CoupledExcitationKernelRepo::GridSpacing
CoupledExcitationKernelRepo::gridSpacing(double omegaIncident) const
{
    scan();

    if (gsCacheValid_ && std::abs(gsCacheOmega_ - omegaIncident) <= 1e-6)
        return gsCache_;

    std::vector<double> fns;
    std::vector<double> betas;
    for (const auto& e : entries_)
    {
        if (std::abs(e.kernel.omegaIncident - omegaIncident) > 1e-6)
            continue;
        fns.push_back(e.kernel.Fn);
        betas.push_back(wrapAngle(e.kernel.betaRel));
    }

    GridSpacing gs;
    gs.dFn = medianNodeSpacing(fns, 1e-6, gs.fnValid);
    gs.dBetaRad = medianNodeSpacing(betas, 1e-4, gs.betaValid);

    gsCacheValid_ = true;
    gsCacheOmega_ = omegaIncident;
    gsCache_ = gs;
    return gs;
}

TDGFExcitationKernelData CoupledExcitationKernelRepo::getKernel(
    double Fn,
    double betaRel,
    double omegaIncident,
    double dtFast) const
{
    scan();

    std::vector<const FileEntry*> candidates;
    for (const auto& e : entries_)
    {
        if (std::abs(e.kernel.omegaIncident - omegaIncident) > 1e-6)
            continue;
        candidates.push_back(&e);
    }
    if (candidates.empty())
        throw std::runtime_error("CoupledExcitationKernelRepo: no candidate kernels at requested omegaIncident.");

    auto metric = [&](const FileEntry* e) -> double
    {
        const double dFn = std::abs(e->kernel.Fn - Fn);
        const double dB = std::abs(wrapAngle(e->kernel.betaRel - betaRel));
        return dFn + 0.1 * dB;
    };

    std::sort(candidates.begin(), candidates.end(), [&](const FileEntry* a, const FileEntry* b)
    {
        return metric(a) < metric(b);
    });

    if (candidates.size() == 1)
        return resampleKernel(candidates.front()->kernel, dtFast);

    const auto& k0 = candidates[0]->kernel;
    const auto& k1 = candidates[1]->kernel;

    const double b0 = wrapAngle(k0.betaRel);
    const double b1 = wrapAngle(k1.betaRel);
    const double denom = std::abs(wrapAngle(b1 - b0));
    if (denom < 1e-8)
        return resampleKernel(k0, dtFast);

    const double alpha = std::clamp(std::abs(wrapAngle(betaRel - b0)) / denom, 0.0, 1.0);
    TDGFExcitationKernelData a0 = resampleKernel(k0, dtFast);
    TDGFExcitationKernelData a1 = resampleKernel(k1, dtFast);

    if (a0.Qlag.rows() != a1.Qlag.rows() || a0.Qlag.cols() != a1.Qlag.cols())
        return a0;

    TDGFExcitationKernelData out = a0;
    out.Fn = Fn;
    out.betaRel = betaRel;
    out.U = (1.0 - alpha) * a0.U + alpha * a1.U;
    out.omegaEncounter = (1.0 - alpha) * a0.omegaEncounter + alpha * a1.omegaEncounter;
    out.Qlag = (1.0 - alpha) * a0.Qlag + alpha * a1.Qlag;
    return out;
}
