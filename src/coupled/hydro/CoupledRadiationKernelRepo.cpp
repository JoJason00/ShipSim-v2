#include "CoupledRadiationKernelRepo.h"
#include <filesystem>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cmath>

CoupledRadiationKernelRepo::CoupledRadiationKernelRepo(
    const ShipConfig& ship,
    const SeakeepingConfig& skCfg,
    const std::string& casePath)
    : ship_(ship), skCfg_(skCfg), casePath_(casePath)
{
}

std::string CoupledRadiationKernelRepo::modesTag(const SeakeepingConfig& skCfg)
{
    std::ostringstream ss;
    for (std::size_t i = 0; i < skCfg.modes.size(); ++i)
    {
        if (i) ss << "_";
        ss << skCfg.modes[i];
    }
    return ss.str();
}

std::string CoupledRadiationKernelRepo::keyDouble(double x)
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

void CoupledRadiationKernelRepo::scan() const
{
    if (scanned_) return;
    scanned_ = true;

    const auto dir = std::filesystem::path(casePath_) / "kernel_cache";
    if (!std::filesystem::exists(dir))
        throw std::runtime_error("CoupledRadiationKernelRepo: kernel_cache directory not found.");

    for (const auto& de : std::filesystem::directory_iterator(dir))
    {
        if (!de.is_regular_file()) continue;
        if (de.path().extension() != ".bin") continue;

        RadiationKernelData k;
        if (!RadiationKernelCache::load(de.path().string(), k))
            continue;

        if (k.DOF != skCfg_.DOF) continue;
        if (k.modes != skCfg_.modes) continue;
        // NE (panel count) may differ from Seakeeping.Panel.NE: use radiation kernels from
        // another mesh resolution when inventories do not match the current case NE.
        if (k.NE != skCfg_.Panel.NE)
        {
            std::cout << "[RadiationKernel] NE mismatch (using library anyway): file NE=" << k.NE
                      << "  case Panel.NE=" << skCfg_.Panel.NE << "  file=" << de.path().filename().string()
                      << "\n";
        }

        std::cout << "[RadiationKernel] read library " << de.path().filename().string()
                  << "  Fn=" << k.Fn
                  << "  dt_offline=" << k.dt
                  << "  TG=" << k.TG
                  << "  Klag_layers=" << k.Klag.size()
                  << "  DOF=" << k.DOF << "\n";

        entries_.push_back({ de.path().string(), std::move(k) });
    }

    if (entries_.empty())
        throw std::runtime_error("CoupledRadiationKernelRepo: no matching radiation kernels found.");

    std::sort(entries_.begin(), entries_.end(), [](const FileEntry& a, const FileEntry& b)
    {
        return a.kernel.Fn < b.kernel.Fn;
    });
}

RadiationKernelData CoupledRadiationKernelRepo::resampleKernel(const RadiationKernelData& src, double dtFast)
{
    if (std::abs(src.dt - dtFast) < 1e-12)
        return src;

    RadiationKernelData out = src;
    const double tMem = src.dt * src.TG;
    const int TGnew = std::max(1, static_cast<int>(std::ceil(tMem / dtFast)));
    out.dt = dtFast;
    out.TG = TGnew;
    out.Klag.assign(static_cast<std::size_t>(TGnew + 1), Eigen::MatrixXd::Zero(src.DOF, src.DOF));
    out.K0 = Eigen::MatrixXd::Zero(src.DOF, src.DOF);

    auto sampleAt = [&](double t) -> Eigen::MatrixXd
    {
        if (t <= 0.0) return src.Klag.front();
        if (t >= src.dt * src.TG) return src.Klag.back();

        const double x = t / src.dt;
        const int i0 = static_cast<int>(std::floor(x));
        const int i1 = std::min(src.TG, i0 + 1);
        const double a = x - i0;
        return (1.0 - a) * src.Klag[static_cast<std::size_t>(i0)]
             + a * src.Klag[static_cast<std::size_t>(i1)];
    };

    for (int i = 0; i <= TGnew; ++i)
        out.Klag[static_cast<std::size_t>(i)] = sampleAt(i * dtFast);

    out.K0 = out.Klag.front();
    return out;
}

double CoupledRadiationKernelRepo::fnSpacingMedian(bool& valid) const
{
    scan();

    if (fnSpacingCacheValid_)
    {
        valid = fnSpacingValid_;
        return fnSpacingMedian_;
    }

    // entries_ is already sorted ascending by Fn (see scan()).
    std::vector<double> gaps;
    gaps.reserve(entries_.size());
    double prev = entries_.empty() ? 0.0 : entries_.front().kernel.Fn;
    for (std::size_t i = 1; i < entries_.size(); ++i)
    {
        const double g = entries_[i].kernel.Fn - prev;
        if (g > 1e-6)
        {
            gaps.push_back(g);
            prev = entries_[i].kernel.Fn;
        }
    }

    fnSpacingCacheValid_ = true;
    if (gaps.empty())
    {
        fnSpacingValid_ = false;
        fnSpacingMedian_ = 0.0;
    }
    else
    {
        std::sort(gaps.begin(), gaps.end());
        const std::size_t m = gaps.size() / 2;
        fnSpacingMedian_ = (gaps.size() % 2 == 1)
                               ? gaps[m]
                               : 0.5 * (gaps[m - 1] + gaps[m]);
        fnSpacingValid_ = true;
    }

    valid = fnSpacingValid_;
    return fnSpacingMedian_;
}

RadiationKernelData CoupledRadiationKernelRepo::getKernel(double Fn, double dtFast) const
{
    scan();

    if (entries_.size() == 1)
        return resampleKernel(entries_.front().kernel, dtFast);

    auto it = std::lower_bound(entries_.begin(), entries_.end(), Fn,
        [](const FileEntry& e, double x) { return e.kernel.Fn < x; });

    if (it == entries_.begin())
        return resampleKernel(it->kernel, dtFast);
    if (it == entries_.end())
        return resampleKernel(entries_.back().kernel, dtFast);

    const auto& hi = *it;
    const auto& lo = *(it - 1);
    if (std::abs(hi.kernel.Fn - lo.kernel.Fn) < 1e-12)
        return resampleKernel(lo.kernel, dtFast);

    RadiationKernelData klo = resampleKernel(lo.kernel, dtFast);
    RadiationKernelData khi = resampleKernel(hi.kernel, dtFast);

    const double a = (Fn - lo.kernel.Fn) / (hi.kernel.Fn - lo.kernel.Fn);

    RadiationKernelData out = klo;
    out.Fn = Fn;
    out.U = (1.0 - a) * klo.U + a * khi.U;
    out.A_inf = (1.0 - a) * klo.A_inf + a * khi.A_inf;
    out.B = (1.0 - a) * klo.B + a * khi.B;
    out.C_prime = (1.0 - a) * klo.C_prime + a * khi.C_prime;
    out.K0 = (1.0 - a) * klo.K0 + a * khi.K0;
    out.Klag.resize(klo.Klag.size());
    for (std::size_t i = 0; i < out.Klag.size(); ++i)
        out.Klag[i] = (1.0 - a) * klo.Klag[i] + a * khi.Klag[i];
    return out;
}
