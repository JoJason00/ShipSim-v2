#include "WaveForceAssembler.h"
#include "../const/Const.h"

#include <cmath>
#include <unordered_map>

void WaveForceAssembler::prepare(const WaveField& field,
                                 double U,
                                 double psi,
                                 const EncounterClassifierConfig& cfg,
                                 IExcitationKernelSource& src,
                                 double dt)
{
    groups_.clear();
    dof_ = 0;

    std::unordered_map<KernelKey, std::size_t> index;
    for (const auto& w : field.components())
    {
        const KernelKey key = encounter_classifier::classify(w, U, psi, cfg);
        auto it = index.find(key);
        std::size_t gi;
        if (it == index.end())
        {
            gi = groups_.size();
            index.emplace(key, gi);
            Group g;
            g.k = src.fetch(key, dt);
            groups_.push_back(std::move(g));
        }
        else
        {
            gi = it->second;
        }
        groups_[gi].comps.push_back(w);
    }

    for (auto& g : groups_)
    {
        if (!g.k || g.k->Qlag.rows() == 0) continue;
        dof_ = std::max(dof_, static_cast<int>(g.k->Qlag.cols()));
        g.etaHist.assign(static_cast<std::size_t>(g.k->Qlag.rows()), 0.0);
    }
}

Eigen::VectorXd WaveForceAssembler::evaluate(double x, double y, double t)
{
    if (dof_ <= 0)
        return Eigen::VectorXd();

    Eigen::VectorXd F = Eigen::VectorXd::Zero(dof_);

    for (auto& g : groups_)
    {
        if (!g.k || g.k->Qlag.rows() == 0 || g.etaHist.empty())
            continue;

        // Pre-sum this group's incident elevation at the ship reference point
        // (same phase convention as kernel generation / WaveField).
        double eta = 0.0;
        for (const auto& w : g.comps)
        {
            const double k = w.omega * w.omega / G;
            const double phase =
                k * (x * std::cos(w.theta) + y * std::sin(w.theta))
                - w.omega * t + w.eps;
            eta += w.a * std::cos(phase);
        }

        g.etaHist.pop_front();
        g.etaHist.push_back(eta);

        // F += Σ_lag Qlag.row(lag)^T * η(t - lag·dt); newest η at back.
        const Eigen::MatrixXd& Q = g.k->Qlag;
        const int nLag = static_cast<int>(Q.rows());
        const int last = static_cast<int>(g.etaHist.size()) - 1;
        for (int lag = 0; lag < nLag; ++lag)
        {
            const double e = g.etaHist[static_cast<std::size_t>(last - lag)];
            if (e == 0.0) continue;
            F.noalias() += Q.row(lag).transpose() * e;
        }
    }

    return F;
}
