#include "CumminsTimeStepper.h"

#include <algorithm>
#include <stdexcept>

namespace cummins
{
    void CumminsTimeStepper::init(StepperConfig cfg)
    {
        cfg_ = std::move(cfg);

        const int D = static_cast<int>(cfg_.Mtotal.rows());

        if (cfg_.Mtotal.cols() != D
            || cfg_.Bconst.rows() != D || cfg_.Bconst.cols() != D
            || cfg_.Ctotal.rows() != D || cfg_.Ctotal.cols() != D)
        {
            throw std::runtime_error("CumminsTimeStepper: system matrices must be DxD square.");
        }
        if (cfg_.dt <= 0.0)
            throw std::runtime_error("CumminsTimeStepper: dt must be positive.");
        if (cfg_.Konline.empty())
            throw std::runtime_error("CumminsTimeStepper: Konline must contain at least K[0].");

        // Newmark implicit operator: A_coef = M + γ·dt·B + β·dt²·C
        Acoef_ = cfg_.Mtotal
            + gamma_ * cfg_.dt * cfg_.Bconst
            + beta_ * cfg_.dt * cfg_.dt * cfg_.Ctotal;

        solver_.compute(Acoef_);
        if (!solver_.isInvertible())
            throw std::runtime_error("CumminsTimeStepper: Newmark operator A_coef is singular.");

        state_.q = Eigen::VectorXd::Zero(D);
        state_.v = Eigen::VectorXd::Zero(D);
        state_.a = Eigen::VectorXd::Zero(D);
        state_.vHist.clear();
        n_ = 0;
        firstStep_ = true;
    }

    void CumminsTimeStepper::setInitialState(const StepperState& s0)
    {
        const int D = dof();

        auto sizeOrZero = [D](const Eigen::VectorXd& x) -> Eigen::VectorXd
        {
            if (x.size() == D) return x;
            return Eigen::VectorXd::Zero(D);
        };

        state_.q = sizeOrZero(s0.q);
        state_.v = sizeOrZero(s0.v);
        state_.a = sizeOrZero(s0.a);
        state_.vHist = s0.vHist;

        for (auto& vv : state_.vHist)
            if (vv.size() != D)
                vv = Eigen::VectorXd::Zero(D);

        trimHistoryToWindow();

        // If caller supplied a non-empty history, treat this as a continuation —
        // the supplied acceleration is honoured and Newmark proceeds from there.
        firstStep_ = state_.vHist.empty();
        n_ = static_cast<int>(state_.vHist.size());
    }

    void CumminsTimeStepper::step(const Eigen::VectorXd& F_exc)
    {
        const int D = dof();
        const double dt = cfg_.dt;
        const int M = maxMemory();

        if (F_exc.size() != D)
            throw std::runtime_error("CumminsTimeStepper::step: F_exc has wrong size.");

        // First step on a fresh start: solve M·a = F − B·v − C·q for a consistent a₀.
        if (firstStep_)
        {
            Eigen::VectorXd rhs = F_exc - cfg_.Bconst * state_.v - cfg_.Ctotal * state_.q;
            if (cfg_.rollIndex >= 0)
                rhs(cfg_.rollIndex) += cfg_.rollVisc.moment(state_.v(cfg_.rollIndex));

            state_.a = cfg_.Mtotal.fullPivLu().solve(rhs);

            state_.vHist.push_back(state_.v);
            trimHistoryToWindow();
            ++n_;
            firstStep_ = false;
            return;
        }

        // Newmark predictor (no acceleration yet).
        const Eigen::VectorXd qPred = state_.q + dt * state_.v
            + dt * dt * (0.5 - beta_) * state_.a;
        const Eigen::VectorXd vPred = state_.v + dt * (1.0 - gamma_) * state_.a;

        // Trapezoid convolution memory: ∫ K(t-τ) v(τ) dτ
        // K0 already folded into Bconst as 0.5·dt·K0; loop covers lag 1..M-1 full weight + lag M half weight.
        Eigen::VectorXd mem = Eigen::VectorXd::Zero(D);
        const int histLen = static_cast<int>(state_.vHist.size());
        const int maxLag = std::min(histLen, M);

        for (int lag = 1; lag < maxLag; ++lag)
        {
            mem.noalias() += dt * cfg_.Konline[static_cast<std::size_t>(lag)]
                * state_.vHist[static_cast<std::size_t>(histLen - lag)];
        }
        if (maxLag >= 1)
        {
            mem.noalias() += 0.5 * dt * cfg_.Konline[static_cast<std::size_t>(maxLag)]
                * state_.vHist[static_cast<std::size_t>(histLen - maxLag)];
        }

        Eigen::VectorXd rhs = F_exc
            - cfg_.Bconst * vPred
            - cfg_.Ctotal * qPred
            - mem;

        if (cfg_.rollIndex >= 0)
            rhs(cfg_.rollIndex) += cfg_.rollVisc.moment(vPred(cfg_.rollIndex));

        const Eigen::VectorXd aNew = solver_.solve(rhs);
        state_.q = qPred + beta_ * dt * dt * aNew;
        state_.v = vPred + gamma_ * dt * aNew;
        state_.a = aNew;

        state_.vHist.push_back(state_.v);
        trimHistoryToWindow();
        ++n_;
    }

    void CumminsTimeStepper::trimHistoryToWindow()
    {
        const int M = maxMemory();
        // Keep the last (M+1) entries — the convolution touches lag ∈ [1, M].
        const int keep = std::max(1, M + 1);
        const int sz = static_cast<int>(state_.vHist.size());
        if (sz > keep)
        {
            state_.vHist.erase(state_.vHist.begin(),
                state_.vHist.begin() + (sz - keep));
        }
    }
}
