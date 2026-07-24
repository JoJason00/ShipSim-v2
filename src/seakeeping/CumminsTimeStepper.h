#pragma once

// CumminsTimeStepper — Newmark-β (β=0.25, γ=0.5) + trapezoidal radiation-memory convolution.
//
// Shared time-stepping core for:
//   - pure seakeeping (LinearCumminsTDGF::solveKernelCase)
//   - coupled seakeeping window (TDGFWindowSeakeepingSolver)
//
// The stepper is stateful. Drive it like this:
//   CumminsTimeStepper s;
//   s.init(cfg);                       // once per case
//   s.setInitialState(s0);             // optional; defaults to zero state
//   for (int n = 0; n < N; ++n) {
//       s.step(F_exc.row(n).transpose());
//       const auto& st = s.state();    // st.q, st.v, st.a, history
//   }
//
// Equation of motion (Cummins, Eq. 3.23):
//   (M + A_inf) q̈ + (B + b_phys + 0.5·dt·K0) q̇ + (C + C') q + ∫ K(t-τ) q̇(τ) dτ = F_exc(t)
//
// K0/2 is folded into Bconst so the trapezoid convolution can skip lag=0.

#include "RollDamping.h"

#include <Eigen/Dense>
#include <vector>

namespace cummins
{
    struct StepperConfig
    {
        Eigen::MatrixXd Mtotal;                 // physical mass + A_inf
        Eigen::MatrixXd Bconst;                 // physical damping + radiation B + 0.5·dt·K0
        Eigen::MatrixXd Ctotal;                 // physical stiffness + C'
        std::vector<Eigen::MatrixXd> Konline;   // K(t) sampled at dt steps, length M+1
        double dt = 0.0;
        int    rollIndex = -1;                  // DOF index of roll (-1 if not present)
        RollViscDamping rollVisc;
    };

    struct StepperState
    {
        Eigen::VectorXd q;                          // generalised displacement
        Eigen::VectorXd v;                          // generalised velocity
        Eigen::VectorXd a;                          // generalised acceleration
        std::vector<Eigen::VectorXd> vHist;         // full velocity history (vHist[n] = v at step n)
    };

    class CumminsTimeStepper
    {
    public:
        CumminsTimeStepper() = default;

        void init(StepperConfig cfg);

        // Set initial state. If s0.vHist is non-empty, the next step() will
        // skip the first-step initialiser (so the existing acceleration s0.a
        // is honoured — used by the coupled window solver for continuation).
        // If s0.vHist is empty, the next step() does the consistent
        // first-step solve (used by pure seakeeping starting from rest).
        void setInitialState(const StepperState& s0);

        // Advance one time step with external force F_exc at this step.
        void step(const Eigen::VectorXd& F_exc);

        const StepperState& state() const { return state_; }
        int    stepIndex()       const { return n_; }
        double dt()              const { return cfg_.dt; }
        int    dof()             const { return static_cast<int>(cfg_.Mtotal.rows()); }
        int    maxMemory()       const
        {
            return static_cast<int>(cfg_.Konline.size()) - 1;
        }

    private:
        void trimHistoryToWindow();

        StepperConfig                       cfg_;
        StepperState                        state_;
        Eigen::FullPivLU<Eigen::MatrixXd>   solver_;
        Eigen::MatrixXd                     Acoef_;
        double                              beta_ = 0.25;
        double                              gamma_ = 0.5;
        int                                 n_ = 0;
        bool                                firstStep_ = true;
    };
}
