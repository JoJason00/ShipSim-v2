#pragma once

// FKImpulseWaveForceProvider — coupled fast-loop wave force from the SAME
// `cases/{case}/fkImpulse/` FK-impulse library the pure-seakeeping path uses.
//
//     F_j(t) = Σ_groups Σ_{lag} K_j^g(lag·dt) · η_g(t - lag·dt) · dt
//
// η_g uses the fixed global sea (WaveField convention):
//   η = Σ_j a_j cos[-k_j(x cosθ_j + y sinθ_j) - ω_j t + ε_j]
// at the earth-frame CG (x_e, y_e) traced along the slow trajectory; ε_j never
// resets when K is refreshed. Past lags use (x_e, y_e) at t−τ via window
// interpolation / short kinematic extrapolation outside the window.
//
// The incident sea (enc.comps) is split into kernel groups:
//   - KernelKey = (β_rel bucket, region). Region from (U, β_rel, ω):
//       head/beam (cos β ≤ 0) → H (one ω-independent kernel);
//       following (cos β > 0) → F1/F2/F3 by ω band.
//   - Library holds 0–180° only: a port-side β_rel reuses its mirror
//     βm = 2π−β and flips the antisymmetric DOFs (sway/roll/yaw).
//   - Components sharing a key are pre-summed into one η_g, convolved once
//     (convolution linearity); different keys add.
// Regular wave = 1 component → 1 group → bit-identical to the legacy path.

#include "IWaveForceProvider.h"
#include "../../config/CoupledConfig.h"
#include "../../config/SeakeepingConfig.h"
#include "../../seakeeping/WaveForceRegion.h"
#include "../../wave/WaveComponent.h"

#include <Eigen/Dense>
#include <string>
#include <vector>

// Which columns of fkImpulse_*.csv enter the convolution kernel.
enum class FKImpulsePart
{
    FkAndDf,   // legacy fkImpulse mode: fk + df
    FkOnly,
    DfOnly     // diffraction impulse only (pairs with directFK)
};

class FKImpulseWaveForceProvider : public IWaveForceProvider
{
public:
    FKImpulseWaveForceProvider(const SeakeepingConfig& skCfg,
                               const std::string& casePath,
                               const std::string& shipName,
                               const HydroRefreshConfig& refreshCfg,
                               FKImpulsePart part = FKImpulsePart::FkAndDf);

    void onWindowBegin(const CoupledSlowState3DOF& slow0,
                       const CoupledEncounterState& enc0,
                       double dtFast,
                       const CoupledFastWindowInfo& win = {}) override;

    Eigen::VectorXd evaluate(const CoupledSlowState3DOF& slow,
                             double t,
                             const CoupledEncounterState& enc,
                             Eigen::RowVectorXd* force6 = nullptr) override;

    void setEtaHistory(const std::vector<double>& hist) override;
    std::vector<double> etaHistory() const override;
    int memoryLength() const override { return memoryLen_; }


private:
    struct FileKernel
    {
        double          Fn     = 0.0;
        double          dirRad = 0.0;
        WaveForceRegion region = WaveForceRegion::Head;
        double          dtSrc  = 0.0;
        int             tMot   = 0;
        Eigen::MatrixXd fk;        // (tMot+1) × 6 causal half
        Eigen::MatrixXd df;        // (tMot+1) × 6 (diffraction, unused yet)
    };

    struct Group
    {
        Eigen::MatrixXd K;                 // nLag × DOF (mode-projected, signed)
        double          dt = 0.0;
        std::vector<WaveComponent> comps;  // absolute θ, ε; a includes ramp
        std::vector<double>        etaHist;
    };

    void scanCaseDir();
    static bool loadFile(const std::string& path, FileKernel& out);
    const FileKernel* findNearest(double Fn, double betaRel,
                                  WaveForceRegion region) const;
    int  buildGroupKernel(const FileKernel& src, double dtFast,
                          bool mirrorPS, Eigen::MatrixXd& Kout) const;
    void rebuildGroups(const CoupledSlowState3DOF& slow0,
                       const CoupledEncounterState& enc0,
                       double dtFast);
    bool needGroupRefresh(const CoupledSlowState3DOF& slow0,
                          const CoupledEncounterState& enc0,
                          double dtFast) const;
    bool needsMidWindowKernelRefresh(const CoupledSlowState3DOF& slow,
                                     const CoupledEncounterState& enc) const;
    static std::size_t groupsSignature(const CoupledSlowState3DOF& slow0,
                                       const CoupledEncounterState& enc0);
    static double componentElevation(double x, double y, double t,
                                     const WaveComponent& w);
    double etaAtGroup(const Group& g, double tWave) const;
    double etaAllComponents(const std::vector<WaveComponent>& comps,
                            double x, double y, double t) const;
    CoupledSlowState3DOF slowAtTime(double tQuery) const;
    static CoupledSlowState3DOF interpSlow(const CoupledSlowState3DOF& a,
                                           const CoupledSlowState3DOF& b,
                                           double alpha);
    static CoupledSlowState3DOF extrapolateSlow(const CoupledSlowState3DOF& ref,
                                                double tRef, double tQuery);

private:
    SeakeepingConfig            skCfg_;
    std::string                 casePath_;
    std::string                 shipName_;
    HydroRefreshConfig          refreshCfg_;
    FKImpulsePart               part_ = FKImpulsePart::FkAndDf;

    std::vector<FileKernel>     kernels_;

    std::vector<Group>          groups_;
    std::vector<double>         etaHist_;
    int                         memoryLen_ = 0;
    double                      dtCached_  = 0.0;
    double                      FnCached_  = 0.0;
    double                      betaRebuildCached_ = 0.0;
    double                      psiRebuildCached_  = 0.0;
    std::size_t                 groupsSig_ = 0;
    bool                        groupsValid_ = false;

    // --- Smooth kernel crossfade on rebuild ------------------------------
    // Cause: each rebuild swaps K_β instantaneously while the past η it is
    // convolved with was acquired under the previous β. The resulting step
    // discontinuity in F has broadband content that pumps the near-resonant
    // roll/heave/pitch DOFs (cf. ξ4 valley never decaying during turning).
    // Fix: keep the previous Group set alongside the new one for a few fast
    // steps after every rebuild and linearly fade from F_old to F_new. The
    // η_hist itself is shared (we never reset it across rebuilds).
    std::vector<Group>          oldGroups_;
    int                         blendStepsLeft_  = 0;
    int                         blendStepsTotal_ = 0;

    CoupledSlowState3DOF        trajSlow0_{};
    CoupledSlowState3DOF        trajSlow1_{};
    double                      trajT0_ = 0.0;
    double                      trajT1_ = 0.0;
    bool                        trajValid_ = false;
};
