#pragma once

#include "KernelKey.h"
#include "EncounterClassifier.h"
#include "TDGFExcitationKernelIO.h"
#include "../wave/WaveField.h"

#include <Eigen/Dense>
#include <vector>
#include <deque>

// Source of excitation kernels keyed by KernelKey. Implement this over the
// existing kernel library (e.g. wrap CoupledExcitationKernelRepo, or the
// pure-seakeeping FK-impulse store). Returning nullptr drops that group.
class IExcitationKernelSource
{
public:
    virtual ~IExcitationKernelSource() = default;
    // dt is the integration step the caller will convolve at; the source is
    // expected to return a kernel resampled to dt.
    virtual const TDGFExcitationKernelData* fetch(const KernelKey& key,
                                                  double dt) = 0;
};

// Groups a WaveField's components by the kernel they need, keeps one summed
// incident-elevation ring history + one kernel per group, and at every step
// returns Σ_groups ( Qlag_g ⊛ η_g ). This is the linear-convolution identity
// K*(η₁+η₂+…)=K*η₁+K*η₂+… made concrete: same-kernel components are pre-summed
// into one convolution; different kernels stay separate and the forces add.
//
// Straight-course (β,U constant) -> partition once via prepare(); maneuvering
// callers re-prepare() when β/U cross a bucket (reuse the adaptive-refresh
// policy already in the coupled side).
class WaveForceAssembler
{
public:
    void prepare(const WaveField& field,
                 double U,
                 double psi,
                 const EncounterClassifierConfig& cfg,
                 IExcitationKernelSource& src,
                 double dt);

    // Excitation force on the active DOF modes at (x, y, t). Size = kernel DOF
    // (0 if not prepared / no kernels).
    Eigen::VectorXd evaluate(double x, double y, double t);

    int dof() const { return dof_; }
    std::size_t groupCount() const { return groups_.size(); }

private:
    struct Group
    {
        std::vector<WaveComponent> comps;     // components sharing one kernel
        const TDGFExcitationKernelData* k = nullptr;
        std::deque<double> etaHist;           // newest at back, length = Qlag rows
    };

    std::vector<Group> groups_;
    int dof_ = 0;
};
