#pragma once

// IWaveForceProvider — abstract source of wave excitation force F_exc(t)
// used by the coupled fast-time-scale Cummins stepper.
//
// Implementations include:
//   ImpulseKernelWaveForceProvider   — Qlag excitation kernel × η
//   DirectPressureFKWaveForceProvider — incident pressure integration (FK)
//   FKImpulseWaveForceProvider       — fkImpulse library (fk±df)
//   DirectFKPlusDfWaveForceProvider  — directFK + df impulse ("directFKdf")
//
// The provider may carry per-window state (e.g. the η history needed by the impulse kernel).
// onWindowBegin() lets the provider decide whether to refresh cached kernels based on slow state.

#include "../common/CoupledTypes.h"
#include <Eigen/Dense>

class IWaveForceProvider
{
public:
    virtual ~IWaveForceProvider() = default;

    // Called once at the start of each fast window. Provider can refresh cached kernels
    // here. dtFast is required because providers may need to resample kernels.
    // win.slowEnd / win.t0 / win.t1 let impulse providers sample η_ref(t−τ) along
    // the predicted slow trajectory inside the window (FK impulse, etc.).
    virtual void onWindowBegin(const CoupledSlowState3DOF& slow0,
                               const CoupledEncounterState& enc0,
                               double dtFast,
                               const CoupledFastWindowInfo& win = {}) = 0;

    // Compute F_exc on active DOF modes for the given slow state / encounter at time t.
    // The returned vector has size skCfg.DOF and is ordered to match skCfg.modes.
    // If force6 is non-null, fill 6DOF force in fixed (surge,sway,heave,roll,pitch,yaw) order.
    virtual Eigen::VectorXd evaluate(const CoupledSlowState3DOF& slow,
                                     double t,
                                     const CoupledEncounterState& enc,
                                     Eigen::RowVectorXd* force6 = nullptr) = 0;

    // History sync hooks. Impulse-kernel mode needs η history continuation across windows;
    // direct-pressure mode is stateless and ignores these.
    virtual void setEtaHistory(const std::vector<double>& /*hist*/) {}
    virtual std::vector<double> etaHistory() const { return {}; }

    // Number of past η samples (lags) the provider needs to allocate. 0 = stateless.
    virtual int memoryLength() const { return 0; }
};
