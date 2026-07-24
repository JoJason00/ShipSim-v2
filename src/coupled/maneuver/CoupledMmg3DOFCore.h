#pragma once

#include "ILowFreqManeuverSolver.h"
#include "../../config/CaseConfig.h"
#include "../../config/MmgConfig.h"
#include <vector>
#include <memory>

class CoupledMmg3DOFCore : public ICoupledLowFreqManeuverSolver
{
public:
    CoupledMmg3DOFCore(const ShipConfig& ship, const MmgConfig& mmg);

    CoupledSlowState3DOF step(
        const CoupledSlowState3DOF& s,
        const CoupledExternalLoads3DOF& ext,
        double deltaCmd,
        double dt) override;

    double referenceSpeedFromFn(double Fn) const override;

private:
    struct State
    {
        double u = 0.0;
        double v = 0.0;
        double r = 0.0;
        double xe = 0.0;
        double ye = 0.0;
        double psi = 0.0;
        double del = 0.0;
    };

    struct Deriv
    {
        double du = 0.0;
        double dv = 0.0;
        double dr = 0.0;
        double dxe = 0.0;
        double dye = 0.0;
        double dpsi = 0.0;
        double ddel = 0.0;
    };

    struct ForceBreakdown
    {
        double U = 0.0;
        double beta = 0.0;
        double r_nd = 0.0;
        double XH = 0.0;
        double YH = 0.0;
        double NH = 0.0;
        double Ru = 0.0;
        double wp = 0.0;
        double Jp = 0.0;
        double KT = 0.0;
        double XP = 0.0;
        double uR = 0.0;
        double vR = 0.0;
        double UR = 0.0;
        double aR = 0.0;
        double FN = 0.0;
        double XR = 0.0;
        double YR = 0.0;
        double NR = 0.0;
    };

private:
    void compute_hullCoeff() {}
    void compute_propellerCoeff() {}
    void compute_rudderCoeff() {}
    void compute_addedCoeff() {}
    void setupModelScaleS175();

    /// When `mmg.coupledMmg3DOF.enabled`, load from `case.json` (see `CoupledMmg3DOFSettings`).
    void initFromCoupledMmg3DOFSettings(const MmgConfig& mmg);
    /// Archived full-scale S175 defaults (175 m equivalent); kept in the .cpp for side-by-side checks.
    void initLegacyHardcodedFullScaleS175();

    ForceBreakdown evalForces(const State& s) const;
    Deriv rhs(const State& s, const CoupledExternalLoads3DOF& ext,
        double deltaCmdRad, double localDt) const;

    /// Integrate the slow-heel roll ODE over `dt` with internal RK4 sub-steps
    /// (so dtSlow need not resolve the ~1.6 s roll period). Forcing is frozen at
    /// the slow-step start: M_ext = zG*(my*vdot + mx*u*r) - zH*YH - zG*YR.
    /// `Kext` is an external roll moment added to M_ext — currently the mean
    /// second-order (wave-drift) heel moment K2 (≈0 in head/following seas).
    /// No-op unless rollEnabled_. Updates phi, p in place.
    void integrateRoll(double& phi, double& p,
        double YH, double YR, double vdot, double u, double r,
        double Kext, double dt) const;
    State rk4Step(const State& s, const CoupledExternalLoads3DOF& ext,
        double dt, double deltaCmdRad) const;

    static double deg2rad(double deg);
    static double clamp(double x, double lo, double hi);
    static State toInner(const CoupledSlowState3DOF& s);
    static CoupledSlowState3DOF toOuter(const State& s, double t, double Fn);

private:
    ShipConfig ShipCfg;
    MmgConfig mmg_;

    // All scalars below are filled in the ctor by either initFromCoupledMmg3DOFSettings()
    // or initLegacyHardcodedFullScaleS175() ?? see CoupledMmg3DOFCore.cpp.

    double rho_{};
    double g_{};
    double pi_{};

    double Lpp_{};
    double Lwl_{};
    double B_{};
    double d_{};
    double Cb_{};

    double m_{};
    double xG_{};  // read from CoupledMmg3DOF.xG_m (unused when body origin is at CG)
    double kzz_{};
    double Izz_{};

    double Dp_{};
    double HR_{};
    double BR_{};
    double AR_{};

    double lambdaR_{};
    double eta_{};

    double xP_nd_{};
    double xR_nd_{};
    double xH_nd_{};
    double lR_nd_{};

    double Xbb_{};
    double Xbr_{};
    double Xrr_{};

    double Yb_{};
    double Yr_{};
    double Ybbb_{};
    double Ybbr_{};
    double Ybrr_{};
    double Yrrr_{};

    double Nb_{};
    double Nr_{};
    double Nbbb_{};
    double Nbbr_{};
    double Nbrr_{};
    double Nrrr_{};

    double mx_nd_{};
    double my_nd_{};
    double Jz_nd_{};
    double Xuu_{};

    double mx_{};
    double my_{};
    double Jz_{};

    double epsilon_{};
    double kappa_{};
    double tR_{};
    double aH_{};
    double gammaR_plus_{};
    double gammaR_minus_{};

    double deltaRateDeg_{};

    double C1_{};
    double wp0_{};
    double tp_{};
    double J0_{};
    double J1_{};
    double J2_{};

    double Fn_{};
    double U0_{};
    double np_rps_{};

    double dt_{};
    double tMax_{};
    double approachTime_{};
    double stopPsiDeg_{};
    double maxDeltaDeg_{};

    // ---- Optional 4-DOF roll (slow heel) ----
    bool   rollEnabled_{ false };
    double GM_{};        // metacentric height [m] (from Ship.Mass.GM)
    double I44_{};       // Ixx + Jxx [kg*m^2]
    double C44roll_{};   // m*g*GM [N*m/rad]
    double rollBlin_{};  // linear roll damping
    double rollBquad_{}; // quadratic roll damping
    double zH_{};        // = zH_over_d * d  [m]
    double zR_{};        // = zR_over_d * d  [m]
    double zG_{};        // = zG_over_d * d  [m]
};
