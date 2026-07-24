#pragma once

#include "ILowFreqManeuverSolver.h"
#include "../../config/CaseConfig.h"
#include "../../config/MmgConfig.h"

/// 3-DOF slow manoeuvring model using the SJTU thesis (Son–Nomoto) hull/prop/rudder formulation.
/// Roll is omitted; hull forces use \(v'=v/V\), \(r'=rL/V\) with \(V=\sqrt{u^2+v^2}\).
class CoupledSjtuMmg3DOFCore : public ICoupledLowFreqManeuverSolver
{
public:
    CoupledSjtuMmg3DOFCore(const ShipConfig& ship, const MmgConfig& mmg);

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
        double V = 0.0;
        double vp = 0.0;
        double rp = 0.0;
        double beta = 0.0;
        double X0 = 0.0;
        double Y0 = 0.0;
        double N0 = 0.0;
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
    void setupAddedMassAndInertia();
    double ittcCf(double speedU) const;
    double resistanceR(double speedU) const;

    ForceBreakdown evalForces(const State& s) const;
    Deriv rhs(const State& s, const CoupledExternalLoads3DOF& ext,
        double deltaCmdRad, double localDt) const;
    State rk4Step(const State& s, const CoupledExternalLoads3DOF& ext,
        double dt, double deltaCmdRad) const;

    static double deg2rad(double deg);
    static double clamp(double x, double lo, double hi);
    static State toInner(const CoupledSlowState3DOF& s);
    static CoupledSlowState3DOF toOuter(const State& s, double t, double Fn);

private:
    ShipConfig ShipCfg;
    CoupledSjtuMmg3DOFSettings s_;

    double pi_{};
    double Izz_{};
    double mx_{};
    double my_{};
    double Jz_{};

    double lambdaR_{};
    double eta_{};

    double U0_{};
};
