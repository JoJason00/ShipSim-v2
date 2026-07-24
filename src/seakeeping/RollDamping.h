#pragma once

#include <string>
#include <vector>
#include <cmath>

struct RollDecaySample
{
    double t = 0.0;   // [s]
    double phi = 0.0; // [rad]
};

struct RollDecayPeak
{
    double t = 0.0;       // [s]
    double phi_abs = 0.0; // absolute peak amplitude [rad]
};

struct DecayPolyCoeff
{
    // Delta phi = a*phi_m + b*phi_m^2 + c*phi_m^3
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;

    double omega_n = 0.0; // [rad/s]
    double T_n = 0.0;     // [s]
};

struct RollViscDamping
{
    double B44_lin = 0.0; // [moment * s]
    double B44_quad = 0.0; // [moment * s^2]
    double B44_cube = 0.0; // [moment * s^3]

    double moment(double phidot) const
    {
        return -(B44_lin * phidot
            + B44_quad * std::abs(phidot) * phidot
            + B44_cube * phidot * phidot * phidot);
    }
};

enum class RollDampingMode
{
    DirectCoefficients,
    FromDecayCsv
};

struct RollDampingDirectConfig
{
    double B44_lin = 0.0;
    double B44_quad = 0.0;
    double B44_cube = 0.0;
};

struct RollDampingDecayCsvConfig
{
    std::string csvPath;
    int polyOrder = 2;      // 2: a,b ; 3: a,b,c
    int minPeakGap = 3;     // min sample gap between extracted peaks
    double IeffOverride = -1.0; // >0 overrides the estimated Ieff
    bool angleInDeg = false;    // whether CSV phi column is in degrees

    // Option A: after the analytic peak-decrement fit, refine (B44_lin, B44_quad)
    // by forward-integrating the 1-DOF roll ODE and minimising the misfit to the
    // measured peak envelope. Removes the equivalent-linearisation bias so the
    // simulated decay tracks the measured one. false keeps the legacy behaviour.
    bool refine = false;
    int refineSkipFirstPeaks = 1; // drop this many leading peaks (release transient)
};

struct RollDampingConfig
{
    bool enabled = false;
    RollDampingMode mode = RollDampingMode::DirectCoefficients;
    RollDampingDirectConfig direct;
    RollDampingDecayCsvConfig decay;
};

struct RollDecayIdentificationResult
{
    std::vector<RollDecayPeak> peaks;
    DecayPolyCoeff poly;
    double Ieff = 0.0;
    RollViscDamping damping;
};

class RollDampingBuilder
{
public:
    static RollViscDamping build(
        const std::string& casePath,
        const RollDampingConfig& cfg,
        double mass,
        double GM);

    static std::vector<RollDecaySample> readDecayCsv(
        const std::string& path,
        bool angleInDeg);

    static std::vector<RollDecayPeak> extractAbsolutePeaks(
        const std::vector<RollDecaySample>& samples,
        int minIndexGap = 3);

    static double estimateNaturalOmega(
        const std::vector<RollDecayPeak>& peaks);

    static DecayPolyCoeff fitPeakDecrement(
        const std::vector<RollDecayPeak>& peaks,
        int order = 2);

    static double estimateIeffFromMassGM(
        double mass,
        double GM,
        double omega_n);

    static RollViscDamping convertToPhysicalDamping(
        const DecayPolyCoeff& poly,
        double Ieff);

    static RollDecayIdentificationResult identifyFromDecay(
        const std::vector<RollDecaySample>& samples,
        double mass,
        double GM,
        int polyOrder = 2,
        int minPeakGap = 3,
        double IeffOverride = -1.0,
        bool refine = false,
        int refineSkipFirstPeaks = 1);

    // Forward-integrate the 1-DOF roll ODE (RK4, identical scheme to
    // runFreeRollDecay) for the given damping and return the absolute-peak
    // envelope (t, |phi|) sampled at internal dt.
    static std::vector<RollDecayPeak> simulateDecayPeaks(
        const RollViscDamping& rv,
        double Ieff,
        double C44,
        double phi0,
        double phidot0,
        double tEnd);

    // Relative-L2 misfit between simulated and measured absolute-peak envelopes,
    // matched by peak index after skipping the first skipFirstPeaks of each.
    static double envelopeMisfit(
        const std::vector<RollDecayPeak>& meas,
        const std::vector<RollDecayPeak>& sim,
        int skipFirstPeaks);

    // Option A: refine (B44_lin, B44_quad) by minimising envelopeMisfit against
    // the measured samples, starting from `init`. B44_cube is kept as-is.
    static RollViscDamping refineByEnvelope(
        const std::vector<RollDecaySample>& samples,
        const RollViscDamping& init,
        double Ieff,
        double C44,
        int skipFirstPeaks);
};