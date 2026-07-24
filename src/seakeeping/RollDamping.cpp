#include "RollDamping.h"
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "const/Const.h"

namespace
{
    static bool tryParseTwoNumbers(const std::string& line, double& x, double& y)
    {
        std::string s = line;
        for (char& c : s)
        {
            if (c == ';' || c == '\t') c = ',';
        }

        std::stringstream ss(s);
        std::string a, b;
        if (!std::getline(ss, a, ',')) return false;
        if (!std::getline(ss, b, ',')) return false;

        try
        {
            x = std::stod(a);
            y = std::stod(b);
            return true;
        }
        catch (...)
        {
            return false;
        }
    }
}

std::vector<RollDecaySample> RollDampingBuilder::readDecayCsv(
    const std::string& path,
    bool angleInDeg)
{
    std::ifstream in(path);
    if (!in.is_open())
        throw std::runtime_error("Cannot open roll decay csv: " + path);

    std::vector<RollDecaySample> out;
    std::string line;

    while (std::getline(in, line))
    {
        if (line.empty()) continue;

        double t = 0.0, phi = 0.0;
        if (!tryParseTwoNumbers(line, t, phi))
            continue;

        if (angleInDeg)
            phi *= (PI / 180.0);

        out.push_back({ t, phi });
    }

    if (out.size() < 10)
        throw std::runtime_error("Too few decay samples in csv: " + path);

    return out;
}

std::vector<RollDecayPeak> RollDampingBuilder::extractAbsolutePeaks(
    const std::vector<RollDecaySample>& samples,
    int minIndexGap)
{
    if (samples.size() < 3)
        throw std::runtime_error("Roll decay samples are too short.");

    std::vector<RollDecayPeak> peaks;
    int lastAccepted = -1000000;

    for (int i = 1; i < static_cast<int>(samples.size()) - 1; ++i)
    {
        const double a0 = std::abs(samples[i - 1].phi);
        const double a1 = std::abs(samples[i].phi);
        const double a2 = std::abs(samples[i + 1].phi);

        const bool isPeak = (a1 >= a0 && a1 > a2);
        if (!isPeak) continue;
        if (i - lastAccepted < minIndexGap) continue;

        peaks.push_back({ samples[i].t, a1 });
        lastAccepted = i;
    }

    if (peaks.size() < 4)
        throw std::runtime_error("Not enough decay peaks extracted.");

    return peaks;
}

double RollDampingBuilder::estimateNaturalOmega(
    const std::vector<RollDecayPeak>& peaks)
{
    if (peaks.size() < 2)
        throw std::runtime_error("Not enough peaks to estimate omega_n.");

    double sumDt = 0.0;
    int cnt = 0;

    for (int i = 1; i < static_cast<int>(peaks.size()); ++i)
    {
        const double dt = peaks[i].t - peaks[i - 1].t;
        if (dt > 0.0)
        {
            sumDt += dt;
            ++cnt;
        }
    }

    if (cnt == 0)
        throw std::runtime_error("Invalid peak intervals.");

    const double meanHalfPeriod = sumDt / static_cast<double>(cnt);
    return PI / meanHalfPeriod;
}

DecayPolyCoeff RollDampingBuilder::fitPeakDecrement(
    const std::vector<RollDecayPeak>& peaks,
    int order)
{
    if (order != 2 && order != 3)
        throw std::runtime_error("polyOrder must be 2 or 3.");

    if (peaks.size() < static_cast<size_t>(order + 2))
        throw std::runtime_error("Not enough peaks for polynomial fitting.");

    const int n = static_cast<int>(peaks.size()) - 1;
    const int m = order;

    Eigen::MatrixXd A(n, m);
    Eigen::VectorXd y(n);

    for (int i = 0; i < n; ++i)
    {
        const double phiPrev = peaks[i].phi_abs;
        const double phiNext = peaks[i + 1].phi_abs;

        const double dphi = phiPrev - phiNext;
        const double phim = 0.5 * (phiPrev + phiNext);

        y(i) = dphi;
        A(i, 0) = phim;
        A(i, 1) = phim * phim;
        if (order == 3)
            A(i, 2) = phim * phim * phim;
    }

    const Eigen::VectorXd x = A.colPivHouseholderQr().solve(y);

    DecayPolyCoeff out;
    out.a = x(0);
    out.b = x(1);
    out.c = (order == 3) ? x(2) : 0.0;
    out.omega_n = estimateNaturalOmega(peaks);
    out.T_n = 2.0 * PI / out.omega_n;
    return out;
}

double RollDampingBuilder::estimateIeffFromMassGM(
    double mass,
    double GM,
    double omega_n)
{
    if (mass <= 0.0 || GM <= 0.0 || omega_n <= 0.0)
        throw std::runtime_error("Invalid mass/GM/omega_n for Ieff estimation.");

    // C44 = m*g*GM
    const double C44 = mass * G * GM;
    return C44 / (omega_n * omega_n);
}

RollViscDamping RollDampingBuilder::convertToPhysicalDamping(
    const DecayPolyCoeff& poly,
    double Ieff)
{
    if (Ieff <= 0.0 || poly.omega_n <= 0.0)
        throw std::runtime_error("Invalid Ieff or omega_n.");

    RollViscDamping out;

    // Delta phi = a*phi_m + b*phi_m^2 + c*phi_m^3
    // a = pi*B1/(2*Ieff*omega_n)
    // b = 4*B2/(3*Ieff)
    // c = 3*pi*omega_n*B3/(8*Ieff)

    out.B44_lin = (2.0 * Ieff * poly.omega_n / PI) * poly.a;
    out.B44_quad = (3.0 * Ieff / 4.0) * poly.b;
    out.B44_cube = (8.0 * Ieff / (3.0 * PI * poly.omega_n)) * poly.c;

    return out;
}

namespace
{
    // Absolute-value local maxima of a (t, phi) trajectory, with a minimum time
    // gap so numerical wiggles are not counted as peaks.
    std::vector<RollDecayPeak> absPeaksFromTraj(
        const std::vector<double>& t,
        const std::vector<double>& phi,
        double minDt)
    {
        std::vector<RollDecayPeak> pk;
        double lastT = -1.0e30;
        for (std::size_t i = 1; i + 1 < phi.size(); ++i)
        {
            const double a0 = std::abs(phi[i - 1]);
            const double a1 = std::abs(phi[i]);
            const double a2 = std::abs(phi[i + 1]);
            if (a1 >= a0 && a1 > a2 && (t[i] - lastT) >= minDt)
            {
                pk.push_back({ t[i], a1 });
                lastT = t[i];
            }
        }
        return pk;
    }
}

std::vector<RollDecayPeak> RollDampingBuilder::simulateDecayPeaks(
    const RollViscDamping& rv,
    double Ieff,
    double C44,
    double phi0,
    double phidot0,
    double tEnd)
{
    if (Ieff <= 0.0 || C44 <= 0.0 || tEnd <= 0.0)
        return {};

    const double omega_n = std::sqrt(C44 / Ieff);
    const double Tn = 2.0 * PI / omega_n;
    const double dt = std::min(0.01, Tn / 200.0);
    const int nStep = static_cast<int>(tEnd / dt) + 1;

    auto rhs = [&](double phi_now, double phidot_now)
    {
        const double Mvisc = rv.moment(phidot_now);
        const double phiddot = (Mvisc - C44 * phi_now) / Ieff;
        return std::pair<double, double>{ phidot_now, phiddot };
    };

    std::vector<double> tHist, phiHist;
    tHist.reserve(nStep);
    phiHist.reserve(nStep);

    double t = 0.0, phi = phi0, phidot = phidot0;
    for (int i = 0; i < nStep; ++i)
    {
        tHist.push_back(t);
        phiHist.push_back(phi);

        auto [k1p, k1v] = rhs(phi, phidot);
        auto [k2p, k2v] = rhs(phi + 0.5 * dt * k1p, phidot + 0.5 * dt * k1v);
        auto [k3p, k3v] = rhs(phi + 0.5 * dt * k2p, phidot + 0.5 * dt * k2v);
        auto [k4p, k4v] = rhs(phi + dt * k3p, phidot + dt * k3v);

        phi    += dt * (k1p + 2.0 * k2p + 2.0 * k3p + k4p) / 6.0;
        phidot += dt * (k1v + 2.0 * k2v + 2.0 * k3v + k4v) / 6.0;
        t      += dt;

        if (!std::isfinite(phi) || std::abs(phi) > 10.0 * std::abs(phi0) + 1.0)
            break; // diverged: stop, misfit will penalise
    }

    return absPeaksFromTraj(tHist, phiHist, 0.3 * Tn);
}

double RollDampingBuilder::envelopeMisfit(
    const std::vector<RollDecayPeak>& meas,
    const std::vector<RollDecayPeak>& sim,
    int skipFirstPeaks)
{
    const int s = std::max(0, skipFirstPeaks);
    const int nm = static_cast<int>(meas.size()) - s;
    const int ns = static_cast<int>(sim.size()) - s;
    if (nm < 3 || ns < 1)
        return 1.0e9; // not enough peaks (e.g. over-damped / diverged)

    const int n = std::min(nm, ns);
    double acc = 0.0;
    int cnt = 0;
    for (int k = 0; k < n; ++k)
    {
        const double am = meas[s + k].phi_abs;
        const double as = sim[s + k].phi_abs;
        if (am <= 1.0e-9) continue;
        const double e = (as - am) / am;
        acc += e * e;
        ++cnt;
    }
    if (cnt == 0)
        return 1.0e9;

    // Penalise a peak-count mismatch (sim decaying too fast/slow loses cycles).
    const double countPenalty = 0.05 * std::abs(nm - ns);
    return std::sqrt(acc / cnt) + countPenalty;
}

RollViscDamping RollDampingBuilder::refineByEnvelope(
    const std::vector<RollDecaySample>& samples,
    const RollViscDamping& init,
    double Ieff,
    double C44,
    int skipFirstPeaks)
{
    if (samples.size() < 4 || Ieff <= 0.0 || C44 <= 0.0)
        return init;

    const double phi0 = samples.front().phi;
    const double phidot0 = 0.0;
    const double tEnd = samples.back().t;

    const double omega_n = std::sqrt(C44 / Ieff);
    const double Tn = 2.0 * PI / omega_n;

    std::vector<double> tm, pm;
    tm.reserve(samples.size());
    pm.reserve(samples.size());
    for (const auto& s : samples) { tm.push_back(s.t); pm.push_back(s.phi); }
    const std::vector<RollDecayPeak> measPk = absPeaksFromTraj(tm, pm, 0.3 * Tn);

    auto evalMisfit = [&](double bl, double bq) -> double
    {
        RollViscDamping rv{ std::max(0.0, bl), std::max(0.0, bq), init.B44_cube };
        const auto simPk = simulateDecayPeaks(rv, Ieff, C44, phi0, phidot0, tEnd);
        return envelopeMisfit(measPk, simPk, skipFirstPeaks);
    };

    // Two-pass grid search: coarse over a generous box around the analytic
    // guess, then a fine box around the coarse optimum. Robust, derivative-free,
    // and cheap (smooth 2-D misfit).
    const int N = 31;
    double loL = 0.0, hiL = std::max(3.0 * init.B44_lin, 0.05);
    double loQ = 0.0, hiQ = std::max(5.0 * init.B44_quad, 0.05);

    double bestL = init.B44_lin, bestQ = init.B44_quad;
    double bestF = evalMisfit(bestL, bestQ);

    for (int pass = 0; pass < 2; ++pass)
    {
        const double stepL = (hiL - loL) / (N - 1);
        const double stepQ = (hiQ - loQ) / (N - 1);
        for (int i = 0; i < N; ++i)
        {
            const double bl = loL + i * stepL;
            for (int j = 0; j < N; ++j)
            {
                const double bq = loQ + j * stepQ;
                const double f = evalMisfit(bl, bq);
                if (f < bestF) { bestF = f; bestL = bl; bestQ = bq; }
            }
        }
        // Shrink box to +/- 2 cells around the current optimum for pass 2.
        loL = std::max(0.0, bestL - 2.0 * stepL);
        hiL = bestL + 2.0 * stepL;
        loQ = std::max(0.0, bestQ - 2.0 * stepQ);
        hiQ = bestQ + 2.0 * stepQ;
    }

    return RollViscDamping{ bestL, bestQ, init.B44_cube };
}

RollDecayIdentificationResult RollDampingBuilder::identifyFromDecay(
    const std::vector<RollDecaySample>& samples,
    double mass,
    double GM,
    int polyOrder,
    int minPeakGap,
    double IeffOverride,
    bool refine,
    int refineSkipFirstPeaks)
{
    RollDecayIdentificationResult out;
    out.peaks = extractAbsolutePeaks(samples, minPeakGap);
    out.poly = fitPeakDecrement(out.peaks, polyOrder);

    if (IeffOverride > 0.0)
        out.Ieff = IeffOverride;
    else
        out.Ieff = estimateIeffFromMassGM(mass, GM, out.poly.omega_n);

    out.damping = convertToPhysicalDamping(out.poly, out.Ieff);

    if (refine)
    {
        const double C44 = mass * G * GM;
        const RollViscDamping before = out.damping;
        const RollViscDamping after =
            refineByEnvelope(samples, before, out.Ieff, C44, refineSkipFirstPeaks);
        out.damping = after;

        std::cout << "[RollDamping] envelope refine (Option A):\n"
                  << "  before: B44_lin=" << before.B44_lin
                  << " B44_quad=" << before.B44_quad
                  << " B44_cube=" << before.B44_cube << "\n"
                  << "  after : B44_lin=" << after.B44_lin
                  << " B44_quad=" << after.B44_quad
                  << " B44_cube=" << after.B44_cube
                  << "  (skipFirstPeaks=" << refineSkipFirstPeaks << ")\n";
    }

    return out;
}

RollViscDamping RollDampingBuilder::build(
    const std::string& casePath,
    const RollDampingConfig& cfg,
    double mass,
    double GM)
{
    if (!cfg.enabled)
        return RollViscDamping{};

    if (cfg.mode == RollDampingMode::DirectCoefficients)
    {
        return RollViscDamping{
            cfg.direct.B44_lin,
            cfg.direct.B44_quad,
            cfg.direct.B44_cube
        };
    }

    if (cfg.mode == RollDampingMode::FromDecayCsv)
    {
        const auto samples = readDecayCsv(casePath + cfg.decay.csvPath, cfg.decay.angleInDeg);

        const auto result = identifyFromDecay(
            samples,
            mass,
            GM,
            cfg.decay.polyOrder,
            cfg.decay.minPeakGap,
            cfg.decay.IeffOverride,
            cfg.decay.refine,
            cfg.decay.refineSkipFirstPeaks
        );

        return result.damping;
    }

    throw std::runtime_error("Unknown RollDampingMode.");
}