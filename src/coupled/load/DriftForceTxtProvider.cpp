#include "DriftForceTxtProvider.h"
#include "../../const/Const.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <limits>
#include <cmath>

namespace
{
    double parseSpeedFromFilename(const std::string& stem)
    {
        // �������� DriftForce_V0 / DriftForce_V2 / DriftForce_V20
        const auto pos = stem.find('V');
        if (pos == std::string::npos || pos + 1 >= stem.size())
            throw std::runtime_error("DriftForceTxtProvider: cannot parse speed from filename: " + stem);

        std::string num;
        for (size_t i = pos + 1; i < stem.size(); ++i)
        {
            const char c = stem[i];
            if (std::isdigit(static_cast<unsigned char>(c)) || c == '.' || c == '-' || c == '+')
                num.push_back(c);
            else
                break;
        }

        if (num.empty())
            throw std::runtime_error("DriftForceTxtProvider: empty speed tag in filename: " + stem);

        // ��ǰ first-pass �ٶ��ļ������ Vxx ��������ٶ� s.U ͬ��λ
        return std::stod(num);
    }

    bool isTxtFile(const std::filesystem::path& p)
    {
        const auto ext = p.extension().string();
        return ext == ".txt" || ext == ".TXT";
    }
}

double DriftForceTxtProvider::wrap360(double deg)
{
    while (deg < 0.0)   deg += 360.0;
    while (deg >= 360.) deg -= 360.0;
    return deg;
}

int DriftForceTxtProvider::lowerBracket(const std::vector<double>& x, double v)
{
    if (x.size() < 2) return 0;
    if (v <= x.front()) return 0;
    if (v >= x.back())  return static_cast<int>(x.size()) - 2;

    auto it = std::upper_bound(x.begin(), x.end(), v);
    int i1 = static_cast<int>(std::distance(x.begin(), it));
    return std::max(0, i1 - 1);
}

double DriftForceTxtProvider::lerp(double a, double b, double t)
{
    return (1.0 - t) * a + t * b;
}

Eigen::Vector4d DriftForceTxtProvider::lerpVec(const Eigen::Vector4d& a, const Eigen::Vector4d& b, double t)
{
    return (1.0 - t) * a + t * b;
}

DriftForceTxtProvider::DriftForceTxtProvider(const std::string& folderPath)
{
    loadFolder(folderPath);
    if (tables_.empty())
        throw std::runtime_error("DriftForceTxtProvider: no valid drift-force txt files found in " + folderPath);
}

void DriftForceTxtProvider::loadFolder(const std::string& folderPath)
{
    namespace fs = std::filesystem;

    if (!fs::exists(folderPath))
        throw std::runtime_error("DriftForceTxtProvider: folder does not exist: " + folderPath);

    for (const auto& ent : fs::directory_iterator(folderPath))
    {
        if (!ent.is_regular_file()) continue;
        if (!isTxtFile(ent.path())) continue;

        const std::string stem = ent.path().stem().string();
        const double Uref = parseSpeedFromFilename(stem);
        loadOneFile(ent.path().string(), Uref);
    }

    std::sort(tables_.begin(), tables_.end(),
        [](const DriftTableAtSpeed& a, const DriftTableAtSpeed& b)
    {
        return a.Uref < b.Uref;
    });
}

void DriftForceTxtProvider::loadOneFile(const std::string& filePath, double Uref)
{
    std::ifstream in(filePath);
    if (!in.is_open())
        throw std::runtime_error("DriftForceTxtProvider: cannot open file " + filePath);

    int nBeta = 0, nOmega = 0;
    in >> nBeta >> nOmega;
    if (!in || nBeta <= 0 || nOmega <= 0)
        throw std::runtime_error("DriftForceTxtProvider: invalid header in " + filePath);

    DriftTableAtSpeed tab;
    tab.Uref = Uref;
    tab.betaDegAxis.resize(nBeta);
    tab.omegaAxis.resize(nOmega);
    tab.qXYN.assign(nBeta, std::vector<Eigen::Vector4d>(nOmega, Eigen::Vector4d::Zero()));

    for (int ib = 0; ib < nBeta; ++ib)
    {
        double betaDeg = 0.0;
        in >> betaDeg;
        if (!in)
            throw std::runtime_error("DriftForceTxtProvider: failed to read beta label in " + filePath);

        tab.betaDegAxis[ib] = betaDeg;

        for (int io = 0; io < nOmega; ++io)
        {
            double omega = 0.0;
            double qx = 0.0, qy = 0.0, qz = 0.0, qk = 0.0, qm = 0.0, qn = 0.0;
            in >> omega >> qx >> qy >> qz >> qk >> qm >> qn;
            if (!in)
                throw std::runtime_error("DriftForceTxtProvider: failed to read data block in " + filePath);

            if (ib == 0)
                tab.omegaAxis[io] = omega;

            // (X, Y, N, K) — 末位 K 为横摇漂移矩（之前被丢弃，现在保留）。
            tab.qXYN[ib][io] = Eigen::Vector4d(qx, qy, qn, qk);
        }
    }

    tables_.push_back(std::move(tab));
}

const DriftForceTxtProvider::DriftTableAtSpeed&
DriftForceTxtProvider::pickNearestSpeedTable(double U) const
{
    if (tables_.empty())
        throw std::runtime_error("DriftForceTxtProvider: no tables loaded.");

    const DriftTableAtSpeed* best = &tables_.front();
    double bestMetric = std::abs(best->Uref - U);

    for (const auto& t : tables_)
    {
        const double m = std::abs(t.Uref - U);
        if (m < bestMetric)
        {
            bestMetric = m;
            best = &t;
        }
    }
    return *best;
}

Eigen::Vector4d DriftForceTxtProvider::evalAtSpeed(double U, double betaDeg, double omega) const
{
    if (tables_.empty())
        return Eigen::Vector4d::Zero();

    // tables_ is sorted ascending by Uref (see loadFolder).
    if (tables_.size() == 1 || U <= tables_.front().Uref)
        return evalAtTable(tables_.front(), betaDeg, omega);
    if (U >= tables_.back().Uref)
        return evalAtTable(tables_.back(), betaDeg, omega);

    // Find the bracketing pair [iLo, iLo+1] with tables_[iLo].Uref <= U < tables_[iLo+1].Uref.
    std::size_t iLo = 0;
    for (std::size_t i = 0; i + 1 < tables_.size(); ++i)
    {
        if (tables_[i].Uref <= U && U < tables_[i + 1].Uref)
        {
            iLo = i;
            break;
        }
    }
    const auto& tL = tables_[iLo];
    const auto& tH = tables_[iLo + 1];
    const double dU = tH.Uref - tL.Uref;
    const double tU = (dU > 1.0e-12) ? std::clamp((U - tL.Uref) / dU, 0.0, 1.0) : 0.0;

    const Eigen::Vector4d qL = evalAtTable(tL, betaDeg, omega);
    const Eigen::Vector4d qH = evalAtTable(tH, betaDeg, omega);
    return lerpVec(qL, qH, tU);
}

Eigen::Vector4d DriftForceTxtProvider::evalAtTable(
    const DriftTableAtSpeed& tab,
    double betaDeg,
    double omega) const
{
    if (tab.betaDegAxis.empty() || tab.omegaAxis.empty())
        return Eigen::Vector4d::Zero();

    // �����ֻ�� 0~180�������ô������ҶԳ���չ�� 0~360
    double signY = 1.0;
    double signN = 1.0;
    double betaQuery = wrap360(betaDeg);

    const double betaMax = tab.betaDegAxis.back();
    if (betaMax <= 180.0 + 1.0e-8)
    {
        if (betaQuery > 180.0)
        {
            betaQuery = 360.0 - betaQuery;
            signY = -1.0;
            signN = -1.0;
        }
    }

    const int ib = lowerBracket(tab.betaDegAxis, betaQuery);
    const int io = lowerBracket(tab.omegaAxis, omega);

    const double b0 = tab.betaDegAxis[ib];
    const double b1 = tab.betaDegAxis[std::min(ib + 1, static_cast<int>(tab.betaDegAxis.size()) - 1)];
    const double w0 = tab.omegaAxis[io];
    const double w1 = tab.omegaAxis[std::min(io + 1, static_cast<int>(tab.omegaAxis.size()) - 1)];

    const double tb = (std::abs(b1 - b0) > 1.0e-12) ? (betaQuery - b0) / (b1 - b0) : 0.0;
    const double tw = (std::abs(w1 - w0) > 1.0e-12) ? (omega - w0) / (w1 - w0) : 0.0;

    const Eigen::Vector4d q00 = tab.qXYN[ib][io];
    const Eigen::Vector4d q01 = tab.qXYN[ib][std::min(io + 1, static_cast<int>(tab.omegaAxis.size()) - 1)];
    const Eigen::Vector4d q10 = tab.qXYN[std::min(ib + 1, static_cast<int>(tab.betaDegAxis.size()) - 1)][io];
    const Eigen::Vector4d q11 = tab.qXYN[std::min(ib + 1, static_cast<int>(tab.betaDegAxis.size()) - 1)]
        [std::min(io + 1, static_cast<int>(tab.omegaAxis.size()) - 1)];

    const Eigen::Vector4d q0 = lerpVec(q00, q01, std::clamp(tw, 0.0, 1.0));
    const Eigen::Vector4d q1 = lerpVec(q10, q11, std::clamp(tw, 0.0, 1.0));
    Eigen::Vector4d q = lerpVec(q0, q1, std::clamp(tb, 0.0, 1.0));

    q(1) *= signY;
    q(2) *= signN;
    q(3) *= signN;   // roll moment K mirrors like Y/N under port/starboard reflection
    return q;
}

CoupledExternalLoads3DOF DriftForceTxtProvider::sample(
    const CoupledSlowState3DOF& s,
    const CoupledEncounterState& enc) const
{
    if (tables_.empty())
        return CoupledExternalLoads3DOF{};

    Eigen::Vector4d acc = Eigen::Vector4d::Zero();

    // Speed is common to all components; q is linearly interpolated between the
    // two speed tables bracketing s.U (per-component β, ω still use the table's
    // own 2D bilerp). Avoids the step jump that nearest-table lookup would
    // introduce when U crosses the midpoint of two adjacent V*.txt tables
    // during e.g. a turning manoeuvre.
    if (!enc.comps.empty())
    {
        // C2: mean second-order drift = linear superposition over components,
        //   F2 = Σ_j q(β_rel,j, ω_j) · a_j²
        // (difference-frequency cross terms neglected — standard mean-drift).
        // Regular wave = 1 component == the scalar path -> bit-identical.
        for (const auto& w : enc.comps)
        {
            const double betaDeg =
                wrap360((w.theta - s.psi) * 180.0 / PI);
            const Eigen::Vector4d q = evalAtSpeed(s.U, betaDeg, w.omega);
            acc += q * (w.a * w.a);
        }
    }
    else
    {
        // Legacy scalar fallback (no component decomposition available).
        const double betaDeg = wrap360(enc.betaRel * 180.0 / PI);
        const Eigen::Vector4d q = evalAtSpeed(s.U, betaDeg, enc.omega);
        acc = q * (enc.waveAmp * enc.waveAmp);
    }

    CoupledExternalLoads3DOF out;
    out.X2 = acc(0);
    out.Y2 = acc(1);
    out.N2 = acc(2);
    out.K2 = acc(3);   // mean second-order roll (heel) moment

    out.X = out.X2;
    out.Y = out.Y2;
    out.N = out.N2;
    out.hasSecondOrder = true;
    return out;
}