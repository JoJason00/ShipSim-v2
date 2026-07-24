#include "WaveSpectrum.h"
#include "../const/Const.h"

#include <cmath>
#include <random>
#include <algorithm>
#include <stdexcept>

namespace
{
    // PM, eq.(4): one-parameter (wind speed at 19.5 m).
    double pmDensity(double Uwind, double w)
    {
        if (w <= 0.0 || Uwind <= 0.0) return 0.0;
        const double a = 0.0081;
        const double beta = 0.74;
        const double g2 = G * G;
        return a * g2 / std::pow(w, 5.0)
             * std::exp(-beta * std::pow(G / (Uwind * w), 4.0));
    }

    // ITTC / ISSC two-parameter, eq.(5).
    double ittcDensity(double Hs, double T1, double w)
    {
        if (w <= 0.0 || T1 <= 0.0) return 0.0;
        const double T1_4 = std::pow(T1, 4.0);
        return 173.0 * Hs * Hs / (T1_4 * std::pow(w, 5.0))
             * std::exp(-691.0 / (T1_4 * std::pow(w, 4.0)));
    }

    // JONSWAP, eq.(6), parameterised by (Hs, Tp, gamma): ωm = 2π/Tp and the
    // energy scale is fixed by the Hasselmann Hs-normalisation so 4√m0 ≈ Hs.
    double jonswapDensity(double Hs, double Tp, double gamma, double w)
    {
        if (w <= 0.0 || Tp <= 0.0) return 0.0;
        const double wm = 2.0 * PI / Tp;
        const double sigma = (w <= wm) ? 0.07 : 0.09;
        const double r = std::exp(-((w - wm) * (w - wm))
                                  / (2.0 * sigma * sigma * wm * wm));
        const double norm = 1.0 - 0.287 * std::log(gamma);
        const double base = (5.0 / 16.0) * Hs * Hs * std::pow(wm, 4.0)
                          * std::pow(w, -5.0)
                          * std::exp(-1.25 * std::pow(wm / w, 4.0));
        return norm * base * std::pow(gamma, r);
    }

    // One Ochi-Hubble system (eq.(7) summand).
    double ohSystem(double Hsj, double wmj, double lamj, double w)
    {
        if (w <= 0.0 || wmj <= 0.0 || lamj <= 0.0) return 0.0;
        const double c = (4.0 * lamj + 1.0) / 4.0 * std::pow(wmj, 4.0);
        const double coef = std::pow(c, lamj) / std::tgamma(lamj);
        return 0.25 * coef * Hsj * Hsj
             / std::pow(w, 4.0 * lamj + 1.0)
             * std::exp(-(4.0 * lamj + 1.0) / 4.0 * std::pow(wmj / w, 4.0));
    }
}

namespace wave_spectrum
{
    double density(const SpectrumParams& p, double w)
    {
        switch (p.type)
        {
        case SpectrumType::PM:
            return pmDensity(p.Uwind, w);
        case SpectrumType::ITTC:
            return ittcDensity(p.Hs, p.T1, w);
        case SpectrumType::JONSWAP:
            return jonswapDensity(p.Hs, p.Tp, p.gamma, w);
        case SpectrumType::OH:
        {
            double s = 0.0;
            if (p.ohMode != 2) s += ohSystem(p.Hs1, p.wm1, p.lam1, w);
            if (p.ohMode != 1) s += ohSystem(p.Hs2, p.wm2, p.lam2, w);
            return s;
        }
        }
        return 0.0;
    }

    double m0(const SpectrumParams& p, int nQuad)
    {
        if (nQuad < 2 || p.omegaHi <= p.omegaLo) return 0.0;
        const double dw = (p.omegaHi - p.omegaLo) / nQuad;
        double s = 0.5 * (density(p, p.omegaLo) + density(p, p.omegaHi));
        for (int i = 1; i < nQuad; ++i)
            s += density(p, p.omegaLo + i * dw);
        return s * dw;
    }

    std::vector<WaveComponent> discretize(const SpectrumParams& p, double theta)
    {
        if (p.nComponents < 1 || p.omegaHi <= p.omegaLo)
            throw std::runtime_error("WaveSpectrum: invalid discretisation range/count.");

        const int M = std::max(2000, 20 * p.nComponents);
        const double dw = (p.omegaHi - p.omegaLo) / M;
        std::vector<double> w(M + 1), C(M + 1);
        w[0] = p.omegaLo;
        C[0] = 0.0;
        double prevS = density(p, p.omegaLo);
        for (int i = 1; i <= M; ++i)
        {
            w[i] = p.omegaLo + i * dw;
            const double curS = density(p, w[i]);
            C[i] = C[i - 1] + 0.5 * (prevS + curS) * dw;
            prevS = curS;
        }
        const double mTot = C[M];
        if (mTot <= 0.0)
            throw std::runtime_error("WaveSpectrum: spectrum has no energy (check params/range).");

        const int N = p.nComponents;
        const double Eband = mTot / N;

        auto omegaAtEnergy = [&](double Et) -> double
        {
            auto it = std::lower_bound(C.begin(), C.end(), Et);
            if (it == C.begin()) return w.front();
            if (it == C.end())   return w.back();
            const int i1 = static_cast<int>(it - C.begin());
            const int i0 = i1 - 1;
            const double f = (Et - C[i0]) / std::max(C[i1] - C[i0], 1e-30);
            return w[i0] + f * (w[i1] - w[i0]);
        };

        std::mt19937_64 rng(p.seed);
        std::uniform_real_distribution<double> uPhase(0.0, 2.0 * PI);
        std::uniform_real_distribution<double> u01(0.0, 1.0);

        const double a = std::sqrt(2.0 * Eband);

        std::vector<WaveComponent> comps;
        comps.reserve(N);
        for (int k = 0; k < N; ++k)
        {
            const double wLo = omegaAtEnergy(k * Eband);
            const double wHi = omegaAtEnergy((k + 1) * Eband);
            const double wj = p.randomFreqInBand
                ? wLo + u01(rng) * (wHi - wLo)
                : omegaAtEnergy((k + 0.5) * Eband);
            comps.push_back(WaveComponent{ a, wj, theta, uPhase(rng) });
        }
        return comps;
    }
}
