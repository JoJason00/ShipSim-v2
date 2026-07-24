#include "ImpulseKernelWaveForceProvider.h"
#include "../common/CoupledMath.h"
#include "../../const/Const.h"

#include <algorithm>
#include <cmath>
#include <iostream>

ImpulseKernelWaveForceProvider::ImpulseKernelWaveForceProvider(
    const SeakeepingConfig& skCfg,
    const CoupledExcitationKernelRepo& excRepo,
    const HydroRefreshConfig& refreshCfg)
    : skCfg_(skCfg), excRepo_(excRepo), refreshCfg_(refreshCfg)
{
}

bool ImpulseKernelWaveForceProvider::needRefresh(
    double Fn, double betaRel, double omegaInc, double dtFast) const
{
    if (!cacheValid_)                              return true;
    if (std::abs(dtFast - dtCached_) > 1e-12)      return true;

    const double FnRef    = std::max(std::abs(FnCached_),    1e-3);
    const double omegaRef = std::max(std::abs(omegaCached_), 1e-3);

    // Fixed-tolerance policy (unchanged default behaviour).
    double fnTolAbs   = refreshCfg_.speedTolRatio * FnRef;
    double betaTolDeg = refreshCfg_.betaTolDeg;

    // Adaptive policy: derive Fn / heading tolerances from the actual library
    // grid spacing at this incident frequency, so a denser library refreshes
    // more often. Any axis with <2 distinct nodes keeps the fixed tolerance.
    if (refreshCfg_.adaptiveRefresh)
    {
        const double f = (refreshCfg_.adaptiveRefreshFactor > 1e-6)
                              ? refreshCfg_.adaptiveRefreshFactor
                              : 0.5;
        const auto gs = excRepo_.gridSpacing(omegaInc);
        if (gs.fnValid)
            fnTolAbs = std::max(1e-4, f * gs.dFn);
        if (gs.betaValid)
            betaTolDeg = std::max(0.1, f * gs.dBetaRad * 180.0 / PI);
    }

    if (std::abs(Fn       - FnCached_)    > fnTolAbs)                             return true;
    if (std::abs(omegaInc - omegaCached_) > refreshCfg_.omegaTolRatio * omegaRef) return true;

    const double dBetaDeg = std::abs(coupled::wrapPiSigned(betaRel - betaCached_)) * 180.0 / PI;
    if (dBetaDeg > betaTolDeg) return true;

    return false;
}

void ImpulseKernelWaveForceProvider::onWindowBegin(
    const CoupledSlowState3DOF& slow0,
    const CoupledEncounterState& enc0,
    double dtFast,
    const CoupledFastWindowInfo& /*win*/)
{
    if (!needRefresh(slow0.Fn, enc0.betaRel, enc0.omega, dtFast))
        return;

    excCached_   = excRepo_.getKernel(slow0.Fn, enc0.betaRel, enc0.omega, dtFast);
    FnCached_    = slow0.Fn;
    betaCached_  = enc0.betaRel;
    omegaCached_ = enc0.omega;
    dtCached_    = dtFast;
    cacheValid_  = true;

    // The offline excitation kernel can be stored at a coarse dt; resampling to
    // dtFast inflates the row count. A wave-excitation impulse response decays
    // to zero physically — truncate where ||Qlag.row(k)|| has fallen below a
    // fraction of the peak. Keep at least 5 encounter periods; never more than
    // 10 encounter periods (radiation-style memory is not needed here).
    if (excCached_.Qlag.rows() > 1)
    {
        const int nRows = static_cast<int>(excCached_.Qlag.rows());
        double peakNorm = 0.0;
        for (int k = 0; k < nRows; ++k)
            peakNorm = std::max(peakNorm, excCached_.Qlag.row(k).norm());

        const double tol = 1e-3 * std::max(peakNorm, 1e-30);

        int cutoff = nRows;
        for (int k = nRows - 1; k >= 1; --k)
        {
            if (excCached_.Qlag.row(k).norm() > tol)
            {
                cutoff = k + 1;
                break;
            }
        }

        const double omegaEnc = std::max(std::abs(enc0.kin.we), 1e-3);
        const double Tenc = 2.0 * PI / omegaEnc;
        const int minLags = std::min(nRows, static_cast<int>(std::ceil(5.0 * Tenc / dtFast)));
        const int maxLags = std::min(nRows, static_cast<int>(std::ceil(10.0 * Tenc / dtFast)));

        cutoff = std::max(cutoff, minLags);
        cutoff = std::min(cutoff, maxLags);

        if (cutoff < nRows)
        {
            Eigen::MatrixXd Qtrim = excCached_.Qlag.topRows(cutoff);
            excCached_.Qlag = std::move(Qtrim);
        }
    }

    const int newMem = std::max<int>(1, static_cast<int>(excCached_.Qlag.rows()));
    if (newMem != memoryLen_)
    {
        coupled::resizeHistoryPreservingEnd(etaHist_, newMem);
        memoryLen_ = newMem;
    }

    std::cout << "[ImpulseKernel] refresh kernel at Fn=" << slow0.Fn
              << " betaRel=" << (enc0.betaRel * 180.0 / PI) << "deg"
              << " omega=" << enc0.omega
              << " memoryLen=" << memoryLen_ << "\n";
}

void ImpulseKernelWaveForceProvider::setEtaHistory(const std::vector<double>& hist)
{
    etaHist_ = hist;
    if (memoryLen_ > 0)
        coupled::resizeHistoryPreservingEnd(etaHist_, memoryLen_);
}

Eigen::VectorXd ImpulseKernelWaveForceProvider::evaluate(
    const CoupledSlowState3DOF& /*slow*/,
    double /*t*/,
    const CoupledEncounterState& enc,
    Eigen::RowVectorXd* force6)
{
    const int D = skCfg_.DOF;
    Eigen::VectorXd F = Eigen::VectorXd::Zero(D);

    coupled::pushRing(etaHist_, enc.waveAmp * std::cos(enc.phaseAtCG));

    if (!cacheValid_ || excCached_.Qlag.rows() == 0 || etaHist_.empty())
        return F;

    const int nLag = std::min<int>(excCached_.Qlag.rows(), static_cast<int>(etaHist_.size()));
    const int last = static_cast<int>(etaHist_.size()) - 1;

    for (int k = 0; k < nLag; ++k)
    {
        const double etaLag = etaHist_[static_cast<std::size_t>(last - k)];
        F.noalias() += excCached_.Qlag.row(k).transpose() * etaLag;
    }

    if (force6) coupled::scatterToSixDof(F, skCfg_, *force6);
    return F;
}
