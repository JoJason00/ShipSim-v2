#pragma once

#include "../../seakeeping/TDGFExcitationKernelIO.h"
#include "../../config/CaseConfig.h"
#include "../../config/SeakeepingConfig.h"
#include <string>
#include <vector>

class CoupledExcitationKernelRepo
{
public:
    CoupledExcitationKernelRepo(const ShipConfig& ship,
                                const SeakeepingConfig& skCfg,
                                const std::string& casePath);

    TDGFExcitationKernelData getKernel(double Fn,
                                       double betaRel,
                                       double omegaIncident,
                                       double dtFast) const;

    struct GridSpacing
    {
        double dFn = 0.0;        ///< median spacing of distinct Fn nodes at this omega
        double dBetaRad = 0.0;   ///< median spacing of distinct betaRel nodes (rad)
        bool fnValid = false;    ///< false if fewer than two distinct Fn nodes
        bool betaValid = false;  ///< false if fewer than two distinct betaRel nodes
    };

    /// Library grid spacing among kernels matching \p omegaIncident (same exact
    /// match rule as getKernel). Used by the adaptive refresh policy. Cached per
    /// omega so it is cheap to call every window.
    GridSpacing gridSpacing(double omegaIncident) const;

private:
    struct FileEntry
    {
        std::string path;
        TDGFExcitationKernelData kernel;
    };

    static TDGFExcitationKernelData resampleKernel(const TDGFExcitationKernelData& src, double dtFast);
    static double wrapAngle(double x);

    void scan() const;

private:
    ShipConfig ship_;
    SeakeepingConfig skCfg_;
    std::string casePath_;
    mutable bool scanned_ = false;
    mutable std::vector<FileEntry> entries_;

    mutable bool gsCacheValid_ = false;
    mutable double gsCacheOmega_ = 0.0;
    mutable GridSpacing gsCache_{};
};
