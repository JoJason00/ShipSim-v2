#pragma once

#include "../../seakeeping/RadiationKernelCache.h"
#include "../../config/CaseConfig.h"
#include "../../config/SeakeepingConfig.h"
#include <string>
#include <vector>

class CoupledRadiationKernelRepo
{
public:
    CoupledRadiationKernelRepo(const ShipConfig& ship,
                               const SeakeepingConfig& skCfg,
                               const std::string& casePath);

    RadiationKernelData getKernel(double Fn, double dtFast) const;

    /// Median spacing of the distinct Fn nodes in the radiation library. Used by
    /// the adaptive refresh policy. \p valid is false when there are fewer than
    /// two distinct Fn nodes (caller should fall back to the fixed tolerance).
    double fnSpacingMedian(bool& valid) const;

private:
    struct FileEntry
    {
        std::string path;
        RadiationKernelData kernel;
    };

    static std::string modesTag(const SeakeepingConfig& skCfg);
    static std::string keyDouble(double x);
    static RadiationKernelData resampleKernel(const RadiationKernelData& src, double dtFast);

    void scan() const;

private:
    ShipConfig ship_;
    SeakeepingConfig skCfg_;
    std::string casePath_;
    mutable bool scanned_ = false;
    mutable std::vector<FileEntry> entries_;

    mutable bool fnSpacingCacheValid_ = false;
    mutable bool fnSpacingValid_ = false;
    mutable double fnSpacingMedian_ = 0.0;
};
