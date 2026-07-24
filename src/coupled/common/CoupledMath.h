#pragma once

// Small shared utilities for the coupled pipeline.
// Header-only on purpose — these are tiny, performance-sensitive operations
// shared across the wave-force providers and a few other call sites.

#include "../../const/Const.h"
#include "../../config/SeakeepingConfig.h"

#include <Eigen/Dense>
#include <algorithm>
#include <vector>

namespace coupled
{
    // Wrap angle into (-π, π].
    inline double wrapPiSigned(double x)
    {
        while (x <= -PI) x += 2.0 * PI;
        while (x >   PI) x -= 2.0 * PI;
        return x;
    }

    // Append `x` to the back of a ring buffer of fixed length, dropping the oldest.
    // No-op when the buffer is empty.
    inline void pushRing(std::vector<double>& hist, double x)
    {
        if (hist.empty()) return;
        std::rotate(hist.begin(), hist.begin() + 1, hist.end());
        hist.back() = x;
    }

    // Resize `hist` to exactly `size` entries while preserving the most-recent
    // values: pads zeros at the front if growing, drops from the front if shrinking.
    inline void resizeHistoryPreservingEnd(std::vector<double>& hist, int size)
    {
        if (size <= 0) { hist.clear(); return; }
        const int sz = static_cast<int>(hist.size());
        if      (sz < size) hist.insert(hist.begin(), size - sz, 0.0);
        else if (sz > size) hist.erase(hist.begin(), hist.begin() + (sz - size));
    }

    // Relative wave heading β_rel → library lookup on [0,π] and port-side mirror.
    // betaRelDeg is the true (-180,180] deg value; betaQueryDeg is the 0–180° table key.
    inline void waveHeadingLookup(double betaRel,
                                  double& lookupDirRad,
                                  bool& mirrorPS,
                                  double& betaRelDeg,
                                  double& betaQueryDeg)
    {
        betaRelDeg = betaRel * 180.0 / PI;
        lookupDirRad = betaRel;
        mirrorPS     = false;
        while (lookupDirRad < 0.0)        lookupDirRad += 2.0 * PI;
        while (lookupDirRad >= 2.0 * PI)  lookupDirRad -= 2.0 * PI;
        betaQueryDeg = lookupDirRad * 180.0 / PI;
        if (betaQueryDeg > 180.0 + 1.0e-9)
        {
            betaQueryDeg = 360.0 - betaQueryDeg;
            mirrorPS     = true;
        }
    }

    // Scatter a DOF-length vector ordered by skCfg.modes into a fixed
    // 6-component (surge, sway, heave, roll, pitch, yaw) row vector.
    inline void scatterToSixDof(const Eigen::VectorXd& F,
                                const SeakeepingConfig& skCfg,
                                Eigen::RowVectorXd& force6)
    {
        force6.setZero(6);
        const int D = skCfg.DOF;
        for (int j = 0; j < D; ++j)
        {
            const int mode = skCfg.modes[j];
            if (mode >= 0 && mode < 6)
                force6(mode) = F(j);
        }
    }
}
