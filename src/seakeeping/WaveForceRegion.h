#pragma once

// WaveForceRegion — classify a wave excitation case into one of four regions
// based on (U, β_rel, ω_inc) per Lu et al. 2024 (J. Ship Mech. 28(7)).
//
// In head / beam / bow-quartering seas (cos β ≤ 0) the dispersion relation
// ω_e(ω_inc) is monotonic; all such cases share region tag "H".
//
// In stern-quartering / following seas (cos β > 0) ω_e(ω_inc) is a parabola
// with peak ω_e^* = g / (4 U cos β) at ω_inc^* = g / (2 U cos β). Lu et al. (2024)
// distinguish three structural bands (F1/F2/F3) for the time-domain kernel.
//
// FK / diffraction impulse **library** files are one per Lu bucket: head and
// bow-quartering / beam (cos β ≤ 0) share tag "H"; following / stern-quartering
// use numeric bucket tags "1","2","3" (F1,F2,F3) — not one file per incident
// frequency within the same bucket.

#include "../const/Const.h"

#include <cmath>
#include <string>

enum class WaveForceRegion
{
    Head = 0,
    F1   = 1,
    F2   = 2,
    F3   = 3
};

namespace wave_force_region
{
    // Returns true when (U, betaRel) puts us in the following-seas branch.
    inline bool inFollowingSeas(double U, double betaRel)
    {
        if (U <= 0.0) return false;
        return std::cos(betaRel) > 1e-6;
    }

    // Critical incident frequencies for following seas. cb must be > 0.
    inline double omegaC1(double U, double cb) { return G / (2.0 * U * cb); }
    inline double omegaC2(double U, double cb) { return G / (U * cb);       }

    // Classify (U, β_rel, ω_inc) → region.
    inline WaveForceRegion classify(double U, double betaRel, double omegaInc)
    {
        if (!inFollowingSeas(U, betaRel))
            return WaveForceRegion::Head;

        const double cb  = std::cos(betaRel);
        const double wc1 = omegaC1(U, cb);
        const double wc2 = omegaC2(U, cb);

        if (omegaInc < wc1) return WaveForceRegion::F1;
        if (omegaInc < wc2) return WaveForceRegion::F2;
        return WaveForceRegion::F3;
    }

    inline std::string tag(WaveForceRegion r)
    {
        switch (r)
        {
        case WaveForceRegion::F1:   return "F1";
        case WaveForceRegion::F2:   return "F2";
        case WaveForceRegion::F3:   return "F3";
        case WaveForceRegion::Head:
        default:                    return "H";
        }
    }

    // Short bucket id for FK impulse filenames / headers (following: 1,2,3).
    inline std::string impulseBucketTag(WaveForceRegion r)
    {
        switch (r)
        {
        case WaveForceRegion::F1:   return "1";
        case WaveForceRegion::F2:   return "2";
        case WaveForceRegion::F3:   return "3";
        case WaveForceRegion::Head:
        default:                    return "H";
        }
    }

    inline WaveForceRegion fromTag(const std::string& s)
    {
        if (s == "F1" || s == "1") return WaveForceRegion::F1;
        if (s == "F2" || s == "2") return WaveForceRegion::F2;
        if (s == "F3" || s == "3") return WaveForceRegion::F3;
        return WaveForceRegion::Head;
    }
}
