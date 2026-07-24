#pragma once

#include "WaveForceRegion.h"
#include <cstddef>
#include <functional>

// Identifies which excitation impulse kernel a wave component needs. Two
// components with the same key share a kernel and are summed into one
// convolution (convolution is linear); different keys -> separate convolutions.
//
//  - fnBucket   : forward-speed bucket (kernel depends on U)
//  - betaBucket : relative-heading bucket (kernel depends on β_rel)
//  - region     : Head for non-following seas (kernel is ω-independent there);
//                  F1/F2/F3 for following / stern-quartering seas, where ω_e(ω)
//                  is non-monotonic and the time-domain kernel splits by band.
struct KernelKey
{
    int fnBucket = 0;
    int betaBucket = 0;
    WaveForceRegion region = WaveForceRegion::Head;

    bool operator==(const KernelKey& o) const
    {
        return fnBucket == o.fnBucket
            && betaBucket == o.betaBucket
            && region == o.region;
    }
};

namespace std
{
    template <> struct hash<KernelKey>
    {
        std::size_t operator()(const KernelKey& k) const noexcept
        {
            std::size_t h = static_cast<std::size_t>(k.fnBucket) * 73856093u;
            h ^= static_cast<std::size_t>(k.betaBucket) * 19349663u;
            h ^= static_cast<std::size_t>(k.region) * 83492791u;
            return h;
        }
    };
}
