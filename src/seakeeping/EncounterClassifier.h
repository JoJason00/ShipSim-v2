#pragma once

#include "KernelKey.h"
#include "../wave/WaveComponent.h"

// Maps a wave component + ship state to the kernel it needs. Bucketisation
// tolerances mirror the coupled-side HydroRefreshConfig so straight-course
// runs collapse all components onto a handful of keys.
struct EncounterClassifierConfig
{
    double Lpp = 1.0;
    double fnTol = 0.03;        // Fn bucket width
    double betaTolRad = 5.0 * 3.141592653589793 / 180.0;
};

namespace encounter_classifier
{
    // theta is absolute wave propagation direction; psi is ship heading (rad).
    KernelKey classify(const WaveComponent& w,
                        double U,
                        double psi,
                        const EncounterClassifierConfig& cfg);
}
