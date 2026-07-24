#include "EncounterClassifier.h"
#include "../const/Const.h"
#include <cmath>

namespace
{
    double wrapPi(double x)
    {
        while (x <= -PI) x += 2.0 * PI;
        while (x >  PI)  x -= 2.0 * PI;
        return x;
    }
}

namespace encounter_classifier
{
    KernelKey classify(const WaveComponent& w,
                        double U,
                        double psi,
                        const EncounterClassifierConfig& cfg)
    {
        const double betaRel = wrapPi(w.theta - psi);

        KernelKey key;
        key.region = wave_force_region::classify(U, betaRel, w.omega);

        const double Fn = (cfg.Lpp > 0.0)
            ? U / std::sqrt(G * cfg.Lpp)
            : 0.0;
        key.fnBucket =
            static_cast<int>(std::lround(Fn / std::max(cfg.fnTol, 1e-9)));
        key.betaBucket =
            static_cast<int>(std::lround(betaRel / std::max(cfg.betaTolRad, 1e-9)));
        return key;
    }
}
