#include "WaveField.h"
#include "../const/Const.h"
#include <cmath>

WaveField WaveField::makeRegular(double amp, double dirRad, double omega,
                                 double eps)
{
    return WaveField(std::vector<WaveComponent>{
        WaveComponent{ amp, omega, dirRad, eps } });
}

WaveField WaveField::makeIrregular(const SpectrumParams& sp, double dirRad)
{
    return WaveField(wave_spectrum::discretize(sp, dirRad));
}

WaveField WaveField::makeCross(const WaveField& a, const WaveField& b)
{
    std::vector<WaveComponent> c = a.comps_;
    c.insert(c.end(), b.comps_.begin(), b.comps_.end());
    return WaveField(std::move(c));
}

double WaveField::elevationAt(double x, double y, double t) const
{
    double eta = 0.0;
    for (const auto& w : comps_)
    {
        const double k = w.omega * w.omega / G;
        const double phase =
            k * (x * std::cos(w.theta) + y * std::sin(w.theta))
            - w.omega * t + w.eps;
        eta += w.a * std::cos(phase);
    }
    return eta;
}
