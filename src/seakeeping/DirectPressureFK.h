#pragma once

#include <Eigen/Dense>
#include "Element.h"

struct DirectPressureFKContext
{
    double amp   = 0.0;   // component amplitude a [m]
    double omega = 0.0;   // incident circular frequency ω (k = ω²/g)
    double beta  = 0.0;   // relative heading β = θ_abs − ψ [rad] (body spatial)
    // Wave-field phase at the BODY ORIGIN, traced along the real trajectory:
    //   Φ0(t) = −k(x_e cosθ_abs + y_e sinθ_abs) − ω t + ε
    // The per-panel phase is Φ0 − k(x_b cosβ + y_b sinβ); integrating in body
    // coordinates yields the force directly in the heave/roll/pitch frame.
    double phiOrigin = 0.0;
};

class DirectPressureFK
{
public:
    static Eigen::RowVectorXd computeForce6(
        const Element& element,
        const DirectPressureFKContext& ctx,
        double t
    );
};