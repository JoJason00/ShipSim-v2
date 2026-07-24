#pragma once

// One regular component of a (possibly irregular / multi-directional) sea.
// Irregular and cross seas are just a list of these — the impulse-response
// machinery never needs to know which sea type produced them.
struct WaveComponent
{
    double a = 0.0;      // amplitude [m]
    double omega = 0.0;  // incident circular frequency [rad/s]
    double theta = 0.0;  // absolute propagation direction [rad]
    double eps = 0.0;    // random phase [rad]
};
