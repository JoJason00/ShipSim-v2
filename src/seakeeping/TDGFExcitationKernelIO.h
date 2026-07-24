#pragma once

#include "WaveForceRegion.h"

#include <Eigen/Dense>
#include <string>
#include <vector>

struct TDGFExcitationKernelData
{
    double Fn = 0.0;
    double U = 0.0;
    double betaRel = 0.0;
    double omegaIncident = 0.0;
    double omegaEncounter = 0.0;
    double dt = 0.0;

    int DOF = 0;
    std::vector<int> modes;

    // Following seas split into three encounter-frequency regions (F1/F2/F3).
    // Head/beam seas stay in region H. Recorded as metadata for sanity checks
    // (the file-name lookup is still by (Fn, betaRel, omegaInc) since those
    // already pin the region uniquely).
    WaveForceRegion region = WaveForceRegion::Head;

    // rows = nLag, cols = DOF
    Eigen::MatrixXd Qlag;
};

class TDGFExcitationKernelIO
{
public:
    static bool save(const std::string& file, const TDGFExcitationKernelData& data);
    static bool load(const std::string& file, TDGFExcitationKernelData& data);
};
