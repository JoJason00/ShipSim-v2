#pragma once

#include "../seakeeping/WaveForceRegion.h"

#include <string>
#include <Eigen/Dense>

namespace FKImpulseKernelIO
{
    struct Params
    {
        std::string casePath;
        std::string shipName;

        double Fn = 0.0;
        double dirRad = 0.0;
        double dt = 0.0;
        int tMot = 0;

        // Following seas split the kernel into three encounter-frequency
        // regions (F1/F2/F3); head / beam seas all use region "H".
        WaveForceRegion region = WaveForceRegion::Head;

        bool includeShipName = true;
    };

    struct Data
    {
        Eigen::MatrixXd fkForce;
        Eigen::MatrixXd dForce;
    };

    // Header-block metadata, parsed by loadFromFile.
    struct LoadInfo
    {
        double          Fn     = 0.0;
        double          dirRad = 0.0;
        WaveForceRegion region = WaveForceRegion::Head;
        double          dt     = 0.0;
        int             tMot   = 0;
    };

    std::string keyDouble(double x);

    std::string makePath(const Params& p);

    // Exact-match load (the file path must be identical to makePath(p)).
    bool load(const Params& p, Data& data);

    void save(const Params& p, const Data& data);

    // Library-style fuzzy search: scan {casePath}/fkImpulse/ for files matching
    // (ship, Fn, Dir, Lu impulse bucket H or 1/2/3) regardless of (dt, tMot).
    // When several candidates match, pick the longest memory window (dt × tMot).
    std::string findByKey(const std::string& casePath,
                          const std::string& shipName,
                          double Fn,
                          double dirRad,
                          WaveForceRegion region);

    // Load by direct file path. Fills both the matrices and the header-block
    // metadata so the caller knows what (dt, tMot) the file is on.
    bool loadFromFile(const std::string& path, LoadInfo& info, Data& data);

    // Resample (src.fkForce, src.dForce) from (srcDt, srcTMot) onto a new
    // (dstDt, dstTMot) grid centred at τ = 0. Linear interpolation in time;
    // when the destination range extends past the source, the kernel is padded
    // with its end-of-window value (assumed to have decayed by then).
    void resample(const Data& src, double srcDt, int srcTMot,
                  double dstDt, int dstTMot, Data& dst);
}