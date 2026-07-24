#pragma once

#include <Eigen/Core>
#include <vector>

// ---------------------------------------------------------------------------
// ChiHistorySmooth
//
// Purpose: smooth panel-wise chi(t) histories used in radiation kernel
// construction, without bloating LinearCumminsTDGF.cpp.
//
// Rationale vs ad-hoc "detect spikes then SG locally":
//   - Global Savitzky¨CGolay (or Gaussian) along discrete time is stable,
//     reproducible, and preserves low-frequency memory content better than
//     aggressive local heuristics when the whole trace is mildly noisy.
//   - Optional robust pre-pass (Hampel) removes isolated outliers before SG.
//
// Typical use: replace smoothStartupHistory(...) with:
//   ChiHistorySmooth::smoothPanelHistory(raw, ChiHistorySmooth::Options::sgDefault());
// ---------------------------------------------------------------------------

namespace ChiHistorySmooth
{

    struct Options
    {
        enum class Method
        {
            None,           // copy-through
            Gaussian,       // discrete Gaussian kernel, reflect boundaries
            SavitzkyGolay,  // FIR from least-squares polynomial fit
            TwoPassGaussianThenSg // mild Gaussian then SG (good default when chi is very jagged)
        };

        Method method = Method::SavitzkyGolay;

        // Savitzky¨CGolay
        int sgWindow = 11;   // must be odd and >= polyDeg + 2
        int sgPolyDeg = 3;   // cubic is a good default

        // Gaussian (in steps; sigma is in units of time indices, not seconds)
        double gaussSigmaSteps = 2.0; // ~ kernel half-width ~ 3*sigma
        int gaussRadius = 0;        // if 0, auto = ceil(3*sigma)

        // Robust pre-pass (per panel, per time index) using local median/MAD.
        bool hampelPrewash = false;
        int hampelWindow = 11; // odd >= 3
        double hampelK = 3.0;  // threshold factor on MAD

        // When true, finite-sample boundaries use reflect padding before convolution.
        bool reflectBoundary = true;

        static Options sgDefault()
        {
            Options o;
            o.method = Method::SavitzkyGolay;
            o.sgWindow = 11;
            o.sgPolyDeg = 3;
            return o;
        }

        static Options gentleTwoPass()
        {
            Options o;
            o.method = Method::TwoPassGaussianThenSg;
            o.gaussSigmaSteps = 1.2;
            o.sgWindow = 9;
            o.sgPolyDeg = 3;
            return o;
        }
    };

    // Each entry of `raw` is chi at one time step; each VectorXd has NE panels.
    // `dt` is kept for API symmetry (future band-pass designs); currently unused.
    std::vector<Eigen::VectorXd> smoothPanelHistory(
        const std::vector<Eigen::VectorXd>& raw,
        double dt,
        const Options& opt);

} // namespace ChiHistorySmooth
