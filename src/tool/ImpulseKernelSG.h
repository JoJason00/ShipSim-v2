// ImpulseKernelSG.h — Savitzky–Golay smoothing for impulse / lag histories.
// Uses tool/SavitzkyGolay.h (namespace sg). Radiation Klag and future FK kernels
// can share smoothMatrixLagHistoryInPlace / writeMatrixLagCompareCsv.

#pragma once

#include "tool/SavitzkyGolay.h"

#include <Eigen/Dense>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace impulse_sg
{

struct Options
{
    int windowLength = 7;
    int polyOrder = 2;
    int preserveLeadingSamples = 4;
};

inline int minOddWindowForPoly(int polyOrder)
{
    int w = polyOrder + 2;
    w = sg::makeOdd(w);
    return w;
}

inline int clampSgWindow(int requested, int n, int polyOrder)
{
    if (n < 3)
        return -1;

    int w = sg::makeOdd(std::max(3, requested));
    const int wMin = minOddWindowForPoly(polyOrder);
    if (w < wMin)
        w = wMin;

    if (w > n)
    {
        w = n;
        if ((w % 2) == 0)
            --w;
        if (w < wMin)
            return -1;
    }
    return w;
}

inline Eigen::VectorXd smoothVectorSame(const Eigen::VectorXd& y, const Options& opt)
{
    const int n = static_cast<int>(y.size());
    if (n < 3)
        return y;

    const int w = clampSgWindow(opt.windowLength, n, opt.polyOrder);
    if (w < 0)
        return y;

    const Eigen::RowVectorXd c = sg::convolutionCoeffs(w, opt.polyOrder, 0);
    Eigen::VectorXd ys = sg::convolveSame(y, c);

    const int lead = std::clamp(opt.preserveLeadingSamples, 0, n);
    for (int m = 0; m < lead; ++m)
        ys(m) = y(m);

    return ys;
}

inline void smoothMatrixLagHistoryInPlace(
    std::vector<Eigen::MatrixXd>& hist,
    const Options& opt)
{
    if (hist.empty())
        return;

    const int rows = static_cast<int>(hist[0].rows());
    const int cols = static_cast<int>(hist[0].cols());
    const int M = static_cast<int>(hist.size());

    for (auto& K : hist)
    {
        if (static_cast<int>(K.rows()) != rows || static_cast<int>(K.cols()) != cols)
            throw std::runtime_error("impulse_sg::smoothMatrixLagHistoryInPlace: size mismatch");
    }

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            Eigen::VectorXd y(M);
            for (int m = 0; m < M; ++m)
                y(m) = hist[static_cast<std::size_t>(m)](i, j);

            const Eigen::VectorXd ys = smoothVectorSame(y, opt);

            for (int m = 0; m < M; ++m)
                hist[static_cast<std::size_t>(m)](i, j) = ys(m);
        }
    }
}

inline bool writeMatrixLagCompareCsv(
    const std::string& path,
    double dt,
    const std::vector<std::pair<int, int>>& modePairs0,
    const std::vector<Eigen::MatrixXd>& before,
    const std::vector<Eigen::MatrixXd>& after)
{
    if (before.size() != after.size() || before.empty())
        return false;

    const int M = static_cast<int>(before.size());
    const int rows = static_cast<int>(before[0].rows());
    const int cols = static_cast<int>(before[0].cols());

    for (int m = 0; m < M; ++m)
    {
        if (static_cast<int>(before[static_cast<std::size_t>(m)].rows()) != rows ||
            static_cast<int>(before[static_cast<std::size_t>(m)].cols()) != cols ||
            static_cast<int>(after[static_cast<std::size_t>(m)].rows()) != rows ||
            static_cast<int>(after[static_cast<std::size_t>(m)].cols()) != cols)
            return false;
    }

    std::vector<std::pair<int, int>> pairs;
    pairs.reserve(modePairs0.size());
    for (const auto& pr : modePairs0)
    {
        if (pr.first >= 0 && pr.first < rows && pr.second >= 0 && pr.second < cols)
            pairs.push_back(pr);
    }
    if (pairs.empty())
        return false;

    namespace fs = std::filesystem;
    const fs::path outPath(path);
    if (outPath.has_parent_path())
        fs::create_directories(outPath.parent_path());

    std::ofstream out(path);
    if (!out)
        return false;

    out << std::setprecision(17);
    out << "lag_index,time";
    for (const auto& pr : pairs)
        out << ",K_" << pr.first << "_" << pr.second << "_raw"
            << ",K_" << pr.first << "_" << pr.second << "_sg";
    out << "\n";

    for (int m = 0; m < M; ++m)
    {
        const double t = (dt > 0.0 && std::isfinite(dt)) ? static_cast<double>(m) * dt : static_cast<double>(m);
        out << m << "," << t;
        for (const auto& pr : pairs)
        {
            const int i = pr.first;
            const int j = pr.second;
            out << "," << before[static_cast<std::size_t>(m)](i, j)
                << "," << after[static_cast<std::size_t>(m)](i, j);
        }
        out << "\n";
    }

    return true;
}

} // namespace impulse_sg