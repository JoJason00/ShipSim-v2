#pragma once

#include <Eigen/Dense>

#include <string>
#include <vector>
#include <utility>

struct RadiationKernelData;

namespace KernelPost
{
    struct KernelScaleInfo
    {
        double L = 0.0;
        double displacement = 0.0;
        double rho = 0.0;
        double g = 0.0;
    };

    struct TailSplitOptions
    {
        bool enabled = true;
        bool addTailToCPrime = true;

        double tailFraction = 0.20;
        int minTailSamples = 30;

        double minMeanToPeak = 0.05;
        double maxOscToMean = 0.70;

        std::vector<int> targetModes = { 2, 4 };
        bool onlyTargetBlock = true;

        std::vector<std::pair<int, int>> forcedModePairs;

        bool verbose = true;
    };

    struct TailSplitResult
    {
        Eigen::MatrixXd Kinf;

        Eigen::MatrixXd peakAbs;
        Eigen::MatrixXd tailMean;
        Eigen::MatrixXd tailOsc;
        Eigen::MatrixXd meanToPeak;
        Eigen::MatrixXd oscToMean;

        Eigen::MatrixXi selected;

        int tailStart = 0;
        int tailEnd = 0;

        bool applied = false;
    };

    TailSplitResult splitTailConstantInPlace(
        RadiationKernelData& kernel,
        const TailSplitOptions& opt,
        const std::string& reportFile = "");

    void writeKernelHistoryCsv(
        const RadiationKernelData& kernel,
        const KernelScaleInfo& scale,
        const std::string& outFile,
        bool writeDimensional = true,
        bool writeNondimensional = true);


    struct TailSplitWorkflowOptions
    {
        bool enabled = true;

        // 是否输出修正前/修正后的 K 时历
        bool writeBeforeHistory = true;
        bool writeAfterHistory = true;

        // 输出 Kdim / Knd
        bool writeDimensional = true;
        bool writeNondimensional = true;

        // 输出目录，例如 filePath/kernel_cache
        std::string outputDir;

        // 文件标签，例如 Fn0p3000
        std::string tag = "case";
    };

    TailSplitResult applyTailSplitWorkflowInPlace(
        RadiationKernelData& kernel,
        const KernelScaleInfo& scale,
        const TailSplitOptions& splitOpt,
        const TailSplitWorkflowOptions& workflowOpt);
}