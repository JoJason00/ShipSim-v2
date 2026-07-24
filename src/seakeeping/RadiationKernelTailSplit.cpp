//#include "RadiationKernelTailSplit.h"
//#include "RadiationKernelTailSplit.h"
////#include "LinearCumminsTDGF.h"   
//#include "RadiationKernelCache.h"
//
//#include <algorithm>
//#include <cmath>
//#include <fstream>
//#include <iomanip>
//#include <iostream>
//#include <stdexcept>
//#include <filesystem>
//
//namespace KernelPost
//{
//    namespace
//    {
//        bool containsMode(const std::vector<int>& modes, int mode)
//        {
//            return std::find(modes.begin(), modes.end(), mode) != modes.end();
//        }
//
//        bool isForcedPair(
//            const std::vector<std::pair<int, int>>& pairs,
//            int rowMode,
//            int colMode)
//        {
//            return std::find(
//                pairs.begin(),
//                pairs.end(),
//                std::make_pair(rowMode, colMode)) != pairs.end();
//        }
//
//        void writeReport(
//            const std::string& file,
//            const RadiationKernelData& kernel,
//            const TailSplitResult& r)
//        {
//            if (file.empty())
//                return;
//
//            std::ofstream out(file);
//            if (!out.is_open())
//            {
//                std::cerr << "[KernelTailSplit] cannot write report: "
//                    << file << "\n";
//                return;
//            }
//
//            const int D = kernel.DOF;
//
//            out << std::setprecision(17);
//
//            out << "# tailStart," << r.tailStart
//                << ",tailEnd," << r.tailEnd
//                << ",dt," << kernel.dt
//                << ",TG," << kernel.TG
//                << "\n";
//
//            out << "row,col,rowMode,colMode,"
//                << "peakAbs,tailMean,tailOsc,meanToPeak,oscToMean,"
//                << "selected,Kinf\n";
//
//            for (int i = 0; i < D; ++i)
//            {
//                for (int j = 0; j < D; ++j)
//                {
//                    out << i << ","
//                        << j << ","
//                        << kernel.modes[i] << ","
//                        << kernel.modes[j] << ","
//                        << r.peakAbs(i, j) << ","
//                        << r.tailMean(i, j) << ","
//                        << r.tailOsc(i, j) << ","
//                        << r.meanToPeak(i, j) << ","
//                        << r.oscToMean(i, j) << ","
//                        << r.selected(i, j) << ","
//                        << r.Kinf(i, j) << "\n";
//                }
//            }
//
//            std::cout << "[KernelTailSplit] report written: "
//                << file << "\n";
//        }
//    }
//
//    TailSplitResult splitTailConstantInPlace(
//        RadiationKernelData& kernel,
//        const TailSplitOptions& opt,
//        const std::string& reportFile)
//    {
//        TailSplitResult r;
//
//        const int D = kernel.DOF;
//        const int TG = std::min(
//            kernel.TG,
//            static_cast<int>(kernel.Klag.size()) - 1);
//
//        r.Kinf = Eigen::MatrixXd::Zero(D, D);
//        r.peakAbs = Eigen::MatrixXd::Zero(D, D);
//        r.tailMean = Eigen::MatrixXd::Zero(D, D);
//        r.tailOsc = Eigen::MatrixXd::Zero(D, D);
//        r.meanToPeak = Eigen::MatrixXd::Zero(D, D);
//        r.oscToMean = Eigen::MatrixXd::Zero(D, D);
//        r.selected = Eigen::MatrixXi::Zero(D, D);
//
//        if (!opt.enabled)
//            return r;
//
//        if (D <= 0)
//            throw std::runtime_error("KernelTailSplit: invalid DOF.");
//
//        if (TG < 5)
//            throw std::runtime_error("KernelTailSplit: kernel history is too short.");
//
//        if (static_cast<int>(kernel.modes.size()) != D)
//            throw std::runtime_error("KernelTailSplit: modes size mismatch.");
//
//        if (static_cast<int>(kernel.Klag.size()) < TG + 1)
//            throw std::runtime_error("KernelTailSplit: Klag size mismatch.");
//
//        for (int m = 0; m <= TG; ++m)
//        {
//            if (kernel.Klag[static_cast<std::size_t>(m)].rows() != D ||
//                kernel.Klag[static_cast<std::size_t>(m)].cols() != D)
//            {
//                throw std::runtime_error("KernelTailSplit: Klag matrix size mismatch.");
//            }
//        }
//
//        if (kernel.C_prime.rows() != D || kernel.C_prime.cols() != D)
//            throw std::runtime_error("KernelTailSplit: C_prime size mismatch.");
//
//        const double tailFraction =
//            std::max(0.01, std::min(0.80, opt.tailFraction));
//
//        int tailCount = static_cast<int>(
//            std::ceil(tailFraction * static_cast<double>(TG + 1)));
//
//        tailCount = std::max(tailCount, opt.minTailSamples);
//        tailCount = std::max(tailCount, 3);
//        tailCount = std::min(tailCount, TG);
//
//        const int i0 = TG - tailCount + 1;
//
//        r.tailStart = i0;
//        r.tailEnd = TG;
//
//        constexpr double eps = 1.0e-30;
//
//        for (int i = 0; i < D; ++i)
//        {
//            const int rowMode = kernel.modes[i];
//
//            for (int j = 0; j < D; ++j)
//            {
//                const int colMode = kernel.modes[j];
//
//                const bool inTargetBlock =
//                    !opt.onlyTargetBlock ||
//                    (containsMode(opt.targetModes, rowMode) &&
//                        containsMode(opt.targetModes, colMode));
//
//                if (!inTargetBlock &&
//                    !isForcedPair(opt.forcedModePairs, rowMode, colMode))
//                {
//                    continue;
//                }
//
//                double peak = 0.0;
//
//                // TimeToFrequency �����߾������Ǵ� lag=1 ��ʼʹ�� Klag��
//                // ��������Ҳ�� 1 ��ʼͳ�ơ�
//                for (int m = 1; m <= TG; ++m)
//                {
//                    const double v =
//                        kernel.Klag[static_cast<std::size_t>(m)](i, j);
//
//                    peak = std::max(peak, std::abs(v));
//                }
//
//                double sum = 0.0;
//                int count = 0;
//
//                for (int m = i0; m <= TG; ++m)
//                {
//                    sum += kernel.Klag[static_cast<std::size_t>(m)](i, j);
//                    ++count;
//                }
//
//                const double mean =
//                    (count > 0) ? sum / static_cast<double>(count) : 0.0;
//
//                double osc = 0.0;
//
//                for (int m = i0; m <= TG; ++m)
//                {
//                    const double v =
//                        kernel.Klag[static_cast<std::size_t>(m)](i, j);
//
//                    osc = std::max(osc, std::abs(v - mean));
//                }
//
//                const double meanRatio =
//                    std::abs(mean) / std::max(eps, peak);
//
//                const double oscRatio =
//                    osc / std::max(eps, std::abs(mean));
//
//                r.peakAbs(i, j) = peak;
//                r.tailMean(i, j) = mean;
//                r.tailOsc(i, j) = osc;
//                r.meanToPeak(i, j) = meanRatio;
//                r.oscToMean(i, j) = oscRatio;
//
//                bool pick =
//                    peak > eps &&
//                    meanRatio >= opt.minMeanToPeak &&
//                    oscRatio <= opt.maxOscToMean;
//
//                if (isForcedPair(opt.forcedModePairs, rowMode, colMode))
//                {
//                    // ǿ�ƴ���ʱ��ȻҪ���Ǵ���ֵ������
//                    pick = peak > eps && meanRatio >= 0.01;
//                }
//
//                if (pick)
//                {
//                    r.selected(i, j) = 1;
//                    r.Kinf(i, j) = mean;
//                }
//            }
//        }
//
//        const double maxKinf = r.Kinf.cwiseAbs().maxCoeff();
//
//        if (maxKinf > 0.0)
//        {
//            // �����д���Լ����m=0 ������ TimeToFrequency��Ҳ���������߾�����
//            // ����ֻ���� m>=1��
//            for (int m = 1; m <= TG; ++m)
//            {
//                kernel.Klag[static_cast<std::size_t>(m)] -= r.Kinf;
//            }
//
//            if (!kernel.Klag.empty())
//            {
//                kernel.Klag[0].setZero();
//            }
//
//            if (opt.addTailToCPrime)
//            {
//                kernel.C_prime += r.Kinf;
//            }
//
//            r.applied = true;
//        }
//
//        if (opt.verbose)
//        {
//            std::cout << "\n[KernelTailSplit]\n"
//                << "  TG = " << TG << "\n"
//                << "  dt = " << kernel.dt << "\n"
//                << "  tail index = [" << r.tailStart
//                << ", " << r.tailEnd << "]\n"
//                << "  addTailToCPrime = "
//                << (opt.addTailToCPrime ? "true" : "false") << "\n";
//
//            bool any = false;
//
//            for (int i = 0; i < D; ++i)
//            {
//                for (int j = 0; j < D; ++j)
//                {
//                    if (r.selected(i, j) == 0)
//                        continue;
//
//                    any = true;
//
//                    std::cout << "  select Kinf("
//                        << kernel.modes[i] << ","
//                        << kernel.modes[j] << ") = "
//                        << r.Kinf(i, j)
//                        << "  mean/peak = "
//                        << r.meanToPeak(i, j)
//                        << "  osc/mean = "
//                        << r.oscToMean(i, j)
//                        << "\n";
//                }
//            }
//
//            if (!any)
//            {
//                std::cout << "  no tail constant selected.\n";
//            }
//
//            std::cout << "\n";
//        }
//
//        writeReport(reportFile, kernel, r);
//
//        return r;
//    }
//
//
//
//    void writeKernelHistoryCsv(
//        const RadiationKernelData& kernel,
//        const KernelScaleInfo& scale,
//        const std::string& outFile,
//        bool writeDimensional,
//        bool writeNondimensional)
//    {
//        if (outFile.empty())
//            return;
//
//        if (!writeDimensional && !writeNondimensional)
//            return;
//
//        const int D = kernel.DOF;
//        const int TG = std::min(
//            kernel.TG,
//            static_cast<int>(kernel.Klag.size()) - 1);
//
//        if (D <= 0 || TG < 0)
//            throw std::runtime_error("writeKernelHistoryCsv: invalid kernel size.");
//
//        if (static_cast<int>(kernel.modes.size()) != D)
//            throw std::runtime_error("writeKernelHistoryCsv: modes size mismatch.");
//
//        if (scale.L <= 0.0 || scale.displacement <= 0.0 ||
//            scale.rho <= 0.0 || scale.g <= 0.0)
//        {
//            throw std::runtime_error("writeKernelHistoryCsv: invalid scale info.");
//        }
//
//        namespace fs = std::filesystem;
//        const fs::path p(outFile);
//        if (p.has_parent_path())
//            fs::create_directories(p.parent_path());
//
//        std::ofstream out(outFile);
//        if (!out.is_open())
//        {
//            throw std::runtime_error("writeKernelHistoryCsv: cannot open " + outFile);
//        }
//
//        const double L = scale.L;
//        const double V = scale.displacement;
//        const double forceScale = scale.rho * scale.g * V;
//        const double momentScale = scale.rho * scale.g * V * L;
//        const double tScale = std::sqrt(scale.g / L);
//
//        auto rowScaleForMode = [&](int mode) -> double
//        {
//            // 0,1,2 ��ƽ�����̣�����������
//            // 3,4,5 ��ת�����̣�������������
//            return (mode >= 0 && mode <= 2) ? forceScale : momentScale;
//        };
//
//        auto colScaleForMode = [&](int mode) -> double
//        {
//            // �������� K_ij * qdot_j��
//            // ƽ��λ������Ҫ�� L��ת����λ���в��� L��
//            return (mode >= 0 && mode <= 2) ? L : 1.0;
//        };
//
//        out << std::setprecision(17);
//
//        out << "step,t,t_nd";
//
//        for (int i = 0; i < D; ++i)
//        {
//            for (int j = 0; j < D; ++j)
//            {
//                const int mi = kernel.modes[i];
//                const int mj = kernel.modes[j];
//
//                if (writeDimensional)
//                    out << ",Kdim_" << mi << mj;
//
//                if (writeNondimensional)
//                    out << ",Knd_" << mi << mj;
//            }
//        }
//
//        out << "\n";
//
//        for (int m = 0; m <= TG; ++m)
//        {
//            const double t = static_cast<double>(m) * kernel.dt;
//            const double t_nd = t * tScale;
//
//            out << m << "," << t << "," << t_nd;
//
//            const Eigen::MatrixXd& K =
//                kernel.Klag[static_cast<std::size_t>(m)];
//
//            for (int i = 0; i < D; ++i)
//            {
//                const int rowMode = kernel.modes[i];
//                const double rowScale = rowScaleForMode(rowMode);
//
//                for (int j = 0; j < D; ++j)
//                {
//                    const int colMode = kernel.modes[j];
//                    const double colScale = colScaleForMode(colMode);
//
//                    const double Kdim = K(i, j);
//                    const double Knd = Kdim * colScale / rowScale;
//
//                    if (writeDimensional)
//                        out << "," << Kdim;
//
//                    if (writeNondimensional)
//                        out << "," << Knd;
//                }
//            }
//
//            out << "\n";
//        }
//
//        std::cout << "[KernelTailSplit] kernel history written: "
//            << outFile << "\n";
//    }
//
//
//    TailSplitResult applyTailSplitWorkflowInPlace(
//        RadiationKernelData& kernel,
//        const KernelScaleInfo& scale,
//        const TailSplitOptions& splitOpt,
//        const TailSplitWorkflowOptions& workflowOpt)
//    {
//        TailSplitResult result;
//
//        if (!workflowOpt.enabled)
//            return result;
//
//        namespace fs = std::filesystem;
//
//        fs::path outDir = workflowOpt.outputDir.empty()
//            ? fs::path(".")
//            : fs::path(workflowOpt.outputDir);
//
//        fs::create_directories(outDir);
//
//        const std::string tag = workflowOpt.tag.empty()
//            ? std::string("case")
//            : workflowOpt.tag;
//
//        const fs::path beforeFile =
//            outDir / ("kernel_history_before_tail_split_" + tag + ".csv");
//
//        const fs::path afterFile =
//            outDir / ("kernel_history_after_tail_split_" + tag + ".csv");
//
//        const fs::path reportFile =
//            outDir / ("kernel_tail_split_" + tag + ".csv");
//
//        if (workflowOpt.writeBeforeHistory)
//        {
//            writeKernelHistoryCsv(
//                kernel,
//                scale,
//                beforeFile.string(),
//                workflowOpt.writeDimensional,
//                workflowOpt.writeNondimensional);
//        }
//
//        result = splitTailConstantInPlace(
//            kernel,
//            splitOpt,
//            reportFile.string());
//
//        if (workflowOpt.writeAfterHistory)
//        {
//            writeKernelHistoryCsv(
//                kernel,
//                scale,
//                afterFile.string(),
//                workflowOpt.writeDimensional,
//                workflowOpt.writeNondimensional);
//        }
//
//        return result;
//    }
//}



#include "RadiationKernelTailSplit.h"
//#include "LinearCumminsTDGF.h"   
#include "RadiationKernelCache.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <filesystem>

namespace KernelPost
{
    namespace
    {
        bool containsMode(const std::vector<int>& modes, int mode)
        {
            return std::find(modes.begin(), modes.end(), mode) != modes.end();
        }

        bool isForcedPair(
            const std::vector<std::pair<int, int>>& pairs,
            int rowMode,
            int colMode)
        {
            return std::find(
                pairs.begin(),
                pairs.end(),
                std::make_pair(rowMode, colMode)) != pairs.end();
        }

        void writeReport(
            const std::string& file,
            const RadiationKernelData& kernel,
            const TailSplitResult& r)
        {
            if (file.empty())
                return;

            std::ofstream out(file);
            if (!out.is_open())
            {
                std::cerr << "[KernelTailSplit] cannot write report: "
                    << file << "\n";
                return;
            }

            const int D = kernel.DOF;

            out << std::setprecision(17);

            out << "# tailStart," << r.tailStart
                << ",tailEnd," << r.tailEnd
                << ",dt," << kernel.dt
                << ",TG," << kernel.TG
                << "\n";

            out << "row,col,rowMode,colMode,"
                << "peakAbs,tailMean,tailOsc,meanToPeak,oscToMean,"
                << "selected,Kinf\n";

            for (int i = 0; i < D; ++i)
            {
                for (int j = 0; j < D; ++j)
                {
                    out << i << ","
                        << j << ","
                        << kernel.modes[i] << ","
                        << kernel.modes[j] << ","
                        << r.peakAbs(i, j) << ","
                        << r.tailMean(i, j) << ","
                        << r.tailOsc(i, j) << ","
                        << r.meanToPeak(i, j) << ","
                        << r.oscToMean(i, j) << ","
                        << r.selected(i, j) << ","
                        << r.Kinf(i, j) << "\n";
                }
            }

            std::cout << "[KernelTailSplit] report written: "
                << file << "\n";
        }
    }

    TailSplitResult splitTailConstantInPlace(
        RadiationKernelData& kernel,
        const TailSplitOptions& opt,
        const std::string& reportFile)
    {
        TailSplitResult r;

        const int D = kernel.DOF;
        const int TG = std::min(
            kernel.TG,
            static_cast<int>(kernel.Klag.size()) - 1);

        r.Kinf = Eigen::MatrixXd::Zero(D, D);
        r.peakAbs = Eigen::MatrixXd::Zero(D, D);
        r.tailMean = Eigen::MatrixXd::Zero(D, D);
        r.tailOsc = Eigen::MatrixXd::Zero(D, D);
        r.meanToPeak = Eigen::MatrixXd::Zero(D, D);
        r.oscToMean = Eigen::MatrixXd::Zero(D, D);
        r.selected = Eigen::MatrixXi::Zero(D, D);

        if (!opt.enabled)
            return r;

        if (D <= 0)
            throw std::runtime_error("KernelTailSplit: invalid DOF.");

        if (TG < 5)
            throw std::runtime_error("KernelTailSplit: kernel history is too short.");

        if (static_cast<int>(kernel.modes.size()) != D)
            throw std::runtime_error("KernelTailSplit: modes size mismatch.");

        if (static_cast<int>(kernel.Klag.size()) < TG + 1)
            throw std::runtime_error("KernelTailSplit: Klag size mismatch.");

        for (int m = 0; m <= TG; ++m)
        {
            if (kernel.Klag[static_cast<std::size_t>(m)].rows() != D ||
                kernel.Klag[static_cast<std::size_t>(m)].cols() != D)
            {
                throw std::runtime_error("KernelTailSplit: Klag matrix size mismatch.");
            }
        }

        if (kernel.C_prime.rows() != D || kernel.C_prime.cols() != D)
            throw std::runtime_error("KernelTailSplit: C_prime size mismatch.");

        const double tailFraction =
            std::max(0.01, std::min(0.80, opt.tailFraction));

        int tailCount = static_cast<int>(
            std::ceil(tailFraction * static_cast<double>(TG + 1)));

        tailCount = std::max(tailCount, opt.minTailSamples);
        tailCount = std::max(tailCount, 3);
        tailCount = std::min(tailCount, TG);

        const int i0 = TG - tailCount + 1;

        r.tailStart = i0;
        r.tailEnd = TG;

        constexpr double eps = 1.0e-30;

        for (int i = 0; i < D; ++i)
        {
            const int rowMode = kernel.modes[i];

            for (int j = 0; j < D; ++j)
            {
                const int colMode = kernel.modes[j];

                const bool inTargetBlock =
                    !opt.onlyTargetBlock ||
                    (containsMode(opt.targetModes, rowMode) &&
                        containsMode(opt.targetModes, colMode));

                if (!inTargetBlock &&
                    !isForcedPair(opt.forcedModePairs, rowMode, colMode))
                {
                    continue;
                }

                double peak = 0.0;

                // TimeToFrequency �����߾������Ǵ� lag=1 ��ʼʹ�� Klag��
                // ��������Ҳ�� 1 ��ʼͳ�ơ�
                for (int m = 1; m <= TG; ++m)
                {
                    const double v =
                        kernel.Klag[static_cast<std::size_t>(m)](i, j);

                    peak = std::max(peak, std::abs(v));
                }

                double sum = 0.0;
                int count = 0;

                for (int m = i0; m <= TG; ++m)
                {
                    sum += kernel.Klag[static_cast<std::size_t>(m)](i, j);
                    ++count;
                }

                const double mean =
                    (count > 0) ? sum / static_cast<double>(count) : 0.0;

                double osc = 0.0;

                for (int m = i0; m <= TG; ++m)
                {
                    const double v =
                        kernel.Klag[static_cast<std::size_t>(m)](i, j);

                    osc = std::max(osc, std::abs(v - mean));
                }

                const double meanRatio =
                    std::abs(mean) / std::max(eps, peak);

                const double oscRatio =
                    osc / std::max(eps, std::abs(mean));

                r.peakAbs(i, j) = peak;
                r.tailMean(i, j) = mean;
                r.tailOsc(i, j) = osc;
                r.meanToPeak(i, j) = meanRatio;
                r.oscToMean(i, j) = oscRatio;

                bool pick =
                    peak > eps &&
                    meanRatio >= opt.minMeanToPeak &&
                    oscRatio <= opt.maxOscToMean;

                if (isForcedPair(opt.forcedModePairs, rowMode, colMode))
                {
                    // ǿ�ƴ���ʱ��ȻҪ���Ǵ���ֵ������
                    pick = peak > eps && meanRatio >= 0.01;
                }

                if (pick)
                {
                    r.selected(i, j) = 1;
                    r.Kinf(i, j) = mean;
                }
            }
        }

        const double maxKinf = r.Kinf.cwiseAbs().maxCoeff();

        if (maxKinf > 0.0)
        {
            // Tail split �󣬲������˶���Ϊ��
            //     K_res(t) = K(t) - Kinf
            //
            // m >= 1 �����߼�������� TimeToFrequency ����Ҫʹ�����䣬
            // �� m = 0 Ҳ�������㡣�������� K ���ݿ��ز���ʱ��
            // t = 0 �����ᱻ��Ϊѹ�ͣ������𲽶κ˺�����������
            for (int m = 1; m <= TG; ++m)
            {
                kernel.Klag[static_cast<std::size_t>(m)] -= r.Kinf;
            }

            if (!kernel.Klag.empty())
            {
                // ���� K(0) �� m>=1 �Ĳ���˴���һ�¡�
                kernel.Klag[0] -= r.Kinf;
                kernel.K0 = kernel.Klag[0];
            }

            if (opt.addTailToCPrime)
            {
                kernel.C_prime += r.Kinf;
            }

            r.applied = true;
        }

        if (opt.verbose)
        {
            std::cout << "\n[KernelTailSplit]\n"
                << "  TG = " << TG << "\n"
                << "  dt = " << kernel.dt << "\n"
                << "  tail index = [" << r.tailStart
                << ", " << r.tailEnd << "]\n"
                << "  addTailToCPrime = "
                << (opt.addTailToCPrime ? "true" : "false") << "\n";

            bool any = false;

            for (int i = 0; i < D; ++i)
            {
                for (int j = 0; j < D; ++j)
                {
                    if (r.selected(i, j) == 0)
                        continue;

                    any = true;

                    std::cout << "  select Kinf("
                        << kernel.modes[i] << ","
                        << kernel.modes[j] << ") = "
                        << r.Kinf(i, j)
                        << "  mean/peak = "
                        << r.meanToPeak(i, j)
                        << "  osc/mean = "
                        << r.oscToMean(i, j)
                        << "\n";
                }
            }

            if (!any)
            {
                std::cout << "  no tail constant selected.\n";
            }

            std::cout << "\n";
        }

        writeReport(reportFile, kernel, r);

        return r;
    }



    void writeKernelHistoryCsv(
        const RadiationKernelData& kernel,
        const KernelScaleInfo& scale,
        const std::string& outFile,
        bool writeDimensional,
        bool writeNondimensional)
    {
        if (outFile.empty())
            return;

        if (!writeDimensional && !writeNondimensional)
            return;

        const int D = kernel.DOF;
        const int TG = std::min(
            kernel.TG,
            static_cast<int>(kernel.Klag.size()) - 1);

        if (D <= 0 || TG < 0)
            throw std::runtime_error("writeKernelHistoryCsv: invalid kernel size.");

        if (static_cast<int>(kernel.modes.size()) != D)
            throw std::runtime_error("writeKernelHistoryCsv: modes size mismatch.");

        if (scale.L <= 0.0 || scale.displacement <= 0.0 ||
            scale.rho <= 0.0 || scale.g <= 0.0)
        {
            throw std::runtime_error("writeKernelHistoryCsv: invalid scale info.");
        }

        namespace fs = std::filesystem;
        const fs::path p(outFile);
        if (p.has_parent_path())
            fs::create_directories(p.parent_path());

        std::ofstream out(outFile);
        if (!out.is_open())
        {
            throw std::runtime_error("writeKernelHistoryCsv: cannot open " + outFile);
        }

        const double L = scale.L;
        const double V = scale.displacement;
        const double forceScale = scale.rho * scale.g * V;
        const double momentScale = scale.rho * scale.g * V * L;
        const double tScale = std::sqrt(scale.g / L);

        auto rowScaleForMode = [&](int mode) -> double
        {
            // 0,1,2 ��ƽ�����̣�����������
            // 3,4,5 ��ת�����̣�������������
            return (mode >= 0 && mode <= 2) ? forceScale : momentScale;
        };

        auto colScaleForMode = [&](int mode) -> double
        {
            // �������� K_ij * qdot_j��
            // ƽ��λ������Ҫ�� L��ת����λ���в��� L��
            return (mode >= 0 && mode <= 2) ? L : 1.0;
        };

        out << std::setprecision(17);

        out << "step,t,t_nd";

        for (int i = 0; i < D; ++i)
        {
            for (int j = 0; j < D; ++j)
            {
                const int mi = kernel.modes[i];
                const int mj = kernel.modes[j];

                if (writeDimensional)
                    out << ",Kdim_" << mi << mj;

                if (writeNondimensional)
                    out << ",Knd_" << mi << mj;
            }
        }

        out << "\n";

        // 非均匀路径下 Klag_times 存了真实节点时间；等距路径下 Klag_times[m] = m*kernel.dt
        // （由 build 或 cache 加载兜底）。两条路径都走 Klag_times，结果与之前对等距完全一致。
        const bool haveTimes =
            kernel.Klag_times.size() == kernel.Klag.size()
            && !kernel.Klag_times.empty();

        for (int m = 0; m <= TG; ++m)
        {
            const double t = haveTimes
                ? kernel.Klag_times[static_cast<std::size_t>(m)]
                : (static_cast<double>(m) * kernel.dt);
            const double t_nd = t * tScale;

            out << m << "," << t << "," << t_nd;

            const Eigen::MatrixXd& K =
                kernel.Klag[static_cast<std::size_t>(m)];

            for (int i = 0; i < D; ++i)
            {
                const int rowMode = kernel.modes[i];
                const double rowScale = rowScaleForMode(rowMode);

                for (int j = 0; j < D; ++j)
                {
                    const int colMode = kernel.modes[j];
                    const double colScale = colScaleForMode(colMode);

                    const double Kdim = K(i, j);
                    const double Knd = Kdim * colScale / rowScale;

                    if (writeDimensional)
                        out << "," << Kdim;

                    if (writeNondimensional)
                        out << "," << Knd;
                }
            }

            out << "\n";
        }

        std::cout << "[KernelTailSplit] kernel history written: "
            << outFile << "\n";
    }


    TailSplitResult applyTailSplitWorkflowInPlace(
        RadiationKernelData& kernel,
        const KernelScaleInfo& scale,
        const TailSplitOptions& splitOpt,
        const TailSplitWorkflowOptions& workflowOpt)
    {
        TailSplitResult result;

        if (!workflowOpt.enabled)
            return result;

        namespace fs = std::filesystem;

        fs::path outDir = workflowOpt.outputDir.empty()
            ? fs::path(".")
            : fs::path(workflowOpt.outputDir);

        fs::create_directories(outDir);

        const std::string tag = workflowOpt.tag.empty()
            ? std::string("case")
            : workflowOpt.tag;

        const fs::path beforeFile =
            outDir / ("kernel_history_before_tail_split_" + tag + ".csv");

        const fs::path afterFile =
            outDir / ("kernel_history_after_tail_split_" + tag + ".csv");

        const fs::path reportFile =
            outDir / ("kernel_tail_split_" + tag + ".csv");

        if (workflowOpt.writeBeforeHistory)
        {
            writeKernelHistoryCsv(
                kernel,
                scale,
                beforeFile.string(),
                workflowOpt.writeDimensional,
                workflowOpt.writeNondimensional);
        }

        result = splitTailConstantInPlace(
            kernel,
            splitOpt,
            reportFile.string());

        // ������ͬ�����������̽����󣬱�֤ K0 �� Klag[0] һ�¡�
        if (!kernel.Klag.empty())
        {
            kernel.K0 = kernel.Klag[0];
        }

        if (workflowOpt.writeAfterHistory)
        {
            writeKernelHistoryCsv(
                kernel,
                scale,
                afterFile.string(),
                workflowOpt.writeDimensional,
                workflowOpt.writeNondimensional);
        }

        return result;
    }
}