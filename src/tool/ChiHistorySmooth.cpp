#include "ChiHistorySmooth.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ChiHistorySmooth
{
    namespace
    {

        bool isOddPositive(int x)
        {
            return x > 0 && (x % 2) == 1;
        }

        Eigen::VectorXd savitzkyGolayConvCoeffs(int window, int polyDeg)
        {
            if (!isOddPositive(window))
                throw std::invalid_argument("ChiHistorySmooth: sgWindow must be a positive odd integer.");

            if (polyDeg < 0 || polyDeg >= window)
                throw std::invalid_argument("ChiHistorySmooth: need polyDeg < sgWindow.");

            const int m = (window - 1) / 2;
            Eigen::MatrixXd A(window, polyDeg + 1);

            for (int k = 0; k < window; ++k)
            {
                const double z = static_cast<double>(k - m);
                double zp = 1.0;
                for (int p = 0; p <= polyDeg; ++p)
                {
                    A(k, p) = zp;
                    zp *= z;
                }
            }

            const Eigen::MatrixXd AtA = A.transpose() * A;
            Eigen::LDLT<Eigen::MatrixXd> ldlt(AtA);
            if (ldlt.info() != Eigen::Success)
                throw std::runtime_error("ChiHistorySmooth: LDLT failed building SG coefficients.");

            const Eigen::VectorXd e0 = Eigen::VectorXd::Unit(polyDeg + 1, 0);

            // b = (A^T A)^{-1} A^T y, smoothed value is b0 = e0^T b = w^T y
            // => w = A * (A^T A)^{-1} e0
            const Eigen::VectorXd invAtA_e0 = ldlt.solve(e0);
            const Eigen::VectorXd w = A * invAtA_e0;

            return w;
        }

        int reflectedIndex(int i, int N)
        {
            if (N <= 1)
                return 0;

            while (i < 0 || i >= N)
            {
                if (i < 0)
                    i = -i - 1;
                if (i >= N)
                    i = 2 * N - i - 1;
            }
            return i;
        }

        double samplePadded(
            const Eigen::VectorXd& y,
            int idx,
            bool reflect)
        {
            const int N = static_cast<int>(y.size());
            if (!reflect)
                return y(std::max(0, std::min(N - 1, idx)));

            return y(reflectedIndex(idx, N));
        }

        Eigen::VectorXd convolveFIR(
            const Eigen::VectorXd& y,
            const Eigen::VectorXd& kernel,
            bool reflect)
        {
            const int N = static_cast<int>(y.size());
            const int K = static_cast<int>(kernel.size());
            const int m = (K - 1) / 2;

            Eigen::VectorXd out(N);
            for (int i = 0; i < N; ++i)
            {
                double acc = 0.0;
                for (int k = 0; k < K; ++k)
                {
                    const int j = i + (k - m);
                    acc += kernel(k) * samplePadded(y, j, reflect);
                }
                out(i) = acc;
            }
            return out;
        }

        Eigen::VectorXd gaussianKernel(double sigmaSteps, int radius)
        {
            if (sigmaSteps <= 0.0)
                throw std::invalid_argument("ChiHistorySmooth: gaussSigmaSteps must be > 0.");

            int r = radius;
            if (r <= 0)
                r = static_cast<int>(std::ceil(3.0 * sigmaSteps));
            r = std::max(1, r);

            const int K = 2 * r + 1;
            Eigen::VectorXd k(K);
            double sum = 0.0;
            const double s2 = 2.0 * sigmaSteps * sigmaSteps;

            for (int i = 0; i < K; ++i)
            {
                const double x = static_cast<double>(i - r);
                const double w = std::exp(-(x * x) / s2);
                k(i) = w;
                sum += w;
            }

            if (sum <= 0.0)
                throw std::runtime_error("ChiHistorySmooth: Gaussian kernel degenerate.");

            k /= sum;
            return k;
        }

        Eigen::VectorXd hampel1D(const Eigen::VectorXd& y, int win, double kThresh)
        {
            if (!isOddPositive(win))
                throw std::invalid_argument("ChiHistorySmooth: hampelWindow must be odd and >= 3.");

            const int N = static_cast<int>(y.size());
            const int m = (win - 1) / 2;
            Eigen::VectorXd out = y;

            std::vector<double> buf(static_cast<std::size_t>(win));
            std::vector<double> dev(static_cast<std::size_t>(win));

            for (int i = 0; i < N; ++i)
            {
                int count = 0;
                for (int j = -m; j <= m; ++j)
                {
                    const int idx = i + j;
                    if (idx < 0 || idx >= N)
                        continue;
                    buf[static_cast<std::size_t>(count)] = y(idx);
                    ++count;
                }
                if (count < 3)
                    continue;

                std::sort(buf.begin(), buf.begin() + count);
                const double med = (count % 2 == 1)
                    ? buf[static_cast<std::size_t>(count / 2)]
                    : 0.5 * (buf[static_cast<std::size_t>(count / 2 - 1)]
                        + buf[static_cast<std::size_t>(count / 2)]);

                int dc = 0;
                for (int t = 0; t < count; ++t)
                {
                    dev[static_cast<std::size_t>(dc)] =
                        std::abs(buf[static_cast<std::size_t>(t)] - med);
                    ++dc;
                }
                std::sort(dev.begin(), dev.begin() + dc);
                const double mad = (dc % 2 == 1)
                    ? dev[static_cast<std::size_t>(dc / 2)]
                    : 0.5 * (dev[static_cast<std::size_t>(dc / 2 - 1)]
                        + dev[static_cast<std::size_t>(dc / 2)]);

                const double sigma = 1.4826 * mad;
                const double thr = kThresh * std::max(1e-30, sigma);

                if (std::abs(y(i) - med) > thr)
                    out(i) = med;
            }

            return out;
        }

    } // namespace

    std::vector<Eigen::VectorXd> smoothPanelHistory(
        const std::vector<Eigen::VectorXd>& raw,
        double dt,
        const Options& opt)
    {
        (void)dt;

        if (raw.empty())
            return raw;

        const int N = static_cast<int>(raw.size());
        const int NE = static_cast<int>(raw.front().size());

        for (const auto& v : raw)
        {
            if (v.size() != NE)
                throw std::runtime_error("ChiHistorySmooth: inconsistent NE across time steps.");
        }

        if (opt.method == Options::Method::None)
            return raw;

        std::vector<Eigen::VectorXd> out(static_cast<std::size_t>(N),
            Eigen::VectorXd::Zero(NE));

        auto smoothSeries = [&](const Eigen::VectorXd& series) -> Eigen::VectorXd
        {
            Eigen::VectorXd y = series;

            if (opt.hampelPrewash)
                y = hampel1D(y, opt.hampelWindow, opt.hampelK);

            switch (opt.method)
            {
            case Options::Method::Gaussian:
            {
                const Eigen::VectorXd ker =
                    gaussianKernel(opt.gaussSigmaSteps, opt.gaussRadius);
                return convolveFIR(y, ker, opt.reflectBoundary);
            }
            case Options::Method::SavitzkyGolay:
            {
                const Eigen::VectorXd ker =
                    savitzkyGolayConvCoeffs(opt.sgWindow, opt.sgPolyDeg);
                return convolveFIR(y, ker, opt.reflectBoundary);
            }
            case Options::Method::TwoPassGaussianThenSg:
            {
                const Eigen::VectorXd gk =
                    gaussianKernel(opt.gaussSigmaSteps, opt.gaussRadius);
                y = convolveFIR(y, gk, opt.reflectBoundary);

                const Eigen::VectorXd sgk =
                    savitzkyGolayConvCoeffs(opt.sgWindow, opt.sgPolyDeg);
                return convolveFIR(y, sgk, opt.reflectBoundary);
            }
            case Options::Method::None:
            default:
                return y;
            }
        };

        for (int p = 0; p < NE; ++p)
        {
            Eigen::VectorXd col(N);
            for (int n = 0; n < N; ++n)
                col(n) = raw[static_cast<std::size_t>(n)](p);

            const Eigen::VectorXd colS = smoothSeries(col);

            for (int n = 0; n < N; ++n)
                out[static_cast<std::size_t>(n)](p) = colS(n);
        }

        return out;
    }

} // namespace ChiHistorySmooth
