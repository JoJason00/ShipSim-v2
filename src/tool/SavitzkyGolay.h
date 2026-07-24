// SavitzkyGolay.h — header-only, Eigen3
// SG 平滑 + SG 一阶时间导数（与 Python scipy.signal.savgol_filter 同构：
//   先对 y 做 deriv=0，再对平滑结果做 deriv=1，导数除以 dt_phys）
//
// 用法：
//   #include "SavitzkyGolay.h"
//   auto r = sg::smoothThenDeriv1(y, dt_const, 31, 41, 3);
//   // r.ySmooth -> 替代 smooth5Binomial 输出；r.dydt -> 替代 differentiateCentral
//
#pragma once

#include <Eigen/Core>
#include <stdexcept>

namespace sg
{

    inline int makeOdd(int w)
    {
        w = std::max(3, w);
        if ((w % 2) == 0) ++w;
        return w;
    }

    // 卷积系数 c，长度 windowLength；与 y 做点积得到窗口中心处的平滑值(deriv=0)
    // 或样本索引方向的一阶导数系数(deriv=1)，再除以物理步长 dt 得 dy/dt。
    inline Eigen::RowVectorXd convolutionCoeffs(int windowLength, int polyOrder, int derivOrder)
    {
        windowLength = makeOdd(windowLength);
        if (polyOrder < 1)
            throw std::invalid_argument("sg::convolutionCoeffs: polyOrder >= 1 required");
        if (derivOrder < 0 || derivOrder > polyOrder)
            throw std::invalid_argument("sg::convolutionCoeffs: 0 <= derivOrder <= polyOrder");
        if (windowLength < polyOrder + 2)
            throw std::invalid_argument("sg::convolutionCoeffs: windowLength too small for polyOrder");

        const int half = windowLength / 2;
        Eigen::MatrixXd A(windowLength, polyOrder + 1);
        for (int i = 0; i < windowLength; ++i)
        {
            const double m = static_cast<double>(i - half);
            double mp = 1.0;
            for (int j = 0; j <= polyOrder; ++j)
            {
                A(i, j) = mp;
                mp *= m;
            }
        }

        const Eigen::MatrixXd AtA = A.transpose() * A;
        const Eigen::MatrixXd ATAinvAT = AtA.ldlt().solve(A.transpose());
        return ATAinvAT.row(derivOrder);
    }

    // 边界：两端复制延拓后卷积，输出与 y 同长度（与 scipy mode="nearest" 类似）
    inline Eigen::VectorXd convolveSame(
        const Eigen::VectorXd& y,
        const Eigen::RowVectorXd& c)
    {
        const int n = static_cast<int>(y.size());
        const int w = static_cast<int>(c.size());
        const int half = w / 2;

        if (n == 0)
            return y;

        Eigen::VectorXd pad(n + 2 * half);
        for (int i = 0; i < half; ++i)
            pad(i) = y(0);
        pad.segment(half, n) = y;
        for (int i = 0; i < half; ++i)
            pad(half + n + i) = y(n - 1);

        Eigen::VectorXd out(n);
        for (int i = 0; i < n; ++i)
        {
            double s = 0.0;
            for (int k = 0; k < w; ++k)
                s += c(k) * pad(i + k);
            out(i) = s;
        }
        return out;
    }

    struct SmoothDeriv1Result
    {
        Eigen::VectorXd ySmooth;
        Eigen::VectorXd dydt;
    };

    // 先平滑再求导；dt 须与 y 的采样间隔一致（你代码里的 ctx.dt_const，单位秒）
    inline SmoothDeriv1Result smoothThenDeriv1(
        const Eigen::VectorXd& y,
        double dt,
        int winSmooth,
        int winDeriv,
        int polyOrder = 3)
    {
        if (!(dt > 0.0) || !std::isfinite(dt))
            throw std::invalid_argument("sg::smoothThenDeriv1: dt must be finite and > 0");

        winSmooth = makeOdd(winSmooth);
        winDeriv = makeOdd(winDeriv);

        const Eigen::RowVectorXd c0 = convolutionCoeffs(winSmooth, polyOrder, 0);
        const Eigen::VectorXd yS = convolveSame(y, c0);

        const Eigen::RowVectorXd c1 = convolutionCoeffs(winDeriv, polyOrder, 1);
        const Eigen::VectorXd dIdx = convolveSame(yS, c1);

        SmoothDeriv1Result r;
        r.ySmooth = yS;
        r.dydt = dIdx / dt;
        return r;
    }

} // namespace sg
