#pragma once

#include <array>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace wigleyi
{
    struct WigleyPoint3D
    {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
    };

    using WigleyQuad = std::array<WigleyPoint3D, 4>;

    inline double WigleyIHalfBreadth(double x, double z, double L, double B, double T)
    {
        if (L <= 0.0 || B <= 0.0 || T <= 0.0)
            throw std::invalid_argument("WigleyIHalfBreadth: L, B, T must be positive.");

        const double xt = 2.0 * x / L;
        const double zt = z / T;

        const double y1 = 1.0 - zt * zt;
        const double y2 = 1.0 - xt * xt;
        const double y3 = 1.0 + 0.2 * xt * xt;
        const double y4 = zt * zt;
        const double y5 = 1.0 - std::pow(zt, 8.0);
        const double y6 = std::pow(y2, 4.0);

        const double y =
            0.5 * B * (y1 * y2 * y3 + y4 * y5 * y6);

        return std::max(0.0, y);
    }

    inline WigleyPoint3D EvalWigleyIPoint(double x, double z, double L, double B, double T)
    {
        return WigleyPoint3D{ x, WigleyIHalfBreadth(x, z, L, B, T), z };
    }

    inline std::vector<std::vector<WigleyPoint3D>>
        BuildWigleyIGridStarboard(double L, double B, double T, int NS, int NP)
    {
        if (NS < 2 || NP < 2)
            throw std::invalid_argument("BuildWigleyIGridStarboard: NS and NP must be >= 2.");

        std::vector<std::vector<WigleyPoint3D>> grid(
            static_cast<std::size_t>(NS),
            std::vector<WigleyPoint3D>(static_cast<std::size_t>(NP)));

        const double dx = L / static_cast<double>(NS - 1);
        const double dz = T / static_cast<double>(NP - 1);

        for (int i = 0; i < NS; ++i)
        {
            const double x = -0.5 * L + i * dx;

            for (int j = 0; j < NP; ++j)
            {
                const double z = -(NP - 1 - j) * dz; // from -T to 0
                grid[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                    EvalWigleyIPoint(x, z, L, B, T);
            }
        }

        return grid;
    }

    inline WigleyQuad MirrorQuadToPort(const WigleyQuad& q)
    {
        WigleyPoint3D p1{ q[0].x, -q[0].y, q[0].z };
        WigleyPoint3D p2{ q[1].x, -q[1].y, q[1].z };
        WigleyPoint3D p3{ q[2].x, -q[2].y, q[2].z };
        WigleyPoint3D p4{ q[3].x, -q[3].y, q[3].z };

        // reverse winding to keep outward normal consistent
        return WigleyQuad{ p1, p4, p3, p2 };
    }

    inline void AppendRefinedPatchOnHull(
        std::vector<WigleyQuad>& quads,
        double x0, double x1,
        double z0, double z1,
        double L, double B, double T,
        int refineS,
        int refineP)
    {
        if (refineS < 1 || refineP < 1)
            throw std::invalid_argument("AppendRefinedPatchOnHull: refineS and refineP must be >= 1.");

        const double dx = (x1 - x0) / static_cast<double>(refineS);
        const double dz = (z1 - z0) / static_cast<double>(refineP);

        for (int ii = 0; ii < refineS; ++ii)
        {
            const double xa = x0 + ii * dx;
            const double xb = x0 + (ii + 1) * dx;

            for (int jj = 0; jj < refineP; ++jj)
            {
                const double za = z0 + jj * dz;
                const double zb = z0 + (jj + 1) * dz;

                // 保持与原始面元一致的绕序:
                // p1(左下) -> p2(右下) -> p3(右上) -> p4(左上)
                const WigleyPoint3D p1 = EvalWigleyIPoint(xa, za, L, B, T);
                const WigleyPoint3D p2 = EvalWigleyIPoint(xb, za, L, B, T);
                const WigleyPoint3D p3 = EvalWigleyIPoint(xb, zb, L, B, T);
                const WigleyPoint3D p4 = EvalWigleyIPoint(xa, zb, L, B, T);

                quads.push_back(WigleyQuad{ p1, p2, p3, p4 });
            }
        }
    }

    //========================================================
    //  新版本：显式给出水线附近面元细化参数
    //
    //  NS, NP   : 整体母网格节点数
    //  WL_NS    : 水线附近每个原始面元沿 x/船长方向切成几份
    //  WL_NP    : 水线附近每个原始面元沿 z/垂向切成几份
    //
    //  这里只对最上面一圈面元（j == NP-2）做局部细化
    //========================================================
    inline std::size_t CountWigleyISeakeepingPanels(
        int NS, int NP, int WL_NS, int WL_NP, bool fullShip)
    {
        if (NS < 2 || NP < 2)
            throw std::invalid_argument("CountWigleyISeakeepingPanels: NS and NP must be >= 2.");
        if (WL_NS < 1 || WL_NP < 1)
            throw std::invalid_argument("CountWigleyISeakeepingPanels: WL_NS and WL_NP must be >= 1.");

        // 除最上面一圈外，其余面元数
        const std::size_t nOther =
            static_cast<std::size_t>((NS - 1) * (NP - 2));

        // 最上面一圈原本有 (NS-1) 个面元
        // 每个变成 WL_NS * WL_NP 个
        const std::size_t nTop =
            static_cast<std::size_t>((NS - 1) * WL_NS * WL_NP);

        const std::size_t nStarboard = nOther + nTop;
        return fullShip ? 2 * nStarboard : nStarboard;
    }

    inline std::vector<WigleyQuad>
        BuildWigleyIQuads(
            double L, double B, double T,
            int NS, int NP,
            int WL_NS, int WL_NP,
            bool fullShip)
    {
        if (NS < 2 || NP < 2)
            throw std::invalid_argument("BuildWigleyIQuads: NS and NP must be >= 2.");
        if (WL_NS < 1 || WL_NP < 1)
            throw std::invalid_argument("BuildWigleyIQuads: WL_NS and WL_NP must be >= 1.");

        auto grid = BuildWigleyIGridStarboard(L, B, T, NS, NP);

        std::vector<WigleyQuad> quads;
        quads.reserve(CountWigleyISeakeepingPanels(NS, NP, WL_NS, WL_NP, fullShip));

        // -------------------------
        // starboard
        // -------------------------
        for (int i = 0; i < NS - 1; ++i)
        {
            for (int j = 0; j < NP - 1; ++j)
            {
                const auto& p1 = grid[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
                const auto& p2 = grid[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(j)];
                const auto& p3 = grid[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(j + 1)];
                const auto& p4 = grid[static_cast<std::size_t>(i)][static_cast<std::size_t>(j + 1)];

                // 最上面一圈：连接倒数第二层到水线 z=0 的这一圈
                if (j == NP - 2)
                {
                    AppendRefinedPatchOnHull(
                        quads,
                        p1.x, p2.x,
                        p1.z, p4.z,
                        L, B, T,
                        WL_NS, WL_NP);
                }
                else
                {
                    quads.push_back(WigleyQuad{ p1, p2, p3, p4 });
                }
            }
        }

        // -------------------------
        // port
        // -------------------------
        if (fullShip)
        {
            const std::size_t nStarboard = quads.size();
            for (std::size_t k = 0; k < nStarboard; ++k)
            {
                quads.push_back(MirrorQuadToPort(quads[k]));
            }
        }

        return quads;
    }

    inline void WritePanelsToFile(
        const std::vector<WigleyQuad>& quads,
        const std::string& fileName,
        int precision = 3)
    {
        std::ofstream out(fileName);
        if (!out.is_open())
            throw std::runtime_error("WritePanelsToFile: cannot open output file: " + fileName);

        out.setf(std::ios::scientific);
        out << std::setprecision(precision);

        for (std::size_t i = 0; i < quads.size(); ++i)
        {
            out << (i + 1) << "\n";

            for (int k = 0; k < 4; ++k)
            {
                out << std::showpos
                    << std::setw(14) << quads[i][k].x
                    << std::setw(16) << quads[i][k].y
                    << std::setw(16) << quads[i][k].z
                    << "\n";
            }

            out << "\n";
        }
    }

    inline void WriteWigleyISeakeepingPanels(
        double L,
        double B,
        double T,
        const std::string& fileName,
        int NS,
        int NP,
        int WL_NS,
        int WL_NP,
        bool fullShip = true,
        int precision = 3)
    {
        const auto quads = BuildWigleyIQuads(L, B, T, NS, NP, WL_NS, WL_NP, fullShip);
        WritePanelsToFile(quads, fileName, precision);
    }

    //========================================================
    //  兼容旧接口
    //  保持你当前这份代码的默认行为：
    //  水线附近一圈默认细化为 2 x 2
    //========================================================
    inline std::size_t CountWigleyISeakeepingPanels(
        int NS, int NP, bool fullShip = true)
    {
        return CountWigleyISeakeepingPanels(NS, NP, 2, 2, fullShip);
    }

    inline std::vector<WigleyQuad>
        BuildWigleyIQuads(
            double L, double B, double T,
            int NS = 22, int NP = 6,
            bool fullShip = true)
    {
        return BuildWigleyIQuads(L, B, T, NS, NP, 2, 2, fullShip);
    }

    inline void WriteWigleyISeakeepingPanels(
        double L,
        double B,
        double T,
        const std::string& fileName,
        int NS = 22,
        int NP = 6,
        bool fullShip = true,
        int precision = 3)
    {
        WriteWigleyISeakeepingPanels(L, B, T, fileName, NS, NP, 2, 2, fullShip, precision);
    }
}