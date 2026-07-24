#pragma once
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "../const/Const.h"

enum class EncounterBranch
{
    HeadOrBeam, // 非随浪，或 U*cos(beta) <= 0
    W1,         // 0 < w < g/(2Uc)
    W2,         // g/(2Uc) < w < g/Uc
    W3          // w > g/Uc
};

struct EncounterInfo
{
    double w_incident = 0.0;   // 入射频率
    double beta = 0.0;         // 波向
    double U = 0.0;            // 航速
    double Uc = 0.0;           // U*cos(beta)

    double split1 = 0.0;       // g / (2 U cos beta)
    double split2 = 0.0;       // g / (U cos beta)

    double we = 0.0;           // 遭遇频率
    EncounterBranch branch = EncounterBranch::HeadOrBeam;

    bool isFollowing() const { return Uc > 0.0; }
};

inline EncounterInfo calcEncounterInfo(double w, double U, double beta)
{
    EncounterInfo info;
    info.w_incident = w;
    info.beta = beta;
    info.U = U;
    info.Uc = U * std::cos(beta);

    constexpr double tol = 1e-12;

    // 顶浪 / 横浪 / 首斜浪：直接用常规公式
    if (info.Uc <= tol)
    {
        info.branch = EncounterBranch::HeadOrBeam;
        info.split1 = 0.0;
        info.split2 = 0.0;
        info.we = w - w * w * info.Uc / G;   // 此时 Uc<=0，公式仍成立
        return info;
    }

    info.split1 = G / (2.0 * info.Uc);
    info.split2 = G / info.Uc;

    // 给分界点留一点容差，避免数值正好卡在切换点
    const double eps = 1e-10 * std::max({ 1.0, std::abs(w), std::abs(info.split1), std::abs(info.split2) });

    if (std::abs(w - info.split2) <= eps)
    {
        throw std::runtime_error("incident frequency is too close to g/(U*cos(beta)); encounter frequency degenerates to zero.");
    }

    if (w < info.split1 - eps)
    {
        info.branch = EncounterBranch::W1;
        info.we = w - w * w * info.Uc / G;
    }
    else if (w < info.split2 - eps)
    {
        info.branch = EncounterBranch::W2;
        info.we = w - w * w * info.Uc / G;
    }
    else
    {
        info.branch = EncounterBranch::W3;
        info.we = -w + w * w * info.Uc / G;
    }

    if (info.we <= tol)
    {
        throw std::runtime_error("encounter frequency is non-positive.");
    }

    return info;
}