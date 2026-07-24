
#include "SeakeepingDOF.h"
#include "Element.h"
#include "../const/Const.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>

int SeakeepingDOF::findModeIndex(const SeakeepingConfig& cfg, int modeId)
{
    for (int i = 0; i < cfg.DOF; ++i)
        if (cfg.modes[i] == modeId) return i;
    return -1;
}

std::vector<std::string> SeakeepingDOF::modeNames(const SeakeepingConfig& cfg)
{
    static const char* kNames[6] = {
        "surge", "sway", "heave", "roll", "pitch", "yaw"
    };

    std::vector<std::string> names;
    names.reserve(cfg.DOF);
    for (int i = 0; i < cfg.DOF; ++i) {
        int m = cfg.modes[i];
        if (m >= 0 && m <= 5) names.push_back(kNames[m]);
        else names.push_back("mode" + std::to_string(m));
    }
    return names;
}

HydrostaticsData SeakeepingDOF::defaultHydrostatics(const SeakeepingConfig& cfg)
{
    (void)cfg;
    HydrostaticsData hs;
    hs.c33 = 628188.0;
    hs.c55 = 2.89305e+07;
    hs.c35 = 0.0;
    hs.c53 = 0.0;
    hs.c44 = 0.0;
    return hs;
}

//旧版本静水回复系数计算

HydrostaticsData SeakeepingDOF::hydrostaticsFromWaterline(
    const ShipConfig& Ship, const SeakeepingConfig& Seakeeping,
    const Element& element, double tolRel)
{
    HydrostaticsData hs{};

    if (Ship.Mass.CG.size() < 3)
        throw std::runtime_error("Hydrostatics: cfg.CG must be [xG,yG,zG].");

    const double xG = Ship.Mass.CG[0];
    const double yG = Ship.Mass.CG[1];
    const double zG = Ship.Mass.CG[2];

    const int nwl = element.n_WL;
    if (nwl <= 1) throw std::runtime_error("Hydrostatics: element.n_WL too small.");

    const double tol = tolRel * std::max(1.0, Ship.Geometry.Length);

    auto same = [&](const Eigen::Vector2d& p, const Eigen::Vector2d& q) {
        return (p - q).cwiseAbs().maxCoeff() <= tol;
    };
    auto dist = [&](const Eigen::Vector2d& p, const Eigen::Vector2d& q) {
        return (p - q).norm();
    };

    struct Seg { Eigen::Vector2d a, b; };
    std::vector<Seg> segs;
    segs.reserve(nwl);
    for (int i = 0; i < nwl; ++i) {
        segs.push_back({
            { element.xpl(i,0), element.ypl(i,0) },
            { element.xpl(i,1), element.ypl(i,1) }
            });
    }

    std::vector<char> used(segs.size(), 0);

    // 1) 璐績鈥滄帓搴忊€濇垚涓€鏉℃姌绾匡細鑳芥帴灏辨帴锛涙帴涓嶄笂灏辨ˉ鎺ュ埌鏈€杩戞绔偣
    std::vector<Eigen::Vector2d> poly;
    poly.reserve(segs.size() * 2 + 8);

    int bridges = 0;

    // 鎵句釜闈為€€鍖栫殑璧峰娈?
    int curSeg = -1;
    for (int i = 0; i < (int)segs.size(); ++i) {
        if (!same(segs[i].a, segs[i].b)) { curSeg = i; break; }
    }
    if (curSeg < 0) return hs; // 鍏ㄦ槸鐐癸紝鐩存帴杩斿洖 0

    used[curSeg] = 1;
    poly.push_back(segs[curSeg].a);
    poly.push_back(segs[curSeg].b);
    Eigen::Vector2d cur = segs[curSeg].b;

    auto anyUnused = [&]() {
        for (char u : used) if (!u) return true;
        return false;
    };

    while (anyUnused()) {
        // 鍏堟壘鑳界洿鎺ユ帴涓婄殑娈碉紙绔偣鍦?tol 鍐咃級
        int best = -1;
        bool flip = false;
        double bestScore = 1e300;

        for (int i = 0; i < (int)segs.size(); ++i) {
            if (used[i]) continue;

            if (same(segs[i].a, cur)) {
                double s = dist(segs[i].a, cur);
                if (s < bestScore) { bestScore = s; best = i; flip = false; }
            }
            if (same(segs[i].b, cur)) {
                double s = dist(segs[i].b, cur);
                if (s < bestScore) { bestScore = s; best = i; flip = true; }
            }
        }

        if (best >= 0) {
            // 鐩存帴鎺?
            used[best] = 1;
            Eigen::Vector2d nxt = flip ? segs[best].a : segs[best].b;
            if (!same(nxt, cur)) poly.push_back(nxt);
            cur = nxt;
            continue;
        }

        // 鎺ヤ笉涓婏細妗ユ帴鍒扳€滄渶杩戠殑鏈敤娈电鐐光€?
        int nearestSeg = -1;
        bool takeA = true;
        double nearestD = 1e300;

        for (int i = 0; i < (int)segs.size(); ++i) {
            if (used[i]) continue;
            double dA = dist(segs[i].a, cur);
            double dB = dist(segs[i].b, cur);
            if (dA < nearestD) { nearestD = dA; nearestSeg = i; takeA = true; }
            if (dB < nearestD) { nearestD = dB; nearestSeg = i; takeA = false; }
        }
        if (nearestSeg < 0) break;

        bridges++;

        // 鍏堝姞涓€鏉℃ˉ鎺ョ洿绾匡細cur -> endpoint
        Eigen::Vector2d endpoint = takeA ? segs[nearestSeg].a : segs[nearestSeg].b;
        if (!same(endpoint, cur)) poly.push_back(endpoint);

        // 鍐嶈蛋杩欎釜娈靛埌鍙︿竴绔?
        Eigen::Vector2d other = takeA ? segs[nearestSeg].b : segs[nearestSeg].a;
        if (!same(other, endpoint)) poly.push_back(other);

        used[nearestSeg] = 1;
        cur = other;
    }

    // 2) 寮鸿闂悎锛氭渶鍚庝竴鐐?-> 绗竴鐐?鎷夌洿绾?
    if (poly.size() < 3) return hs;
    if (!same(poly.front(), poly.back()))
        poly.push_back(poly.front());

    // 3) 鍘绘帀鐩搁偦閲嶅鐐癸紙閬垮厤 0 杈瑰奖鍝嶆暟鍊硷級
    std::vector<Eigen::Vector2d> clean;
    clean.reserve(poly.size());
    clean.push_back(poly.front());
    for (size_t i = 1; i < poly.size(); ++i) {
        if (!same(poly[i], clean.back()))
            clean.push_back(poly[i]);
    }
    if (clean.size() < 4) return hs; // 鑷冲皯寰楁湁 3 鐐?+ 闂悎鐐?

    // 4) shoelace + 浜屾鐭╋紙鎸変綘鍘熸潵鐨勶級
    double A2 = 0.0, Cx6A = 0.0, Cy6A = 0.0, Ixx12 = 0.0, Iyy12 = 0.0;
    for (int i = 0; i < (int)clean.size() - 1; ++i) {
        double x0 = clean[i].x(), y0 = clean[i].y();
        double x1 = clean[i + 1].x(), y1 = clean[i + 1].y();
        double cr = x0 * y1 - x1 * y0;

        A2 += cr;
        Cx6A += (x0 + x1) * cr;
        Cy6A += (y0 + y1) * cr;
        Ixx12 += (y0 * y0 + y0 * y1 + y1 * y1) * cr;
        Iyy12 += (x0 * x0 + x0 * x1 + x1 * x1) * cr;
    }

    if (std::abs(A2) < 1e-12) {
        std::cerr << "[Hydrostatics] Awp ~ 0 (ordering/bridge produced degenerate polygon)\n";
        return hs; // 浣犺鈥滆兘璺戔€濓紝杩欓噷杩斿洖 0
    }

    // 缁熶竴鎴愭闈㈢Н
    if (A2 < 0) { A2 = -A2; Cx6A = -Cx6A; Cy6A = -Cy6A; Ixx12 = -Ixx12; Iyy12 = -Iyy12; }

    hs.Awp = 0.5 * A2;
    hs.xF = Cx6A / (3.0 * A2);
    hs.yF = Cy6A / (3.0 * A2);

    const double Ixx0 = Ixx12 / 12.0;
    const double Iyy0 = Iyy12 / 12.0;
    hs.IT = Ixx0 - hs.Awp * hs.yF * hs.yF;
    hs.IL = Iyy0 - hs.Awp * hs.xF * hs.xF;

    // 计算浮心和排水体积
    computeBuoyancyFromSurfacePanels(element, hs, 0.0);

    // 水线面二次矩平移到 CG
    const double IT_yG = hs.IT + hs.Awp * (hs.yF - yG) * (hs.yF - yG);
    const double IL_xG = hs.IL + hs.Awp * (hs.xF - xG) * (hs.xF - xG);

    const double V = hs.Vdisp;
    const double zB = hs.zB;

    // C33
    hs.c33 = rho * G * hs.Awp;

    // C35, C53
    hs.c35 = -rho * G * hs.Awp * (hs.xF - xG);

    //hs.c35 = 46;
    hs.c53 = hs.c35;

    // C44, C55：补上浮心项
    hs.c44 = rho * G * (V * (zB - zG) + IT_yG);
    hs.c55 = rho * G * (V * (zB - zG) + IL_xG);

    std::cout << "[Hydrostatics@CG]"
        << " Awp=" << hs.Awp
        << " LCF=(" << hs.xF << "," << hs.yF << ")"
        << " B=(" << hs.xB << "," << hs.yB << "," << hs.zB << ")"
        << " V=" << hs.Vdisp
        << " xG=" << xG
        << " yG=" << yG
        << " zG=" << zG
        << " c33=" << hs.c33
        << " c44=" << hs.c44
        << " c55=" << hs.c55
        << " c35=" << hs.c35
        << std::endl;

    return hs;
}


//新修改

//HydrostaticsData SeakeepingDOF::hydrostaticsFromWaterline(
//    const ShipConfig& Ship,
//    const SeakeepingConfig& Seakeeping,
//    const Element& element,
//    double tolRel)
//{
//    (void)Seakeeping;
//
//    HydrostaticsData hs{};
//
//    if (Ship.Mass.CG.size() < 3)
//        throw std::runtime_error("Hydrostatics: cfg.CG must be [xG,yG,zG].");
//
//    const double xG = Ship.Mass.CG[0];
//    const double yG = Ship.Mass.CG[1];
//    const double zG = Ship.Mass.CG[2];
//
//    const int nwl = element.n_WL;
//    if (nwl <= 1)
//        throw std::runtime_error("Hydrostatics: element.n_WL too small.");
//
//    const double tol = tolRel * std::max(1.0, Ship.Geometry.Length);
//
//    auto samePoint = [tol](const Eigen::Vector2d& a, const Eigen::Vector2d& b) -> bool
//    {
//        return (a - b).cwiseAbs().maxCoeff() <= tol;
//    };
//
//    auto dist = [](const Eigen::Vector2d& a, const Eigen::Vector2d& b) -> double
//    {
//        return (a - b).norm();
//    };
//
//    struct Segment
//    {
//        Eigen::Vector2d a;
//        Eigen::Vector2d b;
//    };
//
//    std::vector<Segment> segs;
//    segs.reserve(static_cast<std::size_t>(nwl));
//
//    double meanSegLen = 0.0;
//    int validSegCount = 0;
//
//    for (int i = 0; i < nwl; ++i)
//    {
//        Segment s{
//            Eigen::Vector2d(element.xpl(i, 0), element.ypl(i, 0)),
//            Eigen::Vector2d(element.xpl(i, 1), element.ypl(i, 1))
//        };
//
//        if (samePoint(s.a, s.b))
//            continue;
//
//        meanSegLen += dist(s.a, s.b);
//        ++validSegCount;
//        segs.push_back(s);
//    }
//
//    if (segs.size() < 3)
//        throw std::runtime_error("Hydrostatics: too few valid waterline segments.");
//
//    meanSegLen /= std::max(1, validSegCount);
//
//    std::vector<char> used(segs.size(), 0);
//    std::vector<Eigen::Vector2d> poly;
//    poly.reserve(segs.size() + 2);
//
//    int bridges = 0;
//    double maxBridge = 0.0;
//
//    // 找第一条非退化线段作为起点
//    used[0] = 1;
//    poly.push_back(segs[0].a);
//
//    Eigen::Vector2d cur = segs[0].b;
//    poly.push_back(cur);
//
//    auto hasUnused = [&]() -> bool
//    {
//        for (char u : used)
//        {
//            if (!u)
//                return true;
//        }
//        return false;
//    };
//
//    while (hasUnused())
//    {
//        int best = -1;
//        bool flip = false;
//        double bestScore = 1.0e300;
//
//        // 1. 优先找能直接接上的线段
//        for (int i = 0; i < static_cast<int>(segs.size()); ++i)
//        {
//            if (used[static_cast<std::size_t>(i)])
//                continue;
//
//            const double da = dist(segs[static_cast<std::size_t>(i)].a, cur);
//            const double db = dist(segs[static_cast<std::size_t>(i)].b, cur);
//
//            if (samePoint(segs[static_cast<std::size_t>(i)].a, cur) && da < bestScore)
//            {
//                best = i;
//                flip = false;
//                bestScore = da;
//            }
//
//            if (samePoint(segs[static_cast<std::size_t>(i)].b, cur) && db < bestScore)
//            {
//                best = i;
//                flip = true;
//                bestScore = db;
//            }
//        }
//
//        if (best >= 0)
//        {
//            used[static_cast<std::size_t>(best)] = 1;
//
//            const Segment& s = segs[static_cast<std::size_t>(best)];
//            const Eigen::Vector2d nextPoint = flip ? s.a : s.b;
//
//            if (!samePoint(nextPoint, cur))
//                poly.push_back(nextPoint);
//
//            cur = nextPoint;
//            continue;
//        }
//
//        // 2. 如果没有能直接接上的，桥接到最近未使用线段的最近端点
//        int nearest = -1;
//        bool takeA = true;
//        double nearestDist = 1.0e300;
//
//        for (int i = 0; i < static_cast<int>(segs.size()); ++i)
//        {
//            if (used[static_cast<std::size_t>(i)])
//                continue;
//
//            const Segment& s = segs[static_cast<std::size_t>(i)];
//
//            const double da = dist(s.a, cur);
//            const double db = dist(s.b, cur);
//
//            if (da < nearestDist)
//            {
//                nearest = i;
//                takeA = true;
//                nearestDist = da;
//            }
//
//            if (db < nearestDist)
//            {
//                nearest = i;
//                takeA = false;
//                nearestDist = db;
//            }
//        }
//
//        if (nearest < 0)
//            break;
//
//        ++bridges;
//        maxBridge = std::max(maxBridge, nearestDist);
//
//        const Segment& s = segs[static_cast<std::size_t>(nearest)];
//
//        const Eigen::Vector2d endpoint = takeA ? s.a : s.b;
//        const Eigen::Vector2d other = takeA ? s.b : s.a;
//
//        // 桥接：cur -> endpoint
//        if (!samePoint(endpoint, cur))
//            poly.push_back(endpoint);
//
//        // 走当前线段：endpoint -> other
//        if (!samePoint(other, endpoint))
//            poly.push_back(other);
//
//        used[static_cast<std::size_t>(nearest)] = 1;
//        cur = other;
//    }
//
//    if (poly.size() < 3)
//        throw std::runtime_error("Hydrostatics: failed to build waterline polygon.");
//
//    // 3. 最后闭合 polygon
//    const double closeGap = dist(poly.front(), poly.back());
//
//    if (!samePoint(poly.front(), poly.back()))
//    {
//        ++bridges;
//        maxBridge = std::max(maxBridge, closeGap);
//        poly.push_back(poly.front());
//    }
//
//    // 4. 去掉相邻重复点
//    std::vector<Eigen::Vector2d> clean;
//    clean.reserve(poly.size());
//
//    clean.push_back(poly.front());
//
//    for (std::size_t i = 1; i < poly.size(); ++i)
//    {
//        if (!samePoint(poly[i], clean.back()))
//            clean.push_back(poly[i]);
//    }
//
//    if (clean.size() < 4)
//        throw std::runtime_error("Hydrostatics: cleaned waterline polygon too small.");
//
//    if (!samePoint(clean.front(), clean.back()))
//        clean.push_back(clean.front());
//
//    // 5. shoelace: Awp, LCF, waterplane second moments about origin
//    double A2 = 0.0;
//    double Cx6A = 0.0;
//    double Cy6A = 0.0;
//    double Ixx12 = 0.0;
//    double Iyy12 = 0.0;
//    double Ixy24 = 0.0;
//
//    for (std::size_t i = 0; i + 1 < clean.size(); ++i)
//    {
//        const double x0 = clean[i].x();
//        const double y0 = clean[i].y();
//        const double x1 = clean[i + 1].x();
//        const double y1 = clean[i + 1].y();
//
//        const double cr = x0 * y1 - x1 * y0;
//
//        A2 += cr;
//        Cx6A += (x0 + x1) * cr;
//        Cy6A += (y0 + y1) * cr;
//        Ixx12 += (y0 * y0 + y0 * y1 + y1 * y1) * cr;
//        Iyy12 += (x0 * x0 + x0 * x1 + x1 * x1) * cr;
//
//        // polygon product moment:
//        // Ixy_origin = 1/24 * sum((2*x0*y0 + x0*y1 + x1*y0 + 2*x1*y1) * cr)
//        Ixy24 += (2.0 * x0 * y0 + x0 * y1 + x1 * y0 + 2.0 * x1 * y1) * cr;
//    }
//
//    if (std::abs(A2) < 1.0e-12)
//        throw std::runtime_error("Hydrostatics: waterplane area is nearly zero.");
//
//    if (A2 < 0.0)
//    {
//        A2 = -A2;
//        Cx6A = -Cx6A;
//        Cy6A = -Cy6A;
//        Ixx12 = -Ixx12;
//        Iyy12 = -Iyy12;
//        Ixy24 = -Ixy24;
//    }
//
//    hs.Awp = 0.5 * A2;
//    hs.xF = Cx6A / (3.0 * A2);
//    hs.yF = Cy6A / (3.0 * A2);
//
//    const double Ixx0 = Ixx12 / 12.0;
//    const double Iyy0 = Iyy12 / 12.0;
//    const double Ixy0 = Ixy24 / 24.0;
//
//    // 水线面二次矩，先平移到 LCF
//    hs.IT = Ixx0 - hs.Awp * hs.yF * hs.yF; // ∫(y-yF)^2 dA
//    hs.IL = Iyy0 - hs.Awp * hs.xF * hs.xF; // ∫(x-xF)^2 dA
//
//    const double IxyF = Ixy0 - hs.Awp * hs.xF * hs.yF;
//
//    // 6. 体积和浮心
//    computeBuoyancyFromSurfacePanels(element, hs, 0.0);
//
//    const double V = hs.Vdisp;
//    const double zB = hs.zB;
//
//    // 7. 水线面二次矩平移到 CG
//    const double IT_CG =
//        hs.IT + hs.Awp * (hs.yF - yG) * (hs.yF - yG);
//
//    const double IL_CG =
//        hs.IL + hs.Awp * (hs.xF - xG) * (hs.xF - xG);
//
//    const double Ixy_CG =
//        IxyF + hs.Awp * (hs.xF - xG) * (hs.yF - yG);
//
//    // 8. 静水恢复系数
//    hs.c33 = rho * G * hs.Awp;
//
//    // 按你图里的公式：
//    // C35 = C53 = -rho*g*∫∫ x dA
//    // 其中 x 是相对 CG 的坐标，所以 ∫∫x dA = Awp * (xF - xG)
//    hs.c35 = -rho * G * hs.Awp * (hs.xF - xG);
//    hs.c53 = hs.c35;
//
//    // 若你的坐标系 z 向上，且图中公式采用 zB - zG，则下面保持不变。
//    hs.c44 = rho * G * (V * (zB - zG) + IT_CG);
//    hs.c55 = rho * G * (V * (zB - zG) + IL_CG);
//
//    // 如果 HydrostaticsData 里有这些成员，可以保留；
//    // 如果结构体没有 c34/c43/c45/c54，就删掉这几行。
//    hs.c34 = rho * G * hs.Awp * (hs.yF - yG);
//    hs.c43 = hs.c34;
//
//    // 对非对称水线，roll-pitch coupling 可由 ∫∫xy dA 给出。
//    // 具体正负号还要和你的 roll/pitch 正方向统一。
//    hs.c45 = rho * G * Ixy_CG;
//    hs.c54 = hs.c45;
//
//    const double maxBridgeRatio =
//        maxBridge / std::max(1.0e-12, meanSegLen);
//
//    std::cout << "[Hydrostatics@CG]"
//        << " Awp=" << hs.Awp
//        << " LCF=(" << hs.xF << "," << hs.yF << ")"
//        << " B=(" << hs.xB << "," << hs.yB << "," << hs.zB << ")"
//        << " V=" << hs.Vdisp
//        << " xG=" << xG
//        << " yG=" << yG
//        << " zG=" << zG
//        << " bridges=" << bridges
//        << " maxBridge=" << maxBridge
//        << " maxBridgeRatio=" << maxBridgeRatio
//        << " c33=" << hs.c33
//        << " c35=" << hs.c35
//        << " c44=" << hs.c44
//        << " c55=" << hs.c55
//        << std::endl;
//
//    if (bridges > 0)
//    {
//        std::cout << "[Hydrostatics] waterline was bridged:"
//            << " bridges=" << bridges
//            << " maxBridge=" << maxBridge
//            << " meanSegLen=" << meanSegLen
//            << " ratio=" << maxBridgeRatio
//            << std::endl;
//    }
//
//    return hs;
//}


void SeakeepingDOF::buildSystemMatrices(
    const ShipConfig& Ship, const SeakeepingConfig& Seakeeping,
    const HydrostaticsData& hs,
    Eigen::MatrixXd& M,
    Eigen::MatrixXd& B,
    Eigen::MatrixXd& C)
{
    const int DOF = Seakeeping.DOF;
    M.setZero(DOF, DOF);
    B.setZero(DOF, DOF);
    C.setZero(DOF, DOF);

    for (int i = 0; i < DOF; ++i) {
        const int mode = Seakeeping.modes[i];
        switch (mode) {
        case MODE_SURGE:
        case MODE_SWAY:
        case MODE_HEAVE:
            M(i, i) = Ship.Mass.Mass;
            break;
        case MODE_ROLL:
            M(i, i) = Ship.Mass.Ixx;
            break;
        case MODE_PITCH:
            M(i, i) = Ship.Mass.Iyy;
            break;
        case MODE_YAW:
            M(i, i) = Ship.Mass.Izz;
            break;
        default:
            M(i, i) = 1.0;
            break;
        }
    }

    const int isurge = findModeIndex(Seakeeping, MODE_SURGE);
    const int isway = findModeIndex(Seakeeping, MODE_SWAY);
    const int iheave = findModeIndex(Seakeeping, MODE_HEAVE);
    const int iroll = findModeIndex(Seakeeping, MODE_ROLL);
    const int ipitch = findModeIndex(Seakeeping, MODE_PITCH);
    const int iyaw = findModeIndex(Seakeeping, MODE_YAW);
    (void)isurge; (void)isway; (void)iyaw;

    if (iheave >= 0) C(iheave, iheave) = hs.c33;
    if (iroll >= 0) C(iroll, iroll) = hs.c44;
    if (ipitch >= 0) C(ipitch, ipitch) = hs.c55;

    if (iheave >= 0 && ipitch >= 0) {
        C(iheave, ipitch) = hs.c35;
        C(ipitch, iheave) = hs.c53;
    }
    if (iheave >= 0 && iroll >= 0) {
        C(iheave, iroll) = hs.c34;
        C(iroll, iheave) = hs.c43;
    }
    if (iroll >= 0 && ipitch >= 0) {
        C(iroll, ipitch) = hs.c45;
        C(ipitch, iroll) = hs.c54;
    }
}

void SeakeepingDOF::buildRadiationVn(const SeakeepingConfig& cfg,
    const std::shared_ptr<Element> element,
    const std::vector<double>& y,
    double U,
    Eigen::VectorXd& rVn)
{
    const int DOF = cfg.DOF;
    rVn.setZero(cfg.Panel.NE);

    for (int i = 0; i < DOF; ++i) {
        const int mode = cfg.modes[i];
        rVn.noalias() += y.at(i + DOF) * element->Nvec.col(mode);
    }

    const int ip = findModeIndex(cfg, MODE_PITCH);
    if (ip >= 0) {
        rVn.noalias() += U * y.at(ip) * element->Nvec.col(MODE_HEAVE);
    }

    const int iy = findModeIndex(cfg, MODE_YAW);
    if (iy >= 0) {
        rVn.noalias() -= U * y.at(iy) * element->Nvec.col(MODE_SWAY);
    }
}

void SeakeepingDOF::scaleMotionsForOutput(const SeakeepingConfig& cfg,
    const CaseContextLite& ctx,
    Eigen::MatrixXd& motionsInOut)
{
    const int DOF = cfg.DOF;
    const double rotScale = (ctx.Amp > 0.0 && ctx.W > 0.0)
        ? (G / (ctx.Amp * ctx.W * ctx.W))
        : 1.0;

    for (int i = 0; i < DOF; ++i) {
        const int mode = cfg.modes[i];
        if (mode == MODE_SURGE || mode == MODE_SWAY || mode == MODE_HEAVE) {
            if (ctx.Amp > 0.0) motionsInOut.col(i) *= (1.0 / ctx.Amp);
        }
        else if (mode == MODE_ROLL || mode == MODE_PITCH || mode == MODE_YAW) {
            motionsInOut.col(i) *= rotScale;
        }
    }
}

double SeakeepingDOF::rollDampingFromZeta(const ShipConfig& cfg,
    const HydrostaticsData& hs,
    double zeta)
{
    if (hs.c44 <= 0.0 || cfg.Mass.Ixx <= 0.0) return 0.0;
    return 2.0 * zeta * std::sqrt(hs.c44 * cfg.Mass.Ixx);
}


void  SeakeepingDOF::DumpWaterlineDebugCSV(
    const ShipConfig& Ship, const SeakeepingConfig& Seakeeping,
    const Element& element,
    const std::string& prefix,
    double tolRel)
{
    const int nwl0 = element.n_WL;
    if (nwl0 <= 0) throw std::runtime_error("DumpWaterlineDebugCSV: element.n_WL <= 0");

    const double tol = tolRel * std::max(1.0, Ship.Geometry.Length);

    // ---- local key structs (still inside ONE function) ----
    struct Key {
        long long ix, iy;
        bool operator==(const Key& o) const noexcept { return ix == o.ix && iy == o.iy; }
    };
    struct KeyHash {
        size_t operator()(const Key& k) const noexcept {
            // simple mix
            auto h1 = std::hash<long long>()(k.ix * 1315423911LL);
            auto h2 = std::hash<long long>()(k.iy + 0x9e3779b97f4a7c15ULL);
            return h1 ^ (h2 + (h1 << 6) + (h1 >> 2));
        }
    };
    struct EKey {
        Key a, b; // unordered edge (sorted)
        bool operator==(const EKey& o) const noexcept { return a == o.a && b == o.b; }
    };
    struct EKeyHash {
        size_t operator()(const EKey& e) const noexcept {
            size_t h1 = KeyHash{}(e.a);
            size_t h2 = KeyHash{}(e.b);
            return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
        }
    };

    auto qkey = [&](double x, double y) -> Key {
        return Key{ (long long)std::llround(x / tol), (long long)std::llround(y / tol) };
    };
    auto make_ekey = [&](Key u, Key v) -> EKey {
        if (u.ix > v.ix || (u.ix == v.ix && u.iy > v.iy)) std::swap(u, v);
        return EKey{ u, v };
    };

    struct Seg { double ax, ay, bx, by; int mirrored; };

    // 1) gather segments
    std::vector<Seg> segs;
    segs.reserve(nwl0 * 2);

    double miny = 1e100, maxy = -1e100;
    for (int i = 0; i < nwl0; ++i) {
        double ax = element.xpl(i, 0), ay = element.ypl(i, 0);
        double bx = element.xpl(i, 1), by = element.ypl(i, 1);
        segs.push_back({ ax, ay, bx, by, 0 });
        miny = std::min({ miny, ay, by });
        maxy = std::max({ maxy, ay, by });
    }

    // 2) mirror for halfship (ALWAYS mirror if halfship; skip centerline duplicates)
    if (Seakeeping.Panel.NEType == "halfship") {
        const int old = (int)segs.size();
        segs.reserve(2 * old);
        for (int i = 0; i < old; ++i) {
            auto s = segs[i];
            if (std::abs(s.ay) < tol && std::abs(s.by) < tol) continue; // centerline
            s.ay = -s.ay; s.by = -s.by;
            s.mirrored = 1;
            segs.push_back(s);
        }
    }

    // 3) count degrees and edge multiplicity
    std::unordered_map<Key, int, KeyHash> deg;
    std::unordered_map<EKey, int, EKeyHash> ecount;
    deg.reserve(segs.size() * 2);
    ecount.reserve(segs.size() * 2);

    for (const auto& s : segs) {
        Key ka = qkey(s.ax, s.ay);
        Key kb = qkey(s.bx, s.by);
        deg[ka]++; deg[kb]++;
        ecount[make_ekey(ka, kb)]++;
    }

    int n1 = 0, n2 = 0, n3p = 0;
    for (const auto& kv : deg) {
        if (kv.second == 1) n1++;
        else if (kv.second == 2) n2++;
        else n3p++;
    }

    std::cout << "[WL Debug] nwl0=" << nwl0
        << " segs(after mirror)=" << segs.size()
        << " tol=" << tol
        << " deg1=" << n1 << " deg2=" << n2 << " deg>=3=" << n3p
        << " miny=" << miny << " maxy=" << maxy
        << "\n";

    // 4) write segs_all
    {
        std::ofstream f(prefix + "_segs_all.csv");
        f << "id,ax,ay,bx,by,mirrored,edge_count\n";
        for (int i = 0; i < (int)segs.size(); ++i) {
            const auto& s = segs[i];
            Key ka = qkey(s.ax, s.ay);
            Key kb = qkey(s.bx, s.by);
            int cnt = ecount[make_ekey(ka, kb)];
            f << i << "," << s.ax << "," << s.ay << "," << s.bx << "," << s.by
                << "," << s.mirrored << "," << cnt << "\n";
        }
    }

    // 5) write segs_boundary (edge_count==1)
    {
        std::ofstream f(prefix + "_segs_boundary.csv");
        f << "id,ax,ay,bx,by,mirrored\n";
        int outId = 0;
        for (int i = 0; i < (int)segs.size(); ++i) {
            const auto& s = segs[i];
            Key ka = qkey(s.ax, s.ay);
            Key kb = qkey(s.bx, s.by);
            int cnt = ecount[make_ekey(ka, kb)];
            if (cnt != 1) continue;
            f << outId++ << "," << s.ax << "," << s.ay << "," << s.bx << "," << s.by
                << "," << s.mirrored << "\n";
        }
    }

    // 6) write nodes_bad (deg!=2)
    {
        std::ofstream f(prefix + "_nodes_bad.csv");
        f << "x,y,deg\n";
        for (const auto& kv : deg) {
            if (kv.second == 2) continue;
            // convert quantized key back to coord grid (enough for plotting)
            double x = kv.first.ix * tol;
            double y = kv.first.iy * tol;
            f << x << "," << y << "," << kv.second << "\n";
        }
    }

    std::cout << "  wrote:\n"
        << "   " << prefix << "_segs_all.csv\n"
        << "   " << prefix << "_segs_boundary.csv\n"
        << "   " << prefix << "_nodes_bad.csv\n";
}



void SeakeepingDOF::computeBuoyancyFromSurfacePanels(
    const Element& element,
    HydrostaticsData& hs,
    double zWL)
{
    const Eigen::Vector3d ref(0.0, 0.0, zWL);

    double sumV6 = 0.0;
    Eigen::Vector3d sumMoment = Eigen::Vector3d::Zero();

    auto addTriangle = [&](const Eigen::Vector3d& P0,
        const Eigen::Vector3d& P1,
        const Eigen::Vector3d& P2)
    {
        const Eigen::Vector3d a = P0 - ref;
        const Eigen::Vector3d b = P1 - ref;
        const Eigen::Vector3d c = P2 - ref;

        const double area2 = (b - a).cross(c - a).norm();
        if (area2 < 1.0e-14) return;

        const double v6 = a.dot(b.cross(c));
        if (!std::isfinite(v6)) return;

        sumV6 += v6;

        // 三角面片和参考点构成的虚拟体积单元质心：
        // C = ref + (a + b + c) / 4
        // 所以一阶矩为 V * C。
        // 这里先用 6V 累加，最后统一除。
        sumMoment += v6 * (a + b + c);
    };

    for (int e = 0; e < element.NE; ++e)
    {
        const ElementMatrix& p = element.ElementData->at(e);

        Eigen::Vector3d P0 = p.row(0).transpose();
        Eigen::Vector3d P1 = p.row(1).transpose();
        Eigen::Vector3d P2 = p.row(2).transpose();
        Eigen::Vector3d P3 = p.row(3).transpose();

        // 四边形面元拆成两个三角形。
        // 这不是生成三角网格，只是做体积积分。
        addTriangle(P0, P1, P2);
        addTriangle(P0, P2, P3);
    }

    if (std::abs(sumV6) < 1.0e-12)
    {
        std::cerr << "[BuoyancySurfacePanels] failed: volume is nearly zero. "
            << "Check panel orientation, waterline z, or mesh scale.\n";

        hs.Vdisp = 0.0;
        hs.xB = 0.0;
        hs.yB = 0.0;
        hs.zB = 0.0;
        return;
    }

    const double signedV = sumV6 / 6.0;

    // B = ref + sum[v6 * (a+b+c)] / (4 * sum[v6])
    Eigen::Vector3d B = ref + sumMoment / (4.0 * sumV6);

    hs.Vdisp = std::abs(signedV);
    hs.xB = B.x();
    hs.yB = B.y();
    hs.zB = B.z();

    std::cout << "[BuoyancySurfacePanels]"
        << " V=" << hs.Vdisp
        << " B=(" << hs.xB << ", " << hs.yB << ", " << hs.zB << ")"
        << std::endl;
}

