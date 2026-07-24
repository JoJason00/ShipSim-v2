#include "RAO4DTable.h"
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <iomanip>
#include <unordered_map>
#include <cstring>   // std::strtod
#include <cctype>
#include <iostream>
#include "../const/Const.h"

#if __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#endif

static inline double wrap360(double deg) {
    deg = std::fmod(deg, 360.0);
    if (deg < 0) deg += 360.0;
    return deg;
}
static inline double rad2deg(double x) {
    return x * 180.0 / 3.141592653589793238462643383279502884;
}
static inline bool is_nan_c(const std::complex<double>& z) {
    return std::isnan(z.real()) || std::isnan(z.imag());
}
static inline long long qkey(double x, double scale) {
    return (long long)std::llround(x * scale);
}

// 快速读取一个 CSV field 为 double（支持 nan / 科学计数法）；p 会推进到下一个 field
static inline double parseDoubleField(const char*& p) {
    while (*p == ' ' || *p == '\t') ++p;

    // 空字段
    if (*p == ',' || *p == '\0' || *p == '\r' || *p == '\n') {
        if (*p == ',') ++p;
        return std::numeric_limits<double>::quiet_NaN();
    }

    char* end = nullptr;
    double v = std::strtod(p, &end);

    // 没读到数字：跳过到下一个逗号
    if (end == p) {
        while (*p && *p != ',') ++p;
        if (*p == ',') ++p;
        return std::numeric_limits<double>::quiet_NaN();
    }

    p = end;
    // 跳到下一个字段
    while (*p == ' ' || *p == '\t') ++p;
    if (*p == ',') ++p;
    return v;
}

// 轻量 split 用于 meta/axes 行（不需要很快）
static inline std::vector<std::string> splitCSV(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    cur.reserve(line.size());
    for (char c : line) {
        if (c == ',') { out.push_back(cur); cur.clear(); }
        else if (c != '\r') cur.push_back(c);
    }
    out.push_back(cur);
    return out;
}

static inline void sort_unique(std::vector<double>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

bool RAO4DTable::fileExists(const std::string& path) {
#if __has_include(<filesystem>)
    std::error_code ec;
    return fs::exists(path, ec);
#else
    std::ifstream f(path);
    return f.good();
#endif
}

void RAO4DTable::finalizeAndAllocate(bool fillNaN)
{
    nFn = (int)Fns.size();
    nBe = (int)dir.size();
    nOm = (int)w.size();
    nDof = (int)modeIds.size();
    if (nFn <= 0 || nBe <= 0 || nOm <= 0 || nDof <= 0)
        throw std::runtime_error("RAO4D dims invalid");

    strideOm = (size_t)nDof;
    strideBe = (size_t)nOm * (size_t)nDof;
    strideFn = (size_t)nBe * (size_t)nOm * (size_t)nDof;

    const size_t total = (size_t)nFn * strideFn;

    if (fillNaN) {
        double qnan = std::numeric_limits<double>::quiet_NaN();
        H.assign(total, { qnan, qnan });
    }
    else {
        H.assign(total, { 0.0, 0.0 });
    }
}

void RAO4DTable::setAt(int iFn_, int iBe_, int iOm_, const std::vector<std::complex<double>>& Hk)
{
    if ((int)Hk.size() != nDof) throw std::runtime_error("Hk size mismatch");
    for (int k = 0; k < nDof; ++k) at(iFn_, iBe_, iOm_, k) = Hk[k];
}

static inline void bracket(const std::vector<double>& xs, double x, int& i0, int& i1, double& w) {
    if (xs.size() < 2) throw std::runtime_error("axis too small");
    if (x <= xs.front()) { i0 = 0; i1 = 1; w = 0.0; return; }
    if (x >= xs.back()) { i0 = (int)xs.size() - 2; i1 = (int)xs.size() - 1; w = 1.0; return; }

    auto it = std::upper_bound(xs.begin(), xs.end(), x);
    i1 = (int)std::distance(xs.begin(), it);
    i0 = i1 - 1;
    w = (x - xs[i0]) / (xs[i1] - xs[i0]);
}

static inline void bracket_beta_periodic(const std::vector<double>& dirDegAxis, double betaDeg,
    int& i0, int& i1, double& w) {
    double x = wrap360(betaDeg);
    if (x < dirDegAxis.front()) x += 360.0;

    if (x > dirDegAxis.back()) {
        i0 = (int)dirDegAxis.size() - 1;
        i1 = 0;
        double x0 = dirDegAxis.back();
        double x1 = dirDegAxis.front() + 360.0;
        w = (x - x0) / (x1 - x0);
        return;
    }
    bracket(dirDegAxis, x, i0, i1, w);
}

std::vector<std::complex<double>> RAO4DTable::interp(double Fn, double betaDeg, double omega) const
{
    int f0, f1, b0, b1, o0, o1;
    double wf, wb, wo;

    bracket(Fns, Fn, f0, f1, wf);
    bracket_beta_periodic(dir, betaDeg, b0, b1, wb);
    bracket(w, omega, o0, o1, wo);

    std::vector<std::complex<double>> out(nDof);

    for (int k = 0; k < nDof; ++k) {
        auto V000 = at(f0, b0, o0, k), V100 = at(f1, b0, o0, k);
        auto V010 = at(f0, b1, o0, k), V110 = at(f1, b1, o0, k);
        auto V001 = at(f0, b0, o1, k), V101 = at(f1, b0, o1, k);
        auto V011 = at(f0, b1, o1, k), V111 = at(f1, b1, o1, k);

        auto V00 = V000 * (1 - wf) + V100 * wf;
        auto V10 = V010 * (1 - wf) + V110 * wf;
        auto V01 = V001 * (1 - wf) + V101 * wf;
        auto V11 = V011 * (1 - wf) + V111 * wf;

        auto V0 = V00 * (1 - wb) + V10 * wb;
        auto V1 = V01 * (1 - wb) + V11 * wb;

        out[k] = V0 * (1 - wo) + V1 * wo;
    }
    return out;
}

size_t RAO4DTable::countNaN() const
{
    size_t c = 0;
    for (const auto& z : H) if (is_nan_c(z)) ++c;
    return c;
}

// -------- CSV write/read --------
//
// 文件结构：
// meta/axes 若干行（key,...）
// 空行
// data header: Fn,beta_deg,omega, modeX_re,modeX_im,(modeX_abs,modeX_phase_deg)...
// data rows

void RAO4DTable::writeCSV(const std::string& path, std::vector<double>& non_we, int precision, bool includeAbsPhase) const
{
    if (nFn <= 0 || nBe <= 0 || nOm <= 0 || nDof <= 0 || H.empty())
        throw std::runtime_error("writeCSV: table not initialized");

    std::ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("writeCSV: cannot open file");

    out.setf(std::ios::scientific);
    out << std::setprecision(precision);

    double L = meta.L;
    // meta
    out << "RAO4DCSV,version,1\n";
    out << "dirUnit,Deg\n";                 // 写出统一用 Deg
    out << "freqUnit,RadPerSec\n";          // 写出统一用 RadPerSec
    out << "freqType," << (meta.freqType == RAO4DMeta::FreqType::IncidentOmega ? "IncidentOmega" : "EncounterOmega") << "\n";
    out << "normType," << (meta.normType == RAO4DMeta::NormType::Raw ? "Raw" : "TransA_RotkA") << "\n";
    out << "L," << L << "\n";

    // axes
    out << "modeIds";
    for (int id : modeIds) out << "," << id;
    out << "\n";

    out << "Fns";
    for (double v : Fns) out << "," << v;
    out << "\n";

    out << "dir_deg";
    // 表内若不是 Deg（旧数据），这里也转成 Deg 写出
    if (meta.dirUnit == RAO4DMeta::DirUnit::Deg) {
        for (double v : dir) out << "," << wrap360(v);
    }
    else {
        for (double v : dir) out << "," << wrap360(rad2deg(v));
    }
    out << "\n";

    out << "omega";
    for (double v : w) out << "," << v;
    out << "\n";

    out << "\n";

    // header (Wide)
    out << "Fn,beta_deg,we, non_we,non_nada";
    for (int k = 0; k < nDof; ++k) {
        out << ",mode" << modeIds[k] << "_re"
            << ",mode" << modeIds[k] << "_im";
        if (includeAbsPhase) {
            out << ",mode" << modeIds[k] << "_abs"
                << ",mode" << modeIds[k] << "_phase_deg";
        }
    }
    out << "\n";

    auto betaDegAt = [&](int iBe) {
        double b = dir[iBe];
        if (meta.dirUnit == RAO4DMeta::DirUnit::Rad) b = rad2deg(b);
        return wrap360(b);
    };

    // rows

    int i_we = 0;
    double non_nada;
    double wsquare = 2.0 * PI * G / L;

    double non_omiga = sqrt(G / L);

    for (int iFn = 0; iFn < nFn; ++iFn) {
        for (int iBe = 0; iBe < nBe; ++iBe) {
            const double betaDeg = betaDegAt(iBe);
            for (int iOm = 0; iOm < nOm; ++iOm) {
                non_nada = wsquare / (w[iOm] * w[iOm]);
                out << Fns[iFn] << "," << betaDeg << "," << non_we[i_we] * non_omiga << "," << non_we[i_we] << "," << non_nada;
                i_we++;
                for (int k = 0; k < nDof; ++k) {
                    const auto z = at(iFn, iBe, iOm, k);
                    if (is_nan_c(z)) {
                        out << ",nan,nan";
                        if (includeAbsPhase) out << ",nan,nan";
                    }
                    else {
                        const double re = z.real();
                        const double im = z.imag();
                        out << "," << re << "," << im;
                        if (includeAbsPhase) {
                            const double ab = std::hypot(re, im);
                            const double ph = std::atan2(im, re); // rad
                            out << "," << ab << "," << wrap360(rad2deg(ph));
                        }
                    }
                }
                out << "\n";
            }
        }
    }

    std::cout << "\nwrite RAO done.\n\n";
}

RAO4DTable RAO4DTable::readCSV(const std::string& path)
{
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("readCSV: cannot open file: " + path);

    RAO4DTable t;

    std::string line;

    // ---- meta/axes ----
    std::vector<int>    modeIds;
    std::vector<double> Fns, dirDeg, omega;
    bool sawVersion = false;

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) break;

        auto cols = splitCSV(line);
        if (cols.empty()) continue;

        const std::string& key = cols[0];

        if (key == "RAO4DCSV") {
            // RAO4DCSV,version,1
            if (cols.size() < 3 || cols[1] != "version" || cols[2] != "1")
                throw std::runtime_error("readCSV: unsupported version line");
            sawVersion = true;
        }
        else if (key == "freqType") {
            if (cols.size() < 2) throw std::runtime_error("readCSV: bad freqType");
            t.meta.freqType = (cols[1] == "EncounterOmega")
                ? RAO4DMeta::FreqType::EncounterOmega
                : RAO4DMeta::FreqType::IncidentOmega;
        }
        else if (key == "normType") {
            if (cols.size() < 2) throw std::runtime_error("readCSV: bad normType");
            t.meta.normType = (cols[1] == "Raw")
                ? RAO4DMeta::NormType::Raw
                : RAO4DMeta::NormType::TransA_RotkA;
        }
        else if (key == "L") {
            if (cols.size() < 2) throw std::runtime_error("readCSV: bad L line");
            t.meta.L = std::stod(cols[1]);
        }
        else if (key == "modeIds") {
            modeIds.clear();
            for (size_t i = 1; i < cols.size(); ++i)
                if (!cols[i].empty()) modeIds.push_back(std::stoi(cols[i]));
        }
        else if (key == "Fns") {
            Fns.clear();
            for (size_t i = 1; i < cols.size(); ++i)
                if (!cols[i].empty()) Fns.push_back(std::stod(cols[i]));
        }
        else if (key == "dir_deg") {
            dirDeg.clear();
            for (size_t i = 1; i < cols.size(); ++i)
                if (!cols[i].empty()) dirDeg.push_back(wrap360(std::stod(cols[i])));
        }
        else if (key == "omega") {
            omega.clear();
            for (size_t i = 1; i < cols.size(); ++i)
                if (!cols[i].empty()) omega.push_back(std::stod(cols[i]));
        }
        else {
            // 其他 key 忽略（保留扩展空间）
        }
    }

    if (!sawVersion) throw std::runtime_error("readCSV: missing RAO4DCSV version header");
    if (modeIds.empty() || Fns.empty() || dirDeg.empty() || omega.empty())
        throw std::runtime_error("readCSV: missing axes (modeIds/Fns/dir_deg/omega)");

    // table 内统一用 Deg & rad/s
    t.meta.dirUnit = RAO4DMeta::DirUnit::Deg;
    t.meta.freqUnit = RAO4DMeta::FreqUnit::RadPerSec;

    t.modeIds = std::move(modeIds);
    t.Fns = std::move(Fns);
    t.dir = std::move(dirDeg);
    t.w = std::move(omega);

    sort_unique(t.Fns);
    sort_unique(t.dir);
    sort_unique(t.w);

    t.finalizeAndAllocate(true);

    // ---- data header ----
    if (!std::getline(in, line))
        throw std::runtime_error("readCSV: missing data header");
    if (!line.empty() && line.back() == '\r') line.pop_back();

    // 用表头列数判断：每个 mode 是 2 列(re/im)还是 4 列(re/im/abs/phase_deg)
    auto headerCols = splitCSV(line);
    if (headerCols.size() < 3)
        throw std::runtime_error("readCSV: bad data header");

    const size_t totalCols = headerCols.size();
    const size_t dataCols = totalCols - 3;
    if (dataCols % (size_t)t.nDof != 0)
        throw std::runtime_error("readCSV: header columns not divisible by nDof");

    const size_t perModeCols = dataCols / (size_t)t.nDof;
    if (!(perModeCols == 2 || perModeCols == 4))
        throw std::runtime_error("readCSV: per-mode columns must be 2 or 4");

    // ---- 建 hash map 以 O(1) 定位索引 ----
    const double sFn = 1e9;
    const double sBe = 1e6;
    const double sOm = 1e9;

    std::unordered_map<long long, int> mapFn, mapBe, mapOm;
    mapFn.reserve(t.Fns.size() * 2);
    mapBe.reserve(t.dir.size() * 2);
    mapOm.reserve(t.w.size() * 2);

    for (int i = 0; i < (int)t.Fns.size(); ++i) mapFn[qkey(t.Fns[i], sFn)] = i;
    for (int i = 0; i < (int)t.dir.size(); ++i) mapBe[qkey(wrap360(t.dir[i]), sBe)] = i;
    for (int i = 0; i < (int)t.w.size(); ++i)   mapOm[qkey(t.w[i], sOm)] = i;

    // ---- data rows ----
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        const char* p = line.c_str();

        const double Fn = parseDoubleField(p);
        const double betaDeg = wrap360(parseDoubleField(p));
        const double om = parseDoubleField(p);

        auto itF = mapFn.find(qkey(Fn, sFn));
        auto itB = mapBe.find(qkey(betaDeg, sBe));
        auto itO = mapOm.find(qkey(om, sOm));

        if (itF == mapFn.end() || itB == mapBe.end() || itO == mapOm.end()) {
            // 如果你希望严格：直接 throw
            // throw std::runtime_error("readCSV: row axis not found on grid");
            continue;
        }

        const int iFn = itF->second;
        const int iBe = itB->second;
        const int iOm = itO->second;

        for (int k = 0; k < t.nDof; ++k) {
            const double re = parseDoubleField(p);
            const double im = parseDoubleField(p);

            if (!std::isnan(re) && !std::isnan(im))
                t.at(iFn, iBe, iOm, k) = { re, im };

            if (perModeCols == 4) {
                (void)parseDoubleField(p); // abs
                (void)parseDoubleField(p); // phase_deg
            }
        }
    }

    return t;
}
