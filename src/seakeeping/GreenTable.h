#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

class Greenf;

class GreenTable
{
public:
    struct Params {
        // ---- main controls  ----
        double bmin = 1e-6, bmax = 300.0;
        double mmin = 0.0, mmax = 1.0;

        double mSplit1 = 0.005;
        double mSplit2 = 0.01;
        double mSplit3 = 0.1;

        bool throwOnOutOfRange = true;

        // ---- derived  ----
        int Nm1 = 0, Nm2 = 0, Nm3 = 0, Nm4 = 0;
        int Nb1 = 0, Nb2 = 0, Nb3 = 0, Nb4 = 0;

        // ---- b mapping parameters ----
        double pEndMax = 2.2;   // 段1：endDense（b大更密）
        double pEndMin = 1.0;   // 段2：均匀（务必 1.0）
        double pStartMin = 1.65;  // 段3：startDense（前密）
        double pStartMax = 2.2;   // 段4：startDense（前密更强）

        // ---- m segment shaping (buildMnodes_ 用) ----
        double mPow1 = 1.6;
        double mPow2 = 1.3;
        double mPow3 = 1.1;
        double mPow4 = 1.0;

        // ---- misc (cpp里用到的) ----
        int    mLUT = 512;     // buildMLUT_ 的 LUT 长度
        double m_eps = 1e-6;   // m<=m_eps 时 clamp
        bool   storeFloat = true;

        // ========= autoTune knobs  =========
        // Nm 分配
        int NmTotal = 660;                 // 默认 220+160+160+120
        double nmFrac1 = 220.0 / 660.0;
        double nmFrac2 = 160.0 / 660.0;
        double nmFrac3 = 160.0 / 660.0;      // Nm4 = remainder

        // Nb1 随 bmax 缩放：Nb1 = Nb1Ref*(bmax/bRef)^nbPow
        double bRef = 300.0;
        int    Nb1Ref = 30000;
        double nbPow = 1.0;                // 1=线性；更激进可 2

        // Nb2..4 相对 Nb1 比例（默认保持 30000→6000→3000→2000）
        double nbFrac2 = 0.20;
        double nbFrac3 = 0.10;
        double nbFrac4 = 1.0 / 15.0;

        // 兜底下限（避免中后段太稀）
        int Nb2Min = 1500;
        int Nb3Min = 1200;
        int Nb4Min = 1000;

        // ---- API ----
        static Params AdaptiveDefault(double bmax_ = 300.0) {
            Params p;
            p.bmax = bmax_;
            p.autoTune();
            return p;
        }

        void autoTune() {
            auto clampInt = [](int x, int lo) { return x < lo ? lo : x; };

            // ---- Nm1..4 ----
            int n1 = (int)std::llround(NmTotal * nmFrac1);
            int n2 = (int)std::llround(NmTotal * nmFrac2);
            int n3 = (int)std::llround(NmTotal * nmFrac3);
            int n4 = NmTotal - n1 - n2 - n3;

            n1 = clampInt(n1, 2);
            n2 = clampInt(n2, 2);
            n3 = clampInt(n3, 2);
            n4 = clampInt(n4, 2);

            int sum = n1 + n2 + n3 + n4;
            if (sum != NmTotal) n4 += (NmTotal - sum);

            Nm1 = n1; Nm2 = n2; Nm3 = n3; Nm4 = n4;

            // ---- Nb1..4 ----
            double scale = (bRef > 0.0) ? (bmax / bRef) : 1.0;
            double nb1 = double(Nb1Ref) * std::pow(scale, nbPow);
            Nb1 = (int)std::llround(nb1);
            if (Nb1 < 2) Nb1 = 2;

            Nb2 = (int)std::llround(Nb1 * nbFrac2);
            Nb3 = (int)std::llround(Nb1 * nbFrac3);
            Nb4 = (int)std::llround(Nb1 * nbFrac4);

            if (Nb2 < Nb2Min) Nb2 = Nb2Min;
            if (Nb3 < Nb3Min) Nb3 = Nb3Min;
            if (Nb4 < Nb4Min) Nb4 = Nb4Min;

            // 单调
            if (Nb2 > Nb1) Nb2 = Nb1;
            if (Nb3 > Nb2) Nb3 = Nb2;
            if (Nb4 > Nb3) Nb4 = Nb3;

            // p 值合法
            pEndMin = std::max(1.0, pEndMin);
            pEndMax = std::max(1.0, pEndMax);
            pStartMin = std::max(1.0, pStartMin);
            pStartMax = std::max(1.0, pStartMax);
        }
    };

public:
    GreenTable() = default;
    explicit GreenTable(Params p) : p_(std::move(p)) { p_.autoTune(); }

    void init(const std::string& binPath, double gConst, bool strict);
    void eval(double b, double m, double& outGf, double& outGbd, double& outGmd) const;

    bool ready() const { return ready_; }
    const Params& params() const { return p_; }

private:
#pragma pack(push, 1)
    struct Header
    {
        char     magic[8];     // "GTBLAD2\0"
        uint32_t version;      // 2
        uint32_t flags;        // bit0: storeFloat
        uint32_t Nm;
        uint32_t mLUT;

        double bmin, bmax, mmin, mmax;
        double mSplit1, mSplit2;  // split3 stored in reserved
        uint32_t Nm1, Nm2, Nm3;   // Nm4 stored in reserved
        double mPow1, mPow2, mPow3; // mPow4 stored in reserved
        uint32_t Nb1, Nb2, Nb3;   // Nb4 stored in reserved

        double pEndMin, pEndMax;
        double pStartMin, pStartMax;

        double m_eps;
        double gConst;

        uint64_t totalN;

        uint64_t reserved[8]{};   // store split3/Nm4/mPow4/Nb4 here
    };
#pragma pack(pop)

private:
    Params p_{};
    bool ready_ = false;

    double smin_ = 0.0, smax_ = 0.0;
    double inv_sRange_ = 0.0;

    // m grid
    std::vector<double>   mnode_;
    std::vector<double>   inv_dm_;
    std::vector<uint32_t> mLUT_;

    // per-row
    std::vector<uint32_t> rowNb_;
    std::vector<uint64_t> rowOff_;
    std::vector<uint8_t>  rowBType_; // 0=startDense, 1=endDense
    std::vector<float>    rowBPow_;  // exponent p for that row

    // data arrays
    std::vector<float> Gf_;
    std::vector<float> Gbd_;
    std::vector<float> Gmd_;

private:
    static bool fileExists_(const std::string& path);

    template<class T>
    static void writePod_(std::ofstream& os, const T& v);

    template<class T>
    static void readPod_(std::ifstream& is, T& v);

    void checkBounds_(double b, double m) const;

    // mapping
    static double mapU_(double t, uint8_t type, double p);
    static double invMapU_(double u, uint8_t type, double p);
    double nodeS_(uint32_t k, uint32_t NbRow, uint8_t type, double p) const;

    // build
    void buildGrid_();
    void buildMnodes_();
    void buildMLUT_();

    uint32_t NbForM_(double m) const;
    void bMapForM_(double m, uint8_t& type, double& p) const;

    int locateMfast_(double m) const;
    void interpRow_(int im, double s, double& gf, double& gbd, double& gmd) const;

    void build_(Greenf& exact);
    void save_(const std::string& path, double gConst) const;
    void load_(const std::string& path, double gConst, bool strict);
};
