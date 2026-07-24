#include "GreenTable.h"
#include "greenf.h"

#include <fstream>
#include <stdexcept>
#include <cmath>
#include <cstring>
#include <algorithm>

static inline double clamp01(double x)
{
    return std::min(1.0, std::max(0.0, x));
}

static inline uint64_t packDouble64(double x)
{
    uint64_t u;
    std::memcpy(&u, &x, sizeof(double));
    return u;
}

static inline double unpackDouble64(uint64_t u)
{
    double x;
    std::memcpy(&x, &u, sizeof(double));
    return x;
}

bool GreenTable::fileExists_(const std::string& path)
{
    std::ifstream is(path, std::ios::binary);
    return (bool)is;
}

template<class T>
void GreenTable::writePod_(std::ofstream& os, const T& v)
{
    os.write(reinterpret_cast<const char*>(&v), sizeof(T));
    if (!os) throw std::runtime_error("GreenTable: write failed");
}

template<class T>
void GreenTable::readPod_(std::ifstream& is, T& v)
{
    is.read(reinterpret_cast<char*>(&v), sizeof(T));
    if (!is) throw std::runtime_error("GreenTable: read failed");
}

void GreenTable::checkBounds_(double b, double m) const
{
    if (!p_.throwOnOutOfRange) return;
    if (b < p_.bmin) throw std::out_of_range("GreenTable: b < bmin");
    if (b > p_.bmax) throw std::out_of_range("GreenTable: b > bmax");
    if (m < p_.mmin) throw std::out_of_range("GreenTable: m < mmin");
    if (m > p_.mmax) throw std::out_of_range("GreenTable: m > mmax");
}

// type=0 startDense: u=t^p
// type=1 endDense  : u=1-(1-t)^p
double GreenTable::mapU_(double t, uint8_t type, double p)
{
    t = clamp01(t);
    if (type == 0) return std::pow(t, p);
    return 1.0 - std::pow(1.0 - t, p);
}

// inverse mapping: given u -> t
double GreenTable::invMapU_(double u, uint8_t type, double p)
{
    u = clamp01(u);
    const double invp = 1.0 / p;
    if (type == 0) return std::pow(u, invp);
    return 1.0 - std::pow(1.0 - u, invp);
}

double GreenTable::nodeS_(uint32_t k, uint32_t NbRow, uint8_t type, double p) const
{
    const double t = double(k) / double(NbRow - 1);
    const double u = mapU_(t, type, p);
    return smin_ + (smax_ - smin_) * u;
}

uint32_t GreenTable::NbForM_(double m) const
{
    if (m <= p_.mSplit1) return (uint32_t)std::max(2, p_.Nb1);
    if (m <= p_.mSplit2) return (uint32_t)std::max(2, p_.Nb2);
    if (m <= p_.mSplit3) return (uint32_t)std::max(2, p_.Nb3);
    return (uint32_t)std::max(2, p_.Nb4);
}

void GreenTable::bMapForM_(double m, uint8_t& type, double& p) const
{
    if (m <= p_.mSplit1) { type = 1; p = std::max(1.0, p_.pEndMax); return; }
    if (m <= p_.mSplit2) { type = 0; p = std::max(1.0, p_.pEndMin); return; }
    if (m <= p_.mSplit3) { type = 0; p = std::max(1.0, p_.pStartMin); return; }
    type = 0; p = std::max(1.0, p_.pStartMax);
}

void GreenTable::buildMnodes_()
{
    if (!(p_.bmin > 0.0) || !(p_.bmax > p_.bmin))
        throw std::runtime_error("GreenTable: invalid b range");

    if (!(p_.mmin >= 0.0) || !(p_.mmax > p_.mmin) || !(p_.mmax <= 1.0))
        throw std::runtime_error("GreenTable: invalid m range");

    if (!(p_.mSplit1 > p_.mmin &&
        p_.mSplit2 > p_.mSplit1 &&
        p_.mSplit3 > p_.mSplit2 &&
        p_.mSplit3 < p_.mmax))
        throw std::runtime_error("GreenTable: invalid m splits (need mSplit1<mSplit2<mSplit3)");

    if (p_.Nm1 < 2 || p_.Nm2 < 2 || p_.Nm3 < 2 || p_.Nm4 < 2)
        throw std::runtime_error("GreenTable: Nm per segment too small");

    auto segPow = [](double a, double b, int N, double powp,
        bool includeLeft, bool includeRight,
        std::vector<double>& out)
    {
        const int i0 = includeLeft ? 0 : 1;
        const int i1 = includeRight ? (N - 1) : (N - 2);
        for (int i = i0; i <= i1; ++i) {
            double t = double(i) / double(N - 1);
            double u = std::pow(t, powp);
            out.push_back(a + (b - a) * u);
        }
    };

    mnode_.clear();
    mnode_.reserve((size_t)p_.Nm1 + (size_t)p_.Nm2 + (size_t)p_.Nm3 + (size_t)p_.Nm4 - 3);

    segPow(p_.mmin, p_.mSplit1, p_.Nm1, p_.mPow1, true, true, mnode_);
    segPow(p_.mSplit1, p_.mSplit2, p_.Nm2, p_.mPow2, false, true, mnode_);
    segPow(p_.mSplit2, p_.mSplit3, p_.Nm3, p_.mPow3, false, true, mnode_);
    segPow(p_.mSplit3, p_.mmax, p_.Nm4, p_.mPow4, false, true, mnode_);

    if ((int)mnode_.size() < 2)
        throw std::runtime_error("GreenTable: Nm too small after build");
}

void GreenTable::buildMLUT_()
{
    int L = std::max(64, p_.mLUT);
    mLUT_.assign((size_t)L, 0);

    int im = 0;
    const int Nm = (int)mnode_.size();
    for (int j = 0; j < L; ++j) {
        double u = (L == 1) ? 0.0 : double(j) / double(L - 1);
        double m = p_.mmin + (p_.mmax - p_.mmin) * u;

        while (im < Nm - 2 && m > mnode_[im + 1]) ++im;
        while (im > 0 && m < mnode_[im])     --im;

        mLUT_[j] = (uint32_t)std::clamp(im, 0, Nm - 2);
    }
}

int GreenTable::locateMfast_(double m) const
{
    const int Nm = (int)mnode_.size();
    if (Nm < 2) return 0;

    double u = (m - p_.mmin) / (p_.mmax - p_.mmin);
    u = clamp01(u);

    const int L = (int)mLUT_.size();
    int j = (int)std::floor(u * double(L - 1));
    j = std::clamp(j, 0, L - 1);

    int im = (int)mLUT_[j];
    im = std::clamp(im, 0, Nm - 2);

    while (im > 0 && m < mnode_[im])     --im;
    while (im < Nm - 2 && m > mnode_[im + 1]) ++im;
    return im;
}

void GreenTable::buildGrid_()
{
    buildMnodes_();

    smin_ = p_.bmin * p_.bmin;
    smax_ = p_.bmax * p_.bmax;
    inv_sRange_ = 1.0 / (smax_ - smin_);

    const int Nm = (int)mnode_.size();

    inv_dm_.resize((size_t)Nm - 1);
    for (int im = 0; im < Nm - 1; ++im) {
        double dm = mnode_[im + 1] - mnode_[im];
        inv_dm_[im] = (dm > 0.0) ? (1.0 / dm) : 0.0;
    }

    buildMLUT_();

    rowNb_.resize((size_t)Nm);
    rowBType_.resize((size_t)Nm);
    rowBPow_.resize((size_t)Nm);
    rowOff_.resize((size_t)Nm + 1);
    rowOff_[0] = 0;

    for (int im = 0; im < Nm; ++im) {
        const double m = mnode_[im];

        uint32_t NbRow = NbForM_(m);
        rowNb_[im] = NbRow;

        uint8_t type = 0;
        double  bp = 1.0;
        bMapForM_(m, type, bp);
        rowBType_[im] = type;
        rowBPow_[im] = (float)bp;

        rowOff_[im + 1] = rowOff_[im] + (uint64_t)NbRow;
    }

    const uint64_t totalN = rowOff_.back();
    if (totalN == 0) throw std::runtime_error("GreenTable: totalN == 0");

    Gf_.assign((size_t)totalN, 0.0f);
    Gbd_.assign((size_t)totalN, 0.0f);
    Gmd_.assign((size_t)totalN, 0.0f);
}

void GreenTable::interpRow_(int im, double s, double& gf, double& gbd, double& gmd) const
{
    const uint32_t NbRow = rowNb_[im];
    const uint64_t off = rowOff_[im];

    const uint8_t type = rowBType_[im];
    const double  p = (double)rowBPow_[im];

    double u = (s - smin_) * inv_sRange_;
    u = clamp01(u);

    double t = invMapU_(u, type, p);
    double x = t * double(NbRow - 1);
    int k = (int)std::floor(x);
    if (k < 0) k = 0;
    if (k > (int)NbRow - 2) k = (int)NbRow - 2;

    double s0 = nodeS_((uint32_t)k, NbRow, type, p);
    double s1 = nodeS_((uint32_t)k + 1, NbRow, type, p);
    double tb = (s1 > s0) ? ((s - s0) / (s1 - s0)) : 0.0;
    tb = clamp01(tb);

    const size_t i0 = (size_t)off + (size_t)k;
    const size_t i1 = i0 + 1;

    auto lerp = [](double a, double b, double t) { return a + (b - a) * t; };
    gf = lerp(Gf_[i0], Gf_[i1], tb);
    gbd = lerp(Gbd_[i0], Gbd_[i1], tb);
    gmd = lerp(Gmd_[i0], Gmd_[i1], tb);
}

void GreenTable::eval(double b, double m, double& outGf, double& outGbd, double& outGmd) const
{
    if (!ready_) throw std::runtime_error("GreenTable: not ready");
    checkBounds_(b, m);

    if (!p_.throwOnOutOfRange) {
        b = std::clamp(b, p_.bmin, p_.bmax);
        m = std::clamp(m, p_.mmin, p_.mmax);
    }

    if (m <= p_.m_eps) m = p_.m_eps;

    const double s = b * b;

    int im0 = locateMfast_(m);
    int im1 = im0 + 1;

    double g0, gb0, gm0;
    double g1, gb1, gm1;
    interpRow_(im0, s, g0, gb0, gm0);
    interpRow_(im1, s, g1, gb1, gm1);

    double tm = (m - mnode_[im0]) * inv_dm_[im0];
    tm = clamp01(tm);

    auto lerp = [](double a, double b, double t) { return a + (b - a) * t; };
    outGf = lerp(g0, g1, tm);
    outGbd = lerp(gb0, gb1, tm);
    outGmd = lerp(gm0, gm1, tm);
}

void GreenTable::build_(Greenf& exact)
{
    buildGrid_();

    const int Nm = (int)mnode_.size();

    for (int im = 0; im < Nm; ++im) {
        double m = mnode_[im];
        if (m <= p_.m_eps) m = p_.m_eps;

        const uint32_t NbRow = rowNb_[im];
        const uint64_t off = rowOff_[im];

        const uint8_t type = rowBType_[im];
        const double  p = (double)rowBPow_[im];

        for (uint32_t k = 0; k < NbRow; ++k) {
            double t = double(k) / double(NbRow - 1);
            double u = mapU_(t, type, p);
            double s = smin_ + (smax_ - smin_) * u;
            double b = std::sqrt(s);

            exact.GreenFunctionCal(b, m);

            size_t idx = (size_t)off + (size_t)k;
            Gf_[idx] = (float)exact.Gf;
            Gbd_[idx] = (float)exact.Gbd;
            Gmd_[idx] = (float)exact.Gmd;
        }
    }

    ready_ = true;
}

static inline bool nearlyEq(double a, double b, double tol = 1e-14)
{
    return std::abs(a - b) <= tol;
}

void GreenTable::save_(const std::string& path, double gConst) const
{
    if (!ready_) throw std::runtime_error("GreenTable: not ready");

    Header h{};
    std::memcpy(h.magic, "GTBLAD2\0", 8);
    h.version = 2;
    h.flags = (p_.storeFloat ? 1u : 0u);

    h.Nm = (uint32_t)mnode_.size();
    h.mLUT = (uint32_t)mLUT_.size();

    h.bmin = p_.bmin; h.bmax = p_.bmax;
    h.mmin = p_.mmin; h.mmax = p_.mmax;

    h.mSplit1 = p_.mSplit1;
    h.mSplit2 = p_.mSplit2;

    h.Nm1 = (uint32_t)p_.Nm1;
    h.Nm2 = (uint32_t)p_.Nm2;
    h.Nm3 = (uint32_t)p_.Nm3;

    h.mPow1 = p_.mPow1;
    h.mPow2 = p_.mPow2;
    h.mPow3 = p_.mPow3;

    h.Nb1 = (uint32_t)p_.Nb1;
    h.Nb2 = (uint32_t)p_.Nb2;
    h.Nb3 = (uint32_t)p_.Nb3;

    h.pEndMin = p_.pEndMin;
    h.pEndMax = p_.pEndMax;
    h.pStartMin = p_.pStartMin;
    h.pStartMax = p_.pStartMax;

    h.m_eps = p_.m_eps;
    h.gConst = gConst;
    h.totalN = rowOff_.back();

    // ---- store extra 4-seg params in reserved ----
    h.reserved[0] = packDouble64(p_.mSplit3);
    h.reserved[1] = (uint64_t)(uint32_t)p_.Nm4;
    h.reserved[2] = packDouble64(p_.mPow4);
    h.reserved[3] = (uint64_t)(uint32_t)p_.Nb4;

    std::ofstream os(path, std::ios::binary);
    if (!os) throw std::runtime_error("GreenTable: cannot open for write: " + path);

    writePod_(os, h);

    os.write((const char*)mnode_.data(), (std::streamsize)(mnode_.size() * sizeof(double)));
    os.write((const char*)inv_dm_.data(), (std::streamsize)(inv_dm_.size() * sizeof(double)));
    os.write((const char*)mLUT_.data(), (std::streamsize)(mLUT_.size() * sizeof(uint32_t)));

    os.write((const char*)rowNb_.data(), (std::streamsize)(rowNb_.size() * sizeof(uint32_t)));
    os.write((const char*)rowOff_.data(), (std::streamsize)(rowOff_.size() * sizeof(uint64_t)));
    os.write((const char*)rowBType_.data(), (std::streamsize)(rowBType_.size() * sizeof(uint8_t)));
    os.write((const char*)rowBPow_.data(), (std::streamsize)(rowBPow_.size() * sizeof(float)));

    os.write((const char*)Gf_.data(), (std::streamsize)(Gf_.size() * sizeof(float)));
    os.write((const char*)Gbd_.data(), (std::streamsize)(Gbd_.size() * sizeof(float)));
    os.write((const char*)Gmd_.data(), (std::streamsize)(Gmd_.size() * sizeof(float)));

    if (!os) throw std::runtime_error("GreenTable: write failed: " + path);
}

void GreenTable::load_(const std::string& path, double gConst, bool strict)
{
    std::ifstream is(path, std::ios::binary);
    if (!is) throw std::runtime_error("GreenTable: cannot open for read: " + path);

    Header h{};
    readPod_(is, h);

    if (std::memcmp(h.magic, "GTBLAD2\0", 8) != 0 || h.version != 2)
        throw std::runtime_error("GreenTable: old/unknown table format -> rebuild");

    if (strict && !nearlyEq(h.gConst, gConst, 1e-12))
        throw std::runtime_error("GreenTable: gConst mismatch -> rebuild");

    const double fileSplit3 = unpackDouble64(h.reserved[0]);
    const int    fileNm4 = (int)(uint32_t)h.reserved[1];
    const double filePow4 = unpackDouble64(h.reserved[2]);
    const int    fileNb4 = (int)(uint32_t)h.reserved[3];

    if (!(fileSplit3 > h.mSplit2 && fileSplit3 < h.mmax))
        throw std::runtime_error("GreenTable: file has no valid mSplit3 -> rebuild");

    if (strict) {
        if (!nearlyEq(h.bmin, p_.bmin) || !nearlyEq(h.bmax, p_.bmax) ||
            !nearlyEq(h.mmin, p_.mmin) || !nearlyEq(h.mmax, p_.mmax) ||
            !nearlyEq(h.mSplit1, p_.mSplit1) || !nearlyEq(h.mSplit2, p_.mSplit2) ||
            !nearlyEq(fileSplit3, p_.mSplit3) ||
            h.Nm1 != (uint32_t)p_.Nm1 || h.Nm2 != (uint32_t)p_.Nm2 || h.Nm3 != (uint32_t)p_.Nm3 || fileNm4 != p_.Nm4 ||
            !nearlyEq(h.mPow1, p_.mPow1) || !nearlyEq(h.mPow2, p_.mPow2) || !nearlyEq(h.mPow3, p_.mPow3) || !nearlyEq(filePow4, p_.mPow4) ||
            h.Nb1 != (uint32_t)p_.Nb1 || h.Nb2 != (uint32_t)p_.Nb2 || h.Nb3 != (uint32_t)p_.Nb3 || fileNb4 != p_.Nb4 ||
            !nearlyEq(h.pEndMin, p_.pEndMin) || !nearlyEq(h.pEndMax, p_.pEndMax) ||
            !nearlyEq(h.pStartMin, p_.pStartMin) || !nearlyEq(h.pStartMax, p_.pStartMax) ||
            !nearlyEq(h.m_eps, p_.m_eps) ||
            h.mLUT != (uint32_t)p_.mLUT)
        {
            throw std::runtime_error("GreenTable: params mismatch -> rebuild");
        }
    }

    // adopt from file
    p_.bmin = h.bmin; p_.bmax = h.bmax;
    p_.mmin = h.mmin; p_.mmax = h.mmax;

    p_.mSplit1 = h.mSplit1;
    p_.mSplit2 = h.mSplit2;
    p_.mSplit3 = fileSplit3;

    p_.Nm1 = (int)h.Nm1; p_.Nm2 = (int)h.Nm2; p_.Nm3 = (int)h.Nm3; p_.Nm4 = fileNm4;
    p_.mPow1 = h.mPow1; p_.mPow2 = h.mPow2; p_.mPow3 = h.mPow3; p_.mPow4 = filePow4;

    p_.Nb1 = (int)h.Nb1; p_.Nb2 = (int)h.Nb2; p_.Nb3 = (int)h.Nb3; p_.Nb4 = fileNb4;

    p_.pEndMin = h.pEndMin;
    p_.pEndMax = h.pEndMax;
    p_.pStartMin = h.pStartMin;
    p_.pStartMax = h.pStartMax;

    p_.m_eps = h.m_eps;
    p_.mLUT = (int)h.mLUT;
    p_.storeFloat = (h.flags & 1u) != 0u;

    smin_ = p_.bmin * p_.bmin;
    smax_ = p_.bmax * p_.bmax;
    inv_sRange_ = 1.0 / (smax_ - smin_);

    const uint32_t Nm = h.Nm;

    mnode_.resize(Nm);
    inv_dm_.resize(Nm > 0 ? (Nm - 1) : 0);
    mLUT_.resize(h.mLUT);

    rowNb_.resize(Nm);
    rowOff_.resize((size_t)Nm + 1);
    rowBType_.resize(Nm);
    rowBPow_.resize(Nm);

    is.read((char*)mnode_.data(), (std::streamsize)(mnode_.size() * sizeof(double)));
    is.read((char*)inv_dm_.data(), (std::streamsize)(inv_dm_.size() * sizeof(double)));
    is.read((char*)mLUT_.data(), (std::streamsize)(mLUT_.size() * sizeof(uint32_t)));

    is.read((char*)rowNb_.data(), (std::streamsize)(rowNb_.size() * sizeof(uint32_t)));
    is.read((char*)rowOff_.data(), (std::streamsize)(rowOff_.size() * sizeof(uint64_t)));
    is.read((char*)rowBType_.data(), (std::streamsize)(rowBType_.size() * sizeof(uint8_t)));
    is.read((char*)rowBPow_.data(), (std::streamsize)(rowBPow_.size() * sizeof(float)));

    const uint64_t totalN = h.totalN;
    if (rowOff_.back() != totalN) throw std::runtime_error("GreenTable: totalN mismatch/corrupt file");

    Gf_.resize((size_t)totalN);
    Gbd_.resize((size_t)totalN);
    Gmd_.resize((size_t)totalN);

    is.read((char*)Gf_.data(), (std::streamsize)(Gf_.size() * sizeof(float)));
    is.read((char*)Gbd_.data(), (std::streamsize)(Gbd_.size() * sizeof(float)));
    is.read((char*)Gmd_.data(), (std::streamsize)(Gmd_.size() * sizeof(float)));

    if (!is) throw std::runtime_error("GreenTable: read data failed");

    ready_ = true;
}

void GreenTable::init(const std::string& binPath, double gConst, bool strict)
{
    if (fileExists_(binPath)) {
        try {
            load_(binPath, gConst, strict);
            return;
        }
        catch (...) {
            // rebuild
        }
    }

    Greenf exact;
    build_(exact);
    save_(binPath, gConst);
}
