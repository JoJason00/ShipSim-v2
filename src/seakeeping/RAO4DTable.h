#pragma once
#include <vector>
#include <complex>
#include <string>
#include <cstdint>

struct RAO4DMeta
{
    enum class DirUnit : uint32_t { Deg = 0, Rad = 1 };
    enum class FreqUnit : uint32_t { RadPerSec = 0, Hz = 1 };
    enum class FreqType : uint32_t { IncidentOmega = 0, EncounterOmega = 1 };
    enum class NormType : uint32_t { Raw = 0, TransA_RotkA = 1 };

    DirUnit  dirUnit = DirUnit::Deg;          // 表内 dir 统一用 Deg
    FreqUnit freqUnit = FreqUnit::RadPerSec;   // 表内 w 统一用 rad/s
    FreqType freqType = FreqType::IncidentOmega;
    NormType normType = NormType::TransA_RotkA;

    double L = 0.0;
};

struct RAO4DTable
{
    RAO4DMeta meta;

    std::vector<int> modeIds;       // e.g. [2,3,4]
    std::vector<double> Fns;        // Fn axis
    std::vector<double> dir;        // beta axis in DEG [0,360)
    std::vector<double> w;          // omega axis rad/s

    // H[Fn][beta][omega][dof]
    std::vector<std::complex<double>> H;

    int nFn = 0, nBe = 0, nOm = 0, nDof = 0;
    size_t strideOm = 0, strideBe = 0, strideFn = 0;

    void finalizeAndAllocate(bool fillNaN = true);

    inline std::complex<double>& at(int iFn_, int iBe_, int iOm_, int k) {
        return H[(size_t)iFn_ * strideFn
            + (size_t)iBe_ * strideBe
            + (size_t)iOm_ * strideOm
            + (size_t)k];
    }
    inline const std::complex<double>& at(int iFn_, int iBe_, int iOm_, int k) const {
        return H[(size_t)iFn_ * strideFn
            + (size_t)iBe_ * strideBe
            + (size_t)iOm_ * strideOm
            + (size_t)k];
    }

    void setAt(int iFn_, int iBe_, int iOm_, const std::vector<std::complex<double>>& Hk);

    // betaDeg 输入是 DEG
    std::vector<std::complex<double>> interp(double Fn, double betaDeg, double omega) const;

    // ---- CSV IO (single file: meta + axes + data) ----
    void writeCSV(const std::string& path,  std::vector<double>& non_we, int precision = 12, bool includeAbsPhase = true) const;
    static RAO4DTable readCSV(const std::string& path);

    // ---- diagnostics ----
    size_t countNaN() const;
    bool isComplete() const { return countNaN() == 0; }

    // ---- util ----
    static bool fileExists(const std::string& path);
};
