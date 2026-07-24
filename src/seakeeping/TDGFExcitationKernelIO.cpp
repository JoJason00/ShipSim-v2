#include "TDGFExcitationKernelIO.h"
#include <filesystem>
#include <fstream>
#include <cstring>
#include <cstdint>
#include <stdexcept>

namespace
{
    // TDGFEXC1: original format, no region tag (loaded as region=Head).
    // TDGFEXC2: adds an int32 region tag at the end of the metadata block.
    constexpr const char* kMagicV1 = "TDGFEXC1";
    constexpr const char* kMagicV2 = "TDGFEXC2";
    constexpr std::size_t kMagicLen = 8;

    template <class T>
    void writePod(std::ofstream& os, const T& v)
    {
        os.write(reinterpret_cast<const char*>(&v), static_cast<std::streamsize>(sizeof(T)));
    }

    template <class T>
    void readPod(std::ifstream& is, T& v)
    {
        is.read(reinterpret_cast<char*>(&v), static_cast<std::streamsize>(sizeof(T)));
    }

    void writeIntVector(std::ofstream& os, const std::vector<int>& xs)
    {
        std::uint64_t n = static_cast<std::uint64_t>(xs.size());
        writePod(os, n);
        if (n > 0)
        {
            os.write(reinterpret_cast<const char*>(xs.data()),
                static_cast<std::streamsize>(n * sizeof(int)));
        }
    }

    void readIntVector(std::ifstream& is, std::vector<int>& xs)
    {
        std::uint64_t n = 0;
        readPod(is, n);
        xs.resize(static_cast<std::size_t>(n));
        if (n > 0)
        {
            is.read(reinterpret_cast<char*>(xs.data()),
                static_cast<std::streamsize>(n * sizeof(int)));
        }
    }

    void writeMatrix(std::ofstream& os, const Eigen::MatrixXd& M)
    {
        std::int32_t rows = static_cast<std::int32_t>(M.rows());
        std::int32_t cols = static_cast<std::int32_t>(M.cols());
        writePod(os, rows);
        writePod(os, cols);
        if (rows > 0 && cols > 0)
        {
            os.write(reinterpret_cast<const char*>(M.data()),
                static_cast<std::streamsize>(rows * cols * sizeof(double)));
        }
    }

    void readMatrix(std::ifstream& is, Eigen::MatrixXd& M)
    {
        std::int32_t rows = 0, cols = 0;
        readPod(is, rows);
        readPod(is, cols);
        if (rows < 0 || cols < 0)
            throw std::runtime_error("TDGFExcitationKernelIO: invalid matrix shape.");
        M.resize(rows, cols);
        if (rows > 0 && cols > 0)
        {
            is.read(reinterpret_cast<char*>(M.data()),
                static_cast<std::streamsize>(rows * cols * sizeof(double)));
        }
    }
}

bool TDGFExcitationKernelIO::save(const std::string& file, const TDGFExcitationKernelData& data)
{
    std::filesystem::path p(file);
    if (!p.parent_path().empty())
        std::filesystem::create_directories(p.parent_path());

    std::ofstream os(file, std::ios::binary);
    if (!os.is_open()) return false;

    char magic[16] = {};
    std::memcpy(magic, kMagicV2, kMagicLen);
    os.write(magic, sizeof(magic));

    writePod(os, data.Fn);
    writePod(os, data.U);
    writePod(os, data.betaRel);
    writePod(os, data.omegaIncident);
    writePod(os, data.omegaEncounter);
    writePod(os, data.dt);
    writePod(os, data.DOF);
    writeIntVector(os, data.modes);

    // V2-only: region tag (int32).
    const std::int32_t regionInt = static_cast<std::int32_t>(data.region);
    writePod(os, regionInt);

    writeMatrix(os, data.Qlag);

    return static_cast<bool>(os);
}

bool TDGFExcitationKernelIO::load(const std::string& file, TDGFExcitationKernelData& data)
{
    std::ifstream is(file, std::ios::binary);
    if (!is.is_open()) return false;

    char magic[16] = {};
    is.read(magic, sizeof(magic));
    if (!is) return false;

    const std::string magicStr(magic, magic + kMagicLen);
    const bool isV2 = (magicStr == std::string(kMagicV2));
    const bool isV1 = (magicStr == std::string(kMagicV1));
    if (!isV1 && !isV2) return false;

    readPod(is, data.Fn);
    readPod(is, data.U);
    readPod(is, data.betaRel);
    readPod(is, data.omegaIncident);
    readPod(is, data.omegaEncounter);
    readPod(is, data.dt);
    readPod(is, data.DOF);
    readIntVector(is, data.modes);

    if (isV2)
    {
        std::int32_t regionInt = 0;
        readPod(is, regionInt);
        data.region = static_cast<WaveForceRegion>(regionInt);
    }
    else
    {
        // Old files predate region tagging — default to Head.
        data.region = WaveForceRegion::Head;
    }

    readMatrix(is, data.Qlag);

    return static_cast<bool>(is);
}
