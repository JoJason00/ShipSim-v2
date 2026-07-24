#include "RadiationKernelCache.h"

#include <fstream>
#include <filesystem>
#include <cstdint>
#include <stdexcept>
#include <cstring>
#include <iostream>

namespace
{
    // 新格式 v8：自描述（含船名/Lpp/Klag_times），支持非均匀 K
    constexpr const char* kMagicV8 = "RKCACHE8_NUNI ";   // 14 chars
    // 旧格式 v7：等距 K，无 ship/Lpp/times（仍可读，加载时自动补 times = i*dt）
    constexpr const char* kMagicV7 = "RKCACHE7_VEL34";
    constexpr std::size_t kMagicLen = 14;

    // 防御性上限：任何字段超过这些值就视为文件损坏，直接拒绝
    constexpr std::int32_t kMaxMatrixDim    = 1024;        // DOF 矩阵通常 ≤ 32
    constexpr std::uint64_t kMaxVectorLen   = 10'000'000;  // Klag 长度上限
    constexpr std::uint64_t kMaxStringLen   = 4096;        // ShipName 长度上限

    template<typename T>
    void writePod(std::ofstream& os, const T& v)
    {
        os.write(reinterpret_cast<const char*>(&v), sizeof(T));
    }

    template<typename T>
    void readPod(std::ifstream& is, T& v)
    {
        is.read(reinterpret_cast<char*>(&v), sizeof(T));
    }

    void writeIntVector(std::ofstream& os, const std::vector<int>& v)
    {
        std::uint64_t n = static_cast<std::uint64_t>(v.size());
        writePod(os, n);
        if (n > 0)
        {
            os.write(reinterpret_cast<const char*>(v.data()),
                static_cast<std::streamsize>(n * sizeof(int)));
        }
    }

    void readIntVector(std::ifstream& is, std::vector<int>& v)
    {
        std::uint64_t n = 0;
        readPod(is, n);
        if (!is.good() || n > kMaxVectorLen)
            throw std::runtime_error(
                "RadiationKernelCache: invalid int-vector length (n=" +
                std::to_string(n) + ")");
        v.resize(static_cast<std::size_t>(n));
        if (n > 0)
        {
            is.read(reinterpret_cast<char*>(v.data()),
                static_cast<std::streamsize>(n * sizeof(int)));
        }
    }

    void writeString(std::ofstream& os, const std::string& s)
    {
        std::uint64_t n = static_cast<std::uint64_t>(s.size());
        writePod(os, n);
        if (n > 0)
            os.write(s.data(), static_cast<std::streamsize>(n));
    }

    void readString(std::ifstream& is, std::string& s)
    {
        std::uint64_t n = 0;
        readPod(is, n);
        if (!is.good() || n > kMaxStringLen)
            throw std::runtime_error(
                "RadiationKernelCache: invalid string length (n=" +
                std::to_string(n) + ")");
        s.resize(static_cast<std::size_t>(n));
        if (n > 0)
            is.read(s.data(), static_cast<std::streamsize>(n));
    }

    void writeDoubleVector(std::ofstream& os, const std::vector<double>& v)
    {
        std::uint64_t n = static_cast<std::uint64_t>(v.size());
        writePod(os, n);
        if (n > 0)
        {
            os.write(reinterpret_cast<const char*>(v.data()),
                static_cast<std::streamsize>(n * sizeof(double)));
        }
    }

    void readDoubleVector(std::ifstream& is, std::vector<double>& v)
    {
        std::uint64_t n = 0;
        readPod(is, n);
        if (!is.good() || n > kMaxVectorLen)
            throw std::runtime_error(
                "RadiationKernelCache: invalid double-vector length (n=" +
                std::to_string(n) + ")");
        v.resize(static_cast<std::size_t>(n));
        if (n > 0)
        {
            is.read(reinterpret_cast<char*>(v.data()),
                static_cast<std::streamsize>(n * sizeof(double)));
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

        if (!is.good() || rows < 0 || cols < 0 ||
            rows > kMaxMatrixDim || cols > kMaxMatrixDim)
        {
            throw std::runtime_error(
                "RadiationKernelCache: invalid matrix shape (rows=" +
                std::to_string(rows) + ", cols=" + std::to_string(cols) + ")");
        }

        M.resize(rows, cols);
        if (rows > 0 && cols > 0)
        {
            is.read(reinterpret_cast<char*>(M.data()),
                static_cast<std::streamsize>(rows * cols * sizeof(double)));
        }
    }
}

bool RadiationKernelCache::save(const std::string& file, const RadiationKernelData& data)
{
    std::filesystem::path p(file);
    if (!p.parent_path().empty())
        std::filesystem::create_directories(p.parent_path());

    std::ofstream os(file, std::ios::binary);
    if (!os.is_open()) return false;

    // ---- v8 自描述格式 ----
    char magic[16] = {};
    std::memcpy(magic, kMagicV8, kMagicLen);
    os.write(magic, sizeof(magic));

    // 案例指纹
    writeString(os, data.ShipName);
    writePod(os, data.Fn);
    writePod(os, data.U);
    writePod(os, data.Lpp);
    writePod(os, data.dt);

    writePod(os, data.TG);
    writePod(os, data.NE);
    writePod(os, data.DOF);

    writeIntVector(os, data.modes);

    writeMatrix(os, data.A_inf);
    writeMatrix(os, data.B);
    writeMatrix(os, data.C_prime);
    writeMatrix(os, data.K0);

    // Klag_times 单独写出（与 Klag 等长）
    writeDoubleVector(os, data.Klag_times);

    std::uint64_t nk = static_cast<std::uint64_t>(data.Klag.size());
    writePod(os, nk);
    for (const auto& K : data.Klag)
        writeMatrix(os, K);

    return static_cast<bool>(os);
}

bool RadiationKernelCache::load(const std::string& file, RadiationKernelData& data)
{
    //namespace fs = std::filesystem;

    //fs::path p(file);

    //std::cout << "[RadiationKernelCache] load file = [" << file << "]\n";
    //std::cout << "[RadiationKernelCache] current_path = ["
    //    << fs::current_path().string() << "]\n";
    //std::cout << "[RadiationKernelCache] absolute path = ["
    //    << fs::absolute(p).string() << "]\n";
    //std::cout << "[RadiationKernelCache] exists = "
    //    << fs::exists(p) << "\n";
    //std::cout << "[RadiationKernelCache] is_regular_file = "
    //    << fs::is_regular_file(p) << "\n";

    std::ifstream is(file, std::ios::binary);
    if (!is.is_open()) return false;

    char magic[16] = {};
    is.read(magic, sizeof(magic));
    if (!is) return false;

    const std::string magicStr(magic, magic + kMagicLen);
    const bool isV8 = (magicStr == std::string(kMagicV8));
    const bool isV7 = (magicStr == std::string(kMagicV7));
    if (!isV8 && !isV7)
        return false;

    try
    {
        if (isV8)
        {
            // ---- 自描述格式 ----
            readString(is, data.ShipName);
            readPod(is, data.Fn);
            readPod(is, data.U);
            readPod(is, data.Lpp);
            readPod(is, data.dt);

            readPod(is, data.TG);
            readPod(is, data.NE);
            readPod(is, data.DOF);

            readIntVector(is, data.modes);

            readMatrix(is, data.A_inf);
            readMatrix(is, data.B);
            readMatrix(is, data.C_prime);
            readMatrix(is, data.K0);

            readDoubleVector(is, data.Klag_times);

            std::uint64_t nk = 0;
            readPod(is, nk);
            if (!is.good() || nk > kMaxVectorLen)
                throw std::runtime_error(
                    "RadiationKernelCache: invalid Klag count (nk=" +
                    std::to_string(nk) + ")");
            data.Klag.resize(static_cast<std::size_t>(nk));
            for (auto& K : data.Klag)
                readMatrix(is, K);
        }
        else
        {
            // ---- v7 老格式：等距 K，无 ship/Lpp/times ----
            data.ShipName.clear();
            data.Lpp = 0.0;

            readPod(is, data.Fn);
            readPod(is, data.U);
            readPod(is, data.dt);

            readPod(is, data.TG);
            readPod(is, data.NE);
            readPod(is, data.DOF);

            readIntVector(is, data.modes);

            readMatrix(is, data.A_inf);
            readMatrix(is, data.B);
            readMatrix(is, data.C_prime);
            readMatrix(is, data.K0);

            std::uint64_t nk = 0;
            readPod(is, nk);
            if (!is.good() || nk > kMaxVectorLen)
                throw std::runtime_error(
                    "RadiationKernelCache: invalid Klag count (nk=" +
                    std::to_string(nk) + ")");
            data.Klag.resize(static_cast<std::size_t>(nk));
            for (auto& K : data.Klag)
                readMatrix(is, K);

            // 老格式补 Klag_times = i * dt
            data.Klag_times.resize(data.Klag.size());
            for (std::size_t i = 0; i < data.Klag.size(); ++i)
                data.Klag_times[i] = static_cast<double>(i) * data.dt;
        }
    }
    catch (const std::exception& e)
    {
        // 任何读取异常（损坏文件 / 不匹配格式 / OOM）都视为 cache 不可用，
        // 上游会自动走重建路径，安全保底。
        std::cerr << "[RadiationKernelCache] load failed: " << e.what()
                  << "  (file: " << file << ")\n";
        return false;
    }
    catch (...)
    {
        std::cerr << "[RadiationKernelCache] load failed with unknown exception"
                  << "  (file: " << file << ")\n";
        return false;
    }

    return static_cast<bool>(is);
}
