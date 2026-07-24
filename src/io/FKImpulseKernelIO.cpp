#include "FKImpulseKernelIO.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>

namespace FKImpulseKernelIO
{
    std::string keyDouble(double x)
    {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(4) << x;

        std::string s = ss.str();
        for (char& c : s)
        {
            if (c == '.') c = 'p';
            if (c == '-') c = 'm';
        }

        return s;
    }

    std::string makePath(const Params& p)
    {
        std::filesystem::path dir =
            std::filesystem::path(p.casePath) / "fkImpulse";

        std::ostringstream name;
        name << "fkImpulse_";

        if (p.includeShipName && !p.shipName.empty())
        {
            name << p.shipName << "_";
        }

        name << "Fn"   << keyDouble(p.Fn)
             << "_Dir" << keyDouble(p.dirRad)
             << "_Reg" << wave_force_region::impulseBucketTag(p.region)
             << "_dt"  << keyDouble(p.dt)
             << "_tMot" << p.tMot
             << ".csv";

        return (dir / name.str()).string();
    }

    void save(const Params& p, const Data& data)
    {
        const std::string outFile = makePath(p);

        std::filesystem::create_directories(
            std::filesystem::path(outFile).parent_path());

        const int nImp = 2 * p.tMot + 1;

        if (data.fkForce.rows() < nImp || data.fkForce.cols() < 6)
            throw std::runtime_error("FKImpulseKernelIO::save: invalid fkForce size.");

        if (data.dForce.rows() < nImp || data.dForce.cols() < 6)
            throw std::runtime_error("FKImpulseKernelIO::save: invalid dForce size.");

        std::ofstream out(outFile);
        if (!out.is_open())
            throw std::runtime_error("FKImpulseKernelIO::save: cannot create " + outFile);

        out << std::setprecision(17);

        out << "# FK impulse kernel cache\n";
        out << "# Fn,"     << p.Fn << "\n";
        out << "# Dir,"    << p.dirRad << "\n";
        out << "# Region," << wave_force_region::impulseBucketTag(p.region) << "\n";
        out << "# dt,"     << p.dt << "\n";
        out << "# tMot,"   << p.tMot << "\n";
        out << "# rows,"   << nImp << "\n";

        out << "i,tau";

        for (int j = 0; j < 6; ++j)
            out << ",fk" << j;

        for (int j = 0; j < 6; ++j)
            out << ",df" << j;

        for (int j = 0; j < 6; ++j)
            out << ",total" << j;

        out << "\n";

        for (int i = 0; i < nImp; ++i)
        {
            const double tau = (i - p.tMot) * p.dt;

            out << i << "," << tau;

            for (int j = 0; j < 6; ++j)
                out << "," << data.fkForce(i, j);

            for (int j = 0; j < 6; ++j)
                out << "," << data.dForce(i, j);

            for (int j = 0; j < 6; ++j)
                out << "," << data.fkForce(i, j) + data.dForce(i, j);

            out << "\n";
        }

        std::cout << "[FKImpulse] saved: " << outFile << "\n";
    }

    bool load(const Params& p, Data& data)
    {
        // Exact-match path uses the canonical filename built from p.
        LoadInfo info{};
        return loadFromFile(makePath(p), info, data);
    }

    bool loadFromFile(const std::string& inFile, LoadInfo& info, Data& data)
    {
        std::ifstream in(inFile);
        if (!in.is_open())
            return false;

        int tMot = -1;
        int rows = -1;
        double Fn = 0.0, dirRad = 0.0, dt = 0.0;
        WaveForceRegion region = WaveForceRegion::Head;

        std::string line;
        // Parse header lines: `# key,value` until the column header `i,tau,...`.
        while (std::getline(in, line))
        {
            if (line.empty()) continue;
            if (line[0] == '#')
            {
                const auto comma = line.find(',');
                if (comma == std::string::npos) continue;
                int keyStart = 1;
                while (keyStart < static_cast<int>(comma)
                       && std::isspace(static_cast<unsigned char>(line[keyStart])))
                    ++keyStart;
                const std::string key = line.substr(keyStart, comma - keyStart);
                const std::string val = line.substr(comma + 1);
                if      (key == "Fn")     Fn     = std::stod(val);
                else if (key == "Dir")    dirRad = std::stod(val);
                else if (key == "Region") region = wave_force_region::fromTag(val);
                else if (key == "dt")     dt     = std::stod(val);
                else if (key == "tMot")   tMot   = std::stoi(val);
                else if (key == "rows")   rows   = std::stoi(val);
                continue;
            }
            if (line.rfind("i,tau", 0) == 0) break;
        }

        if (tMot <= 0 || rows != 2 * tMot + 1 || dt <= 0.0)
        {
            std::cout << "[FKImpulse] invalid header in " << inFile << "\n";
            return false;
        }

        const int nImp = 2 * tMot + 1;
        Eigen::MatrixXd fk = Eigen::MatrixXd::Zero(nImp, 6);
        Eigen::MatrixXd df = Eigen::MatrixXd::Zero(nImp, 6);
        int loadedRows = 0;

        while (std::getline(in, line))
        {
            if (line.empty()) continue;

            std::stringstream ss(line);
            std::string cell;
            std::vector<double> values;
            while (std::getline(ss, cell, ','))
                if (!cell.empty()) values.push_back(std::stod(cell));

            // i,tau,fk0..fk5,df0..df5,total0..total5
            if (values.size() < 14)
                throw std::runtime_error("FKImpulseKernelIO::loadFromFile: invalid row in " + inFile);

            const int i = static_cast<int>(values[0]);
            if (i < 0 || i >= nImp) continue;

            for (int j = 0; j < 6; ++j)
            {
                fk(i, j) = values[2 + j];
                df(i, j) = values[8 + j];
            }
            ++loadedRows;
        }

        if (loadedRows != nImp)
        {
            std::cout << "[FKImpulse] cache row mismatch, ignore: " << inFile << "\n";
            return false;
        }

        info.Fn     = Fn;
        info.dirRad = dirRad;
        info.region = region;
        info.dt     = dt;
        info.tMot   = tMot;
        data.fkForce = std::move(fk);
        data.dForce  = std::move(df);
        return true;
    }

    std::string findByKey(const std::string& casePath,
                          const std::string& shipName,
                          double Fn,
                          double dirRad,
                          WaveForceRegion region)
    {
        const std::filesystem::path dir =
            std::filesystem::path(casePath) / "fkImpulse";
        if (!std::filesystem::exists(dir)) return {};

        std::string  bestPath;
        double       bestMemory = -1.0;            // longest dt·tMot wins on ties

        for (const auto& de : std::filesystem::directory_iterator(dir))
        {
            if (!de.is_regular_file())                          continue;
            if (de.path().extension() != ".csv")                continue;
            const std::string fn = de.path().filename().string();
            if (fn.rfind("fkImpulse_", 0) != 0)                 continue;

            std::ifstream in(de.path());
            if (!in.is_open()) continue;

            // Header-only sniff (no body parse) for speed.
            int tMot = -1;
            double thisFn = 0.0, thisDir = 0.0, thisDt = 0.0;
            WaveForceRegion thisReg = WaveForceRegion::Head;
            std::string line;
            while (std::getline(in, line))
            {
                if (line.empty()) continue;
                if (line[0] != '#') break;
                const auto comma = line.find(',');
                if (comma == std::string::npos) continue;
                int ks = 1;
                while (ks < static_cast<int>(comma)
                       && std::isspace(static_cast<unsigned char>(line[ks]))) ++ks;
                const std::string key = line.substr(ks, comma - ks);
                const std::string val = line.substr(comma + 1);
                if      (key == "Fn")     thisFn  = std::stod(val);
                else if (key == "Dir")    thisDir = std::stod(val);
                else if (key == "Region") thisReg = wave_force_region::fromTag(val);
                else if (key == "dt")     thisDt  = std::stod(val);
                else if (key == "tMot")   tMot    = std::stoi(val);
            }

            const bool shipOk = shipName.empty() || fn.find(shipName) != std::string::npos;
            if (!shipOk)                                        continue;
            if (std::abs(thisFn  - Fn)     > 1e-6)              continue;
            if (std::abs(thisDir - dirRad) > 1e-4)              continue;
            if (thisReg != region)                              continue;
            if (tMot <= 0 || thisDt <= 0.0)                     continue;

            const double mem = thisDt * tMot;
            if (mem > bestMemory)
            {
                bestMemory = mem;
                bestPath   = de.path().string();
            }
        }
        return bestPath;
    }

    void resample(const Data& src, double srcDt, int srcTMot,
                  double dstDt, int dstTMot, Data& dst)
    {
        const int nDst = 2 * dstTMot + 1;
        dst.fkForce = Eigen::MatrixXd::Zero(nDst, 6);
        dst.dForce  = Eigen::MatrixXd::Zero(nDst, 6);

        const int srcRows = 2 * srcTMot + 1;
        if (srcRows <= 0 || srcDt <= 0.0) return;

        for (int idst = 0; idst < nDst; ++idst)
        {
            const double tau = (idst - dstTMot) * dstDt;        // dst τ
            // Convert to fractional row index in source.
            double x = tau / srcDt + srcTMot;
            if (x <= 0.0)            x = 0.0;
            if (x >= srcRows - 1)    x = srcRows - 1;
            const int    i0 = static_cast<int>(std::floor(x));
            const int    i1 = std::min(srcRows - 1, i0 + 1);
            const double a  = x - i0;
            for (int j = 0; j < 6; ++j)
            {
                dst.fkForce(idst, j) = (1.0 - a) * src.fkForce(i0, j) + a * src.fkForce(i1, j);
                dst.dForce(idst,  j) = (1.0 - a) * src.dForce(i0,  j) + a * src.dForce(i1,  j);
            }
        }
    }
}