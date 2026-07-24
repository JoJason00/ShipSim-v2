#include "../src/seakeeping/Greenf.h"
#include "../src/seakeeping/TDGFProvider.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <tool/findRoot.h>

namespace
{
    constexpr uint32_t MAGIC = 0x54444746; // "TDGF"
    constexpr uint32_t VERSION = 1;

    TDGFValue exactEval(double tau, double mu)
    {
        Greenf gf;
        gf.GreenFunctionCal(tau, mu);

        TDGFValue v;
        v.F = gf.Gf;
        v.Ft = gf.Gbd;
        v.Fm = gf.Gmd;
        return v;
    }

    std::vector<double> buildMuNodes()
    {
        std::vector<double> nodes;

        for (double x = 1.0e-6; x < 0.03; x += 0.001)
            nodes.push_back(x);

        for (double x = 0.03; x < 0.10; x += 0.0025)
            nodes.push_back(x);

        for (double x = 0.10; x < 0.30; x += 0.01)
            nodes.push_back(x);

        for (double x = 0.30; x <= 1.0000001; x += 0.025)
            nodes.push_back(std::min(x, 1.0));

        nodes.push_back(1.0);

        std::sort(nodes.begin(), nodes.end());
        nodes.erase(std::unique(nodes.begin(), nodes.end(),
            [](double a, double b)
        {
            return std::abs(a - b) < 1.0e-12;
        }), nodes.end());

        return nodes;
    }

    std::vector<double> buildTauNodes(double tauMax)
    {
        std::vector<double> nodes;

        for (double x = 1.0e-6; x < 6.0; x += 0.25)
            nodes.push_back(x);

        for (double x = 6.0; x < 12.0; x += 0.20)
            nodes.push_back(x);

        for (double x = 12.0; x < 28.0; x += 0.10)
            nodes.push_back(x);

        for (double x = 28.0; x <= tauMax + 1.0e-12; x += 0.05)
            nodes.push_back(x);

        if (nodes.empty() || nodes.back() < tauMax)
            nodes.push_back(tauMax);

        std::sort(nodes.begin(), nodes.end());
        nodes.erase(std::unique(nodes.begin(), nodes.end(),
            [](double a, double b)
        {
            return std::abs(a - b) < 1.0e-12;
        }), nodes.end());

        return nodes;
    }

    TDGFNode buildStateAt(double tau, double mu)
    {
        const double h = 1.0e-4 * std::max(1.0, tau);

        TDGFValue c = exactEval(tau, mu);

        TDGFNode node;
        node.F = c.F;
        node.Ft = c.Ft;

        if (tau - h <= 1.0e-6)
        {
            TDGFValue p1 = exactEval(tau + h, mu);
            TDGFValue p2 = exactEval(tau + 2.0 * h, mu);

            node.Ftt = (p1.Ft - c.Ft) / h;
            node.Fttt = (p2.Ft - 2.0 * p1.Ft + c.Ft) / (h * h);
        }
        else
        {
            TDGFValue p = exactEval(tau + h, mu);
            TDGFValue m = exactEval(tau - h, mu);

            node.Ftt = (p.Ft - m.Ft) / (2.0 * h);
            node.Fttt = (p.Ft - 2.0 * c.Ft + m.Ft) / (h * h);
        }

        return node;
    }

    void saveTable(
        const std::string& outFile,
        const std::vector<double>& tauNodes,
        const std::vector<double>& muNodes,
        const std::vector<TDGFNode>& table)
    {
        std::ofstream out(outFile, std::ios::binary);
        if (!out)
            throw std::runtime_error("cannot create table file: " + outFile);

        const uint64_t nt = static_cast<uint64_t>(tauNodes.size());
        const uint64_t nm = static_cast<uint64_t>(muNodes.size());

        out.write(reinterpret_cast<const char*>(&MAGIC), sizeof(MAGIC));
        out.write(reinterpret_cast<const char*>(&VERSION), sizeof(VERSION));
        out.write(reinterpret_cast<const char*>(&nt), sizeof(nt));
        out.write(reinterpret_cast<const char*>(&nm), sizeof(nm));

        out.write(reinterpret_cast<const char*>(tauNodes.data()),
            static_cast<std::streamsize>(tauNodes.size() * sizeof(double)));

        out.write(reinterpret_cast<const char*>(muNodes.data()),
            static_cast<std::streamsize>(muNodes.size() * sizeof(double)));

        out.write(reinterpret_cast<const char*>(table.data()),
            static_cast<std::streamsize>(table.size() * sizeof(TDGFNode)));

        if (!out)
            throw std::runtime_error("write table failed: " + outFile);
    }
}

int main(int argc, char** argv)
{
    try
    {
        //=======================
       //  ¼ÓÔØ¹¤¿ö
       //=======================
        std::vector<std::string> cases = { "wigleyI","wigleyI_Lu","S175_1424","S175_1262","S175_924" ,"S175_848", "S175_1722" };
        std::cout << "chose a case to run:\n";
        for (size_t i = 0; i < cases.size(); ++i)
            std::cout << "  " << i << ":\t " << cases[i] << "\n";
        std::cout << "input case number:\t";
        int caseNum;
        std::cin >> caseNum;

        std::string run_case = (argc >= 2) ? argv[1] : cases[caseNum];
        auto root = findProjectRootFromExe(argv[0]);
        std::string casePath = (root / "cases" / run_case).string() + "/";

        std::string outFile = casePath + "tdgf_ode_table.bin";
        double tauMax = 300.0;

        if (argc >= 2)
            outFile = argv[1];

        if (argc >= 3)
            tauMax = std::stod(argv[2]);

        std::vector<double> tauNodes = buildTauNodes(tauMax);
        std::vector<double> muNodes = buildMuNodes();

        std::vector<TDGFNode> table;
        table.resize(tauNodes.size() * muNodes.size());

        const int nt = static_cast<int>(tauNodes.size());
        const int nm = static_cast<int>(muNodes.size());

        std::cout << "Building TDGF ODE table...\n";
        std::cout << "tau nodes: " << nt << "\n";
        std::cout << "mu nodes:  " << nm << "\n";
        std::cout << "tau max:   " << tauMax << "\n";

        for (int it = 0; it < nt; ++it)
        {
            std::cout << "tau node " << (it + 1) << " / " << nt
                << ", tau = " << std::setprecision(10) << tauNodes[it] << "\n";

            for (int im = 0; im < nm; ++im)
            {
                table[static_cast<std::size_t>(it * nm + im)] =
                    buildStateAt(tauNodes[it], muNodes[im]);
            }
        }

        saveTable(outFile, tauNodes, muNodes, table);

        std::cout << "Saved table to: " << outFile << "\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}