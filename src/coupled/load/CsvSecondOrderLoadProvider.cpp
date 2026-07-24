#include "CsvSecondOrderLoadProvider.h"
#include "../../const/Const.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>

CsvSecondOrderLoadProvider::CsvSecondOrderLoadProvider(const std::string& csvPath)
{
    std::ifstream in(csvPath);
    if (!in.is_open())
        throw std::runtime_error("CsvSecondOrderLoadProvider: cannot open file " + csvPath);

    std::string line;
    std::getline(in, line); // header: Fn,betaDeg,X2,Y2,N2
    while (std::getline(in, line))
    {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell;
        Row r{};
        std::getline(ss, cell, ','); r.Fn = std::stod(cell);
        std::getline(ss, cell, ','); r.betaDeg = std::stod(cell);
        std::getline(ss, cell, ','); r.X2 = std::stod(cell);
        std::getline(ss, cell, ','); r.Y2 = std::stod(cell);
        std::getline(ss, cell, ','); r.N2 = std::stod(cell);
        rows_.push_back(r);
    }
}

CoupledExternalLoads3DOF CsvSecondOrderLoadProvider::sample(
    const CoupledSlowState3DOF& s,
    const CoupledEncounterState& enc) const
{
    if (rows_.empty()) return CoupledExternalLoads3DOF{};

    const double betaDeg = enc.betaRel * 180.0 / PI;
    auto best = rows_.front();
    double bestMetric = std::numeric_limits<double>::max();
    for (const auto& r : rows_)
    {
        const double m = std::abs(r.Fn - s.Fn) + 0.05 * std::abs(r.betaDeg - betaDeg);
        if (m < bestMetric)
        {
            bestMetric = m;
            best = r;
        }
    }

    CoupledExternalLoads3DOF out;
    out.X2 = best.X2;
    out.Y2 = best.Y2;
    out.N2 = best.N2;
    out.X = out.X2;
    out.Y = out.Y2;
    out.N = out.N2;
    out.hasSecondOrder = true;
    return out;
}
