#pragma once

#include "ISecondOrderLoadProvider.h"
#include <string>
#include <vector>

class CsvSecondOrderLoadProvider : public ISecondOrderLoadProvider
{
public:
    explicit CsvSecondOrderLoadProvider(const std::string& csvPath);

    CoupledExternalLoads3DOF sample(
        const CoupledSlowState3DOF& s,
        const CoupledEncounterState& enc) const override;

private:
    struct Row
    {
        double Fn = 0.0;
        double betaDeg = 0.0;
        double X2 = 0.0;
        double Y2 = 0.0;
        double N2 = 0.0;
    };

    std::vector<Row> rows_;
};
