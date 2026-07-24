#pragma once
#include "../common/CoupledTypes.h"
#include "../environment/WaveEnvironment.h"

// 快时间耐波统一接口
class IWindowSeakeepingSolver
{
public:
    virtual ~IWindowSeakeepingSolver() = default;

    virtual CoupledWindowResult solveWindow(
        const CoupledWindowRequest& req,
        const CoupledWaveEnvironment& env) = 0;

    double Fn_now{ 0.0 };
};