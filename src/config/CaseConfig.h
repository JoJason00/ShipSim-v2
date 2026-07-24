#pragma once

#include "WaveConfig.h"
#include "SeakeepingConfig.h"
#include "MmgConfig.h"
#include "CoupledConfig.h"
#include <vector>
#include <optional>

struct ShipGeometry
{
    double Length{};
    double Draft{};
    double Displacement{};
    double Breadth{};
    double CB{};

    // ������ˮǰ��һ��ʱ����ˮ�Ĳ�ֵΪTrim��Trim>0��ʾ���׳�ˮ���Trim<0��ʾ��β��ˮ����
    double Trim{};
};

struct ShipMass
{
    double Mass{};
    double Ixx{}, Iyy{}, Izz{};
    std::vector<double> CG;
    double GM{};
};

struct ShipConfig
{
    std::string Name;
    ShipGeometry Geometry;
    ShipMass Mass;
};



struct CaseConfig
{
    ShipConfig Ship;
    bool enable_seakeeping{};
    bool enable_maneuvering{};
    bool enable_coupling{};
    std::optional<SeakeepingConfig> Seakeeping;
    std::optional<MmgConfig> Mmg;
    std::optional<CoupledConfig> Coupled;
};