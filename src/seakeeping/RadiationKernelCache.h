#pragma once

#include "Eigen/Dense"
#include <string>
#include <vector>

struct RadiationKernelData
{
    double Fn = 0.0;
    double U = 0.0;
    double dt = 0.0;     // offline kernel dt（等距时即真实步长；非均匀时仅作"代表性 dt"，存 dt_fine 或 0）

    int TG = 0;          // offline memory steps（= Klag.size()-1，向后兼容字段）
    int NE = 0;
    int DOF = 0;

    std::vector<int> modes;

    // 案例指纹：用于自描述与命名一致性校验（新格式才填）
    std::string ShipName;
    double      Lpp = 0.0;

    // 论文 3.4 / 3.21 / 3.22 对应的刚体系数
    Eigen::MatrixXd A_inf;               // A (常值附加质量项)
    Eigen::MatrixXd B;                   // B (航速附加阻尼常数)
    Eigen::MatrixXd C_prime;             // C' (航速流体恢复力常数)
    std::vector<Eigen::MatrixXd> Klag;   // K(t_m), m = 0..N-1
    std::vector<double>          Klag_times;   // 与 Klag 等长；t[0]=0，严格递增
                                               //（等距时 t[i]=i*dt，非均匀时直接是节点时间）

    // 兼容旧接口：这里存 K(0+)，不再存 A_inf
    Eigen::MatrixXd K0;
};

class RadiationKernelCache
{
public:
    static bool save(const std::string& file, const RadiationKernelData& data);
    static bool load(const std::string& file, RadiationKernelData& data);
};
