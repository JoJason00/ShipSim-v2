//#pragma once
//
//#include "../config/CaseConfig.h"
//#include <string>
//#include <vector>
//#include <memory>
//#include "Eigen/Dense"
//#include <wave/WaveBase.h>
//#include <wave/RegularWave.h>
//
//using Vec6 = Eigen::Matrix<double, 6, 1>;
//using Mat6 = Eigen::Matrix<double, 6, 6>;
//using Vec3 = Eigen::Vector3d;
//using Mat3 = Eigen::Matrix<double, 3, 3>;
//
//struct CaseContext_mmg {
//    int i_wave{};
//
//    // 为了不破坏后续波浪接口，先保留
//    std::shared_ptr<WaveBase> wave;
//    std::shared_ptr<RegularWave> reg;
//
//    double dt{};
//    double Fn{};
//    double U{};
//};
//
//struct MmgMassND {
//    // Fortran 里的:
//    // Mm(1)=M0(1)/(0.5*rou*Lpp*Lpp*drau)
//    // Mm(4)=M0(4)/(0.5*rou*drau*Lpp^4)
//    // m_x, m_y, J_zz 也做相同无量纲化
//    double Mm1{};
//    double Mm4{};
//    double mx{};
//    double my{};
//    double Jz{};
//};
//
//struct HullForceND {
//    double X{};
//    double Y{};
//    double N{};
//};
//
//struct RudderForceDim {
//    double X{};
//    double Y{};
//    double N{};
//};
//
//class Mmg
//{
//public:
//    Mmg(const ShipConfig& Ship, std::string casePath, const MmgConfig& Mmg);
//    void run();
//    void runTurningCase(const CaseContext_mmg& ctx, const bool direction);
//    void runZigzagCase(const CaseContext_mmg& ctx);
//
//private:
//    void compute_hullCoeff();
//    void compute_propellerCoeff();
//    void compute_rudderCoeff();
//    void compute_addedCoeff();
//    void compute_normMass();
//
//    CaseContext_mmg buildCaseContext(const int i_wave, const double Fn);
//
//    void stepRK4Turning(Vec6& Y, double t, double dt, const CaseContext_mmg& ctx, const double rudder_angle);
//    Vec6 stateDerivativeTurning(double t, const Vec6& Y, const CaseContext_mmg& ctx, const double rudder_angle);
//
//    // 按 Fortran 静水操纵逻辑拆分
//    HullForceND HullForceFortranND(const Vec6& Y) const;
//    double ResistanceForce(const Vec6& Y) const;
//    double PropellerForce(const Vec6& Y, double& wP, double& Jp, double& KT) const;
//    RudderForceDim RudderForce(const Vec6& Y, const double rudder_angle, const double Jp, const double KT) const;
//
//private:
//    std::string         filePath;
//    ShipConfig          ShipCfg;
//
//    // 先保留波浪接口，静水验证阶段不使用波浪力
//    std::vector<std::shared_ptr<WaveBase>> waves;
//
//    MmgTimeConfig Time;
//    hullCoeff Hull;
//    propellerCoeff Propeller;
//    rudderCoeff Rudder;
//    addedCoeff Added;
//
//    MmgMassND MassND;
//
//    std::vector<double> Fns;
//
//    TuringCfg Turning;
//    ZigzagCfg Zigzag;
//
//    bool turn_on_turningcase{};
//    bool turn_on_zigzagcase{};
//};