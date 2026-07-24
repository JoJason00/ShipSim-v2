#pragma once

#include "ISecondOrderLoadProvider.h"
#include <Eigen/Dense>
#include <string>
#include <vector>

class DriftForceTxtProvider : public ISecondOrderLoadProvider
{
public:
    explicit DriftForceTxtProvider(const std::string& folderPath);

    CoupledExternalLoads3DOF sample(
        const CoupledSlowState3DOF& s,
        const CoupledEncounterState& enc) const override;

private:
    struct DriftTableAtSpeed
    {
        double Uref = 0.0; // ���ļ��� Vxx ������������ǰ�ȼٶ��� s.U ͬ��λ
        std::vector<double> betaDegAxis;   // ���� 0,15,...,180
        std::vector<double> omegaAxis;     // 20 �������Ƶ�ʵ�
        std::vector<std::vector<Eigen::Vector4d>> qXYN; // [ibeta][iomega] = (X,Y,N,K)
    };

private:
    static double wrap360(double deg);
    static int lowerBracket(const std::vector<double>& x, double v);
    static double lerp(double a, double b, double t);
    static Eigen::Vector4d lerpVec(const Eigen::Vector4d& a, const Eigen::Vector4d& b, double t);

    void loadFolder(const std::string& folderPath);
    void loadOneFile(const std::string& filePath, double Uref);
    Eigen::Vector4d evalAtTable(const DriftTableAtSpeed& tab, double betaDeg, double omega) const;
    const DriftTableAtSpeed& pickNearestSpeedTable(double U) const;
    // Linearly interpolate q(β, ω) over the two speed tables bracketing U.
    // Falls back to evalAtTable() when U is outside the loaded speed range
    // or only one table is loaded.
    Eigen::Vector4d evalAtSpeed(double U, double betaDeg, double omega) const;

private:
    std::vector<DriftTableAtSpeed> tables_;
};