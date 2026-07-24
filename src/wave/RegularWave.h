#pragma once

#include "WaveBase.h"
#include "../config/WaveConfig.h"
#include "Encounter.h"

class RegularWave : public WaveBase
{
public:
    RegularWave(const RegularWaveConfig& regularwave);

    void loadData(const fkpData& Data) override;
    double Eta(double t) const override;
    void Exciting(double tn, FKphi& fkphi) override;

    double direction() override;
    double getAmp() override;
    double getFreq() override;
    double initialPhase();   // config.phase0 [rad], default 0
    
private:
    double w = 0.0;
    double we = 0.0;
    EncounterInfo encounter_;
    const RegularWaveConfig config;

    std::array<std::array<double, 2>, 10> wt;

    void Headseas (double tn, int i, double& eVx, double& eVy, double& eVz, double& ePr);
    void FollowingSeas(double tn, int i, double& eVx, double& eVy, double& eVz, double& ePr);
    void Cerrorfun(double& x, double& y, std::array<double,2>& Riw);
    void Taylor   (std::array<double, 2>& Riw, double& x, double& y);
    void Gaussh   (std::array<double, 2>& Riw, double& x, double& y);
    void Quadrant (std::array<double, 2>& Riw, double& x, double& y, double& sgnx, double& sgny);
    void Quadrant1(double k, double& x, double& y, double& Rew1, double& Imw1, double& Rew, double& Imw);
};
