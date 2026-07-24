#include "DirectFKPlusDfWaveForceProvider.h"

#include <filesystem>
#include <fstream>
#include <iostream>

DirectFKPlusDfWaveForceProvider::DirectFKPlusDfWaveForceProvider(
    const ShipConfig& ship,
    const SeakeepingConfig& skCfg,
    const std::string& casePath,
    const std::string& shipName,
    const HydroRefreshConfig& refreshCfg,
    std::shared_ptr<Element> element)
    : fkDirect_(ship, skCfg, std::move(element)),
      dfImpulse_(skCfg, casePath, shipName, refreshCfg, FKImpulsePart::DfOnly),
      casePath_(casePath)
{
    std::cout << "[directFKdf] FK = direct pressure integration; "
              << "diffraction = fkImpulse df convolution\n";
}

void DirectFKPlusDfWaveForceProvider::onWindowBegin(
    const CoupledSlowState3DOF& slow0,
    const CoupledEncounterState& enc0,
    double dtFast,
    const CoupledFastWindowInfo& win)
{
    fkDirect_.onWindowBegin(slow0, enc0, dtFast, win);
    dfImpulse_.onWindowBegin(slow0, enc0, dtFast, win);
}

void DirectFKPlusDfWaveForceProvider::setEtaHistory(const std::vector<double>& hist)
{
    fkDirect_.setEtaHistory(hist);
    dfImpulse_.setEtaHistory(hist);
}

std::vector<double> DirectFKPlusDfWaveForceProvider::etaHistory() const
{
    return dfImpulse_.etaHistory();
}

Eigen::VectorXd DirectFKPlusDfWaveForceProvider::evaluate(
    const CoupledSlowState3DOF& slow,
    double t,
    const CoupledEncounterState& enc,
    Eigen::RowVectorXd* force6)
{
    Eigen::RowVectorXd F6fk;
    Eigen::RowVectorXd* pFk = force6 ? &F6fk : nullptr;

    Eigen::VectorXd Ffk = fkDirect_.evaluate(slow, t, enc, pFk);

    Eigen::RowVectorXd F6df;
    Eigen::RowVectorXd* pDf = force6 ? &F6df : nullptr;

    Eigen::VectorXd Fdf = dfImpulse_.evaluate(slow, t, enc, pDf);

    if (force6)
        *force6 = F6fk + F6df;

    // ===== 临时诊断：把 |F_fk|、|F_df|、ratio 写到 coupled_turning/directFKdf_force_dump.csv =====
    // 用完删除/注释这段即可。落点：<casePath>/coupled_turning/directFKdf_force_dump.csv
    {
        static std::ofstream dump([this]() {
            const std::string dir = casePath_ + "coupled_turning/";
            std::filesystem::create_directories(dir);
            return dir + "directFKdf_force_dump.csv";
        }());
        static bool wroteHeader = false;
        if (!wroteHeader)
        {
            dump << "t,|F_fk|,|F_df|,ratio_df_over_fk";
            for (int j = 0; j < Ffk.size(); ++j)
                dump << ",F_fk" << j << ",F_df" << j;
            dump << "\n";
            wroteHeader = true;
        }
        const double nFk = Ffk.norm();
        const double nDf = Fdf.norm();
        dump << t << "," << nFk << "," << nDf << ","
             << (nFk > 1e-30 ? nDf / nFk : 0.0);
        for (int j = 0; j < Ffk.size(); ++j)
            dump << "," << Ffk(j) << "," << Fdf(j);
        dump << "\n";
    }
    // ===== 临时诊断结束 =====

    return Ffk + Fdf;
}
