#pragma once

#include "IWindowSeakeepingSolver.h"
#include "../hydro/CoupledRadiationKernelRepo.h"
#include "../../config/CaseConfig.h"
#include "../../config/SeakeepingConfig.h"
#include "../../const/Const.h"
#include "../../seakeeping/RollDamping.h"
#include "../../seakeeping/SeakeepingDOF.h"

#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <array>
#include <seakeeping/LinearCumminsTDGF.h>

#include "../../io/CaseLoader.h"
#include "../../seakeeping/Element.h"

#include <filesystem>
#include <iostream>

namespace {
    static double clamp01(double x)
    {
        return std::max(0.0, std::min(1.0, x));
    }
}


class FittedWaveForceWindowSolver : public IWindowSeakeepingSolver
{
public:
    FittedWaveForceWindowSolver(
        const ShipConfig& ship,
        const SeakeepingConfig& skCfg,
        const std::string& casePath,
        const CoupledRadiationKernelRepo& radRepo,
        const std::string& forceTableCsv)
        : ship_(ship),
        skCfg_(skCfg),
        casePath_(casePath),
        radRepo_(radRepo)
    {
        // �ؼ������ϴ��Ͳ���ʽ�ֳ���ˮ������������ defaultHydrostatics
        loadHydrostaticsByOldMethod();

        // ��ȡ��ֵ��λ���ݿ�
        loadForceTable(forceTableCsv);

        // ��ҡճ�����ᣬ���ϴ��Ͳ�����һ��
        rollVisc_ = RollDampingBuilder::build(
            casePath_,
            skCfg_.RollDamping,
            ship_.Mass.Mass,
            ship_.Mass.GM);

        std::cout << "[FittedWaveForceWindowSolver] legacy seakeeping mode\n";
        std::cout << "  c33 = " << hs_.c33 << "\n";
        std::cout << "  c44 = " << hs_.c44 << "\n";
        std::cout << "  c55 = " << hs_.c55 << "\n";
        std::cout << "  c35 = " << hs_.c35 << "\n";
    }

    CoupledWindowResult solveWindow(
        const CoupledWindowRequest& req,
        const CoupledWaveEnvironment& env) override
    {
        CoupledWindowResult out;

        const int D = skCfg_.DOF;
        const int N = std::max(1, req.nFast);
        const double dt = req.dtFast;

        out.loads = CoupledWaveLoadsWindow{};

        out.tHist = Eigen::VectorXd::Zero(N);
        out.qHist = Eigen::MatrixXd::Zero(N, D);
        out.vHist = Eigen::MatrixXd::Zero(N, D);

        // �����֮ǰ�Ѿ��� CoupledWindowResult �������Щ�����ֶΣ��ͱ���
        out.forceHist = Eigen::MatrixXd::Zero(N, D);
        out.force6Hist = Eigen::MatrixXd::Zero(N, 6);
        out.betaRelDegHist = Eigen::VectorXd::Zero(N);
        out.betaQueryDegHist = Eigen::VectorXd::Zero(N);
        out.mirroredHist = Eigen::VectorXd::Zero(N);
        out.thetaHist = Eigen::VectorXd::Zero(N);
        out.omegaEncounterHist = Eigen::VectorXd::Zero(N);
        out.waveAmpHist = Eigen::VectorXd::Zero(N);
        out.etaCGHist = Eigen::VectorXd::Zero(N);

        // 1. ��ȡ/��ֵ����ˣ�radRepo �ڲ��Ѿ��� dtFast �ز���
        const double Fn0 = currentFnFromSpeed(req.slow0.U);
        const double Fn1 = currentFnFromSpeed(req.slow1_pred.U);

        // �����ֻ�溽��/Fn �л������������л���
        // һ���������ڶ���Ϊ����ƽ���ٶȶ�Ӧ�ĺˡ�
        const double FnRad = 0.5 * (Fn0 + Fn1);

        const RadiationKernelData rad =
            radRepo_.getKernel(FnRad, dt);

        // 2. ���ϴ��Ͳ���ʽ��װ M/B/C
        Eigen::MatrixXd Mphys, Bphys, Cphys;
        SeakeepingDOF::buildSystemMatrices(
            ship_,
            skCfg_,
            hs_,
            Mphys,
            Bphys,
            Cphys);

        // 3. �� Cummins ����
        const Eigen::MatrixXd Mtotal = Mphys + rad.A_inf;
        const Eigen::MatrixXd Bconst = Bphys + rad.B;
        const Eigen::MatrixXd Ctotal = Cphys + rad.C_prime;

        const double beta = 0.25;
        const double gamma = 0.5;

        const Eigen::MatrixXd Acoef =
            Mtotal
            + gamma * dt * Bconst
            + beta * dt * dt * Ctotal;

        Eigen::FullPivLU<Eigen::MatrixXd> solver(Acoef);
        if (!solver.isInvertible())
            throw std::runtime_error("FittedWaveForceWindowSolver: Acoef is singular.");

        Eigen::FullPivLU<Eigen::MatrixXd> massSolver(Mtotal);
        if (!massSolver.isInvertible())
            throw std::runtime_error("FittedWaveForceWindowSolver: Mtotal is singular.");

        const int iRoll = SeakeepingDOF::findModeIndex(skCfg_, 3);

        Eigen::VectorXd q = Eigen::VectorXd::Zero(D);
        Eigen::VectorXd v = Eigen::VectorXd::Zero(D);
        Eigen::VectorXd a = Eigen::VectorXd::Zero(D);

        if (req.fast0.initialized)
        {
            q = req.fast0.q;
            v = req.fast0.v;
            a = req.fast0.a;

            resizeOrZero(q, D);
            resizeOrZero(v, D);
            resizeOrZero(a, D);
        }

        // 4. �細�ڱ����ٶ���ʷ������ Klag ����������
        std::vector<Eigen::VectorXd> vMemHist = req.fast0.vMemHist;

        const int maxMem = std::max(1, static_cast<int>(rad.Klag.size()) - 1);
        trimVelocityHistory(vMemHist, maxMem, D);

        // 5. ��һ�����������ʼ���ٶ�
        if (!req.fast0.initialized)
        {
            const ForceEval fe0 = fittedForceAt(req.slow0, req.t0, env);

            Eigen::VectorXd rhs0 =
                fe0.Fmode - Bconst * v - Ctotal * q;

            if (iRoll >= 0 && iRoll < D)
                rhs0(iRoll) += rollVisc_.moment(v(iRoll));

            a = massSolver.solve(rhs0);
        }

        // 6. ֻ��ӡһ�ξ��󣬼�������ź�����
        static bool printed = false;
        if (!printed)
        {
            printed = true;

            std::cout << "\n[FittedWaveForceWindowSolver Matrix]\n";
            std::cout << "Mphys =\n" << Mphys << "\n";
            std::cout << "A_inf =\n" << rad.A_inf << "\n";
            std::cout << "Mtotal =\n" << Mtotal << "\n";

            std::cout << "Bphys =\n" << Bphys << "\n";
            std::cout << "rad.B =\n" << rad.B << "\n";
            std::cout << "Bconst =\n" << Bconst << "\n";

            std::cout << "Cphys =\n" << Cphys << "\n";
            std::cout << "C_prime =\n" << rad.C_prime << "\n";
            std::cout << "Ctotal =\n" << Ctotal << "\n";

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esM(
                0.5 * (Mtotal + Mtotal.transpose()));
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esB(
                0.5 * (Bconst + Bconst.transpose()));
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esC(
                0.5 * (Ctotal + Ctotal.transpose()));

            std::cout << "eig(Mtotal) = "
                << esM.eigenvalues().transpose() << "\n";
            std::cout << "eig(Bconst) = "
                << esB.eigenvalues().transpose() << "\n";
            std::cout << "eig(Ctotal) = "
                << esC.eigenvalues().transpose() << "\n\n";
        }

        // 7. Newmark ���֣��� Klag �������
        for (int n = 0; n < N; ++n)
        {
            const double alpha =
                static_cast<double>(n + 1) / static_cast<double>(N);

            const CoupledSlowState3DOF sSlow =
                interpSlow(req.slow0, req.slow1_pred, alpha);

            const double t = req.t0 + (n + 1) * dt;

            const ForceEval fe = fittedForceAt(sSlow, t, env);
            const Eigen::VectorXd Fnow = fe.Fmode;

            if (Fnow.size() != D)
                throw std::runtime_error("FittedWaveForceWindowSolver: Fnow size mismatch.");

            const Eigen::VectorXd qPred =
                q + dt * v + dt * dt * (0.5 - beta) * a;

            const Eigen::VectorXd vPred =
                v + dt * (1.0 - gamma) * a;

            // ���������sum dt*Klag[lag]*v(t-lag)
            Eigen::VectorXd mem = Eigen::VectorXd::Zero(D);

            const int available = static_cast<int>(vMemHist.size());
            const int maxLag = std::min(
                std::min(available, maxMem),
                static_cast<int>(rad.Klag.size()) - 1);

            for (int lag = 1; lag <= maxLag; ++lag)
            {
                const int idx = available - lag;
                if (idx < 0)
                    break;

                mem.noalias() += dt
                    * rad.Klag[static_cast<std::size_t>(lag)]
                    * vMemHist[static_cast<std::size_t>(idx)];
            }

            Eigen::VectorXd rhs =
                Fnow
                - Bconst * vPred
                - Ctotal * qPred
                - mem;

            if (iRoll >= 0 && iRoll < D)
                rhs(iRoll) += rollVisc_.moment(vPred(iRoll));

            const Eigen::VectorXd aNew = solver.solve(rhs);

            const Eigen::VectorXd qNew =
                qPred + beta * dt * dt * aNew;

            const Eigen::VectorXd vNew =
                vPred + gamma * dt * aNew;

            q = qNew;
            v = vNew;
            a = aNew;

            vMemHist.push_back(v);
            trimVelocityHistory(vMemHist, maxMem, D);

            out.tHist(n) = t;
            out.qHist.row(n) = q.transpose();
            out.vHist.row(n) = v.transpose();

            out.forceHist.row(n) = fe.Fmode.transpose();
            out.force6Hist.row(n) = fe.F6.transpose();

            out.betaRelDegHist(n) = fe.betaRelDeg;
            out.betaQueryDegHist(n) = fe.betaQueryDeg;
            out.mirroredHist(n) = fe.mirrored;
            out.thetaHist(n) = fe.theta;
            out.omegaEncounterHist(n) = fe.omegaEncounter;
            out.waveAmpHist(n) = fe.waveAmp;
            out.etaCGHist(n) = fe.etaCG;
        }

        out.etaMean = out.qHist.colwise().mean().transpose();

        out.etaRms = Eigen::VectorXd::Zero(D);
        for (int j = 0; j < D; ++j)
        {
            const Eigen::ArrayXd x = out.qHist.col(j).array();
            out.etaRms(j) = std::sqrt((x * x).mean());
        }

        out.fastLast.q = q;
        out.fastLast.v = v;
        out.fastLast.a = a;
        out.fastLast.initialized = true;

        // �ؼ������ٶ���ʷ������һ���촰
        out.fastLast.vMemHist = vMemHist;

        // ���ڲ���ʵʱ�������������� etaHist ���
        out.fastLast.etaHist.clear();

        return out;
    }

private:
    struct ForceRow
    {
        double Fn = 0.0;
        double betaDeg = 0.0;          // ���ݿ�ֻ�� 0~180
        double omegaIncident = 0.0;
        double omegaEncounter = 0.0;
        double waveAmp = 1.0;

        std::array<double, 6> A{};
        std::array<double, 6> P{};
    };

    struct ForceEval
    {
        Eigen::VectorXd Fmode;   // DOF ˳�򣬶�Ӧ skCfg_.modes
        Eigen::VectorXd F6;      // 6DOF ˳��F0~F5

        double betaRelDeg = 0.0;
        double betaQueryDeg = 0.0;
        double mirrored = 0.0;
        double theta = 0.0;
        double omegaEncounter = 0.0;
        double waveAmp = 0.0;
        double etaCG = 0.0;
    };

private:
    ShipConfig ship_;
    SeakeepingConfig skCfg_;
    std::string casePath_;
    const CoupledRadiationKernelRepo& radRepo_;

    // ���������ϴ��Ͳ���ʽ�ֳ�����õ���ˮ����
    HydrostaticsData hs_;

    // ��������ҡճ�����ᣬ���ϴ��Ͳ�����һ��
    RollViscDamping rollVisc_;

    std::vector<ForceRow> rows_;

private:
    double nearestOmegaInTable(double omegaIncident) const
    {
        if (rows_.empty())
            throw std::runtime_error("FittedWaveForceWindowSolver: empty force table.");

        double bestOmega = rows_.front().omegaIncident;
        double bestErr = std::abs(bestOmega - omegaIncident);

        for (const auto& r : rows_)
        {
            const double e = std::abs(r.omegaIncident - omegaIncident);
            if (e < bestErr)
            {
                bestErr = e;
                bestOmega = r.omegaIncident;
            }
        }

        return bestOmega;
    }

    std::vector<double> collectFnAxisAtOmega(double omegaIncident) const
    {
        const double omegaUse = nearestOmegaInTable(omegaIncident);
        const double omegaTol = 1.0e-3;

        std::vector<double> fns;

        for (const auto& r : rows_)
        {
            if (std::abs(r.omegaIncident - omegaUse) < omegaTol)
                fns.push_back(r.Fn);
        }

        std::sort(fns.begin(), fns.end());
        fns.erase(std::unique(fns.begin(), fns.end(),
            [](double a, double b)
        {
            return std::abs(a - b) < 1.0e-10;
        }),
            fns.end());

        if (fns.empty())
            throw std::runtime_error("FittedWaveForceWindowSolver: no Fn axis found.");

        return fns;
    }

    std::complex<double> interpBetaAtFixedFnOmega(
        double Fn,
        double omegaIncident,
        double betaDeg,
        int mode) const
    {
        const double omegaUse = nearestOmegaInTable(omegaIncident);
        const double omegaTol = 1.0e-3;

        std::vector<ForceRow> cand;

        for (const auto& r : rows_)
        {
            if (std::abs(r.Fn - Fn) < 1.0e-8 &&
                std::abs(r.omegaIncident - omegaUse) < omegaTol)
            {
                cand.push_back(r);
            }
        }

        if (cand.empty())
        {
            throw std::runtime_error(
                "FittedWaveForceWindowSolver: no rows at fixed Fn/Omega.");
        }

        std::sort(cand.begin(), cand.end(),
            [](const ForceRow& a, const ForceRow& b)
        {
            return a.betaDeg < b.betaDeg;
        });

        if (cand.size() == 1)
            return complexAmpPhase(cand.front(), mode);

        if (betaDeg <= cand.front().betaDeg)
            return complexAmpPhase(cand.front(), mode);

        if (betaDeg >= cand.back().betaDeg)
            return complexAmpPhase(cand.back(), mode);

        for (std::size_t i = 0; i + 1 < cand.size(); ++i)
        {
            const ForceRow& r0 = cand[i];
            const ForceRow& r1 = cand[i + 1];

            if (betaDeg >= r0.betaDeg && betaDeg <= r1.betaDeg)
            {
                const double den = r1.betaDeg - r0.betaDeg;
                const double tb =
                    (std::abs(den) > 1.0e-12)
                    ? (betaDeg - r0.betaDeg) / den
                    : 0.0;

                const std::complex<double> z0 = complexAmpPhase(r0, mode);
                const std::complex<double> z1 = complexAmpPhase(r1, mode);

                return (1.0 - clamp01(tb)) * z0 + clamp01(tb) * z1;
            }
        }

        return complexAmpPhase(cand.back(), mode);
    }

    std::complex<double> interpForceComplex(
        double FnNow,
        double omegaIncident,
        double betaDeg,
        int mode) const
    {
        const std::vector<double> fnAxis =
            collectFnAxisAtOmega(omegaIncident);

        double Fn0 = 0.0;
        double Fn1 = 0.0;
        double tf = 0.0;

        bracketAxis(fnAxis, FnNow, Fn0, Fn1, tf);

        const std::complex<double> z0 =
            interpBetaAtFixedFnOmega(Fn0, omegaIncident, betaDeg, mode);

        const std::complex<double> z1 =
            interpBetaAtFixedFnOmega(Fn1, omegaIncident, betaDeg, mode);

        return (1.0 - tf) * z0 + tf * z1;
    }

    double currentFnFromSpeed(double U) const
    {
        const double L = ship_.Geometry.Length;
        if (L <= 0.0)
            return 0.0;

        return U / std::sqrt(G * L);
    }

    static void resizeOrZero(Eigen::VectorXd& x, int n)
    {
        if (x.size() != n)
            x = Eigen::VectorXd::Zero(n);
    }

    static double wrap360(double deg)
    {
        while (deg < 0.0) deg += 360.0;
        while (deg >= 360.0) deg -= 360.0;
        return deg;
    }

    static void bracketAxis(
        const std::vector<double>& axis,
        double x,
        double& x0,
        double& x1,
        double& tx)
    {
        if (axis.empty())
            throw std::runtime_error("bracketAxis: empty axis.");

        if (axis.size() == 1)
        {
            x0 = axis.front();
            x1 = axis.front();
            tx = 0.0;
            return;
        }

        if (x <= axis.front())
        {
            x0 = axis.front();
            x1 = axis.front();
            tx = 0.0;
            return;
        }

        if (x >= axis.back())
        {
            x0 = axis.back();
            x1 = axis.back();
            tx = 0.0;
            return;
        }

        auto it = std::upper_bound(axis.begin(), axis.end(), x);
        const int i1 = static_cast<int>(std::distance(axis.begin(), it));
        const int i0 = i1 - 1;

        x0 = axis[static_cast<std::size_t>(i0)];
        x1 = axis[static_cast<std::size_t>(i1)];

        const double den = x1 - x0;
        tx = (std::abs(den) > 1.0e-12) ? (x - x0) / den : 0.0;
        tx = clamp01(tx);
    }

    static int mirrorSignOfMode(int mode)
    {
        // ���ҶԳƴ��壺
        // surge/heave/pitch ͬ�ţ�sway/roll/yaw ����
        if (mode == 1 || mode == 3 || mode == 5)
            return -1;
        return +1;
    }

    static CoupledSlowState3DOF interpSlow(
        const CoupledSlowState3DOF& a,
        const CoupledSlowState3DOF& b,
        double alpha)
    {
        CoupledSlowState3DOF s;
        const double om = 1.0 - alpha;

        s.t = om * a.t + alpha * b.t;
        s.u = om * a.u + alpha * b.u;
        s.v = om * a.v + alpha * b.v;
        s.r = om * a.r + alpha * b.r;
        s.xe = om * a.xe + alpha * b.xe;
        s.ye = om * a.ye + alpha * b.ye;
        s.psi = om * a.psi + alpha * b.psi;
        s.delta = om * a.delta + alpha * b.delta;
        s.U = om * a.U + alpha * b.U;
        s.betaDrift = om * a.betaDrift + alpha * b.betaDrift;
        s.Fn = om * a.Fn + alpha * b.Fn;
        s.rDotRadPerS2 = om * a.rDotRadPerS2 + alpha * b.rDotRadPerS2;

        return s;
    }


    static std::vector<std::string> splitCsv(const std::string& line)
    {
        std::vector<std::string> cells;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ','))
            cells.push_back(cell);

        return cells;
    }

    static double snapBetaDeg(double beta)
    {
        const double targets[] = { 0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0 };

        for (double b : targets)
        {
            if (std::abs(beta - b) < 0.5)
                return b;
        }

        return beta;
    }

    void loadForceTable(const std::string& csv)
    {
        std::ifstream in(csv);
        if (!in.is_open())
            throw std::runtime_error("FittedWaveForceWindowSolver: cannot open " + csv);

        rows_.clear();

        std::string line;
        std::getline(in, line); // header

        while (std::getline(in, line))
        {
            if (line.empty())
                continue;

            const auto c = splitCsv(line);

            // Fn,betaDeg,omegaIncident,omegaEncounter,waveAmp,A0,P0,...,A5,P5
            if (c.size() < 17)
                continue;

            ForceRow r;
            r.Fn = std::stod(c[0]);
            r.betaDeg = snapBetaDeg(std::stod(c[1]));
            r.omegaIncident = std::stod(c[2]);
            r.omegaEncounter = std::stod(c[3]);
            r.waveAmp = std::stod(c[4]);

            int k = 5;
            for (int mode = 0; mode < 6; ++mode)
            {
                r.A[mode] = std::stod(c[k++]);
                r.P[mode] = std::stod(c[k++]);
            }

            r.betaDeg = wrap360(r.betaDeg);
            if (r.betaDeg > 180.0)
                r.betaDeg = 360.0 - r.betaDeg;

            rows_.push_back(r);
        }

        if (rows_.empty())
            throw std::runtime_error("FittedWaveForceWindowSolver: empty force table.");

        std::sort(rows_.begin(), rows_.end(),
            [](const ForceRow& a, const ForceRow& b)
        {
            if (std::abs(a.Fn - b.Fn) > 1.0e-12)
                return a.Fn < b.Fn;

            if (std::abs(a.omegaIncident - b.omegaIncident) > 1.0e-12)
                return a.omegaIncident < b.omegaIncident;

            return a.betaDeg < b.betaDeg;
        });

        std::vector<double> fns;
        std::vector<double> betas;
        std::vector<double> oms;

        for (const auto& r : rows_)
        {
            fns.push_back(r.Fn);
            betas.push_back(r.betaDeg);
            oms.push_back(r.omegaIncident);
        }

        auto sortUnique = [](std::vector<double>& x)
        {
            std::sort(x.begin(), x.end());
            x.erase(std::unique(x.begin(), x.end(),
                [](double a, double b)
            {
                return std::abs(a - b) < 1.0e-10;
            }),
                x.end());
        };

        sortUnique(fns);
        sortUnique(betas);
        sortUnique(oms);

        std::cout << "[FirstOrderForceTable loaded]\n";
        std::cout << "  rows = " << rows_.size() << "\n";

        std::cout << "  Fn axis:";
        for (double x : fns) std::cout << " " << x;
        std::cout << "\n";

        std::cout << "  beta axis:";
        for (double x : betas) std::cout << " " << x;
        std::cout << "\n";

        std::cout << "  omega axis:";
        for (double x : oms) std::cout << " " << x;
        std::cout << "\n";
    }

    ForceRow nearestFnOmegaRow(double Fn, double omegaIncident) const
    {
        const ForceRow* best = &rows_.front();
        double bestMetric = 1.0e100;

        for (const auto& r : rows_)
        {
            const double m =
                std::abs(r.Fn - Fn)
                + 0.2 * std::abs(r.omegaIncident - omegaIncident);

            if (m < bestMetric)
            {
                bestMetric = m;
                best = &r;
            }
        }

        return *best;
    }

    std::vector<ForceRow> collectSameFnOmega(
        double Fn,
        double omegaIncident) const
    {
        const ForceRow base = nearestFnOmegaRow(Fn, omegaIncident);

        std::vector<ForceRow> cand;
        for (const auto& r : rows_)
        {
            if (std::abs(r.Fn - base.Fn) < 1.0e-10 &&
                std::abs(r.omegaIncident - base.omegaIncident) < 1.0e-10)
            {
                cand.push_back(r);
            }
        }

        if (cand.empty())
            cand.push_back(base);

        std::sort(cand.begin(), cand.end(),
            [](const ForceRow& a, const ForceRow& b)
        {
            return a.betaDeg < b.betaDeg;
        });

        return cand;
    }

    std::complex<double> complexAmpPhase(
        const ForceRow& r,
        int mode) const
    {
        return std::complex<double>(
            r.A[mode] * std::cos(r.P[mode]),
            r.A[mode] * std::sin(r.P[mode]));
    }

    std::complex<double> interpComplexByBeta(
        const std::vector<ForceRow>& cand,
        double betaDeg,
        int mode) const
    {
        if (cand.empty())
            return { 0.0, 0.0 };

        if (cand.size() == 1)
            return complexAmpPhase(cand.front(), mode);

        if (betaDeg <= cand.front().betaDeg)
            return complexAmpPhase(cand.front(), mode);

        if (betaDeg >= cand.back().betaDeg)
            return complexAmpPhase(cand.back(), mode);

        for (std::size_t i = 0; i + 1 < cand.size(); ++i)
        {
            const ForceRow& r0 = cand[i];
            const ForceRow& r1 = cand[i + 1];

            if (betaDeg >= r0.betaDeg && betaDeg <= r1.betaDeg)
            {
                const double den = r1.betaDeg - r0.betaDeg;
                const double tb =
                    (std::abs(den) > 1.0e-12)
                    ? (betaDeg - r0.betaDeg) / den
                    : 0.0;

                const std::complex<double> z0 = complexAmpPhase(r0, mode);
                const std::complex<double> z1 = complexAmpPhase(r1, mode);

                return (1.0 - tb) * z0 + tb * z1;
            }
        }

        return complexAmpPhase(cand.back(), mode);
    }

    ForceEval fittedForceAt(
        const CoupledSlowState3DOF& s,
        double t,
        const CoupledWaveEnvironment& env) const
    {
        ForceEval ret;
        ret.Fmode = Eigen::VectorXd::Zero(skCfg_.DOF);
        ret.F6 = Eigen::VectorXd::Zero(6);

        const auto enc = env.evaluateEncounter(s, t);

        const double FnNow =
            currentFnFromSpeed(s.U);

        double betaRelDeg = wrap360(enc.betaRel * 180.0 / PI);
        double betaQueryDeg = betaRelDeg;

        bool mirrored = false;
        if (betaQueryDeg > 180.0)
        {
            betaQueryDeg = 360.0 - betaQueryDeg;
            mirrored = true;
        }

        // ���Խ׶ν���������Ƶ����λ������ phaseAtCG ���Ÿ���
        const double theta = enc.kin.we * t;

        double ampRef = 1.0;
        {
            double bestMetric = 1.0e100;
            const double omegaUse = nearestOmegaInTable(enc.omega);

            for (const auto& r : rows_)
            {
                const double m =
                    std::abs(r.omegaIncident - omegaUse)
                    + 10.0 * std::abs(r.Fn - FnNow)
                    + 0.01 * std::abs(r.betaDeg - betaQueryDeg);

                if (m < bestMetric)
                {
                    bestMetric = m;
                    ampRef = r.waveAmp;
                }
            }
        }

        const double ampScale =
            (std::abs(ampRef) > 1.0e-12)
            ? enc.waveAmp / ampRef
            : 1.0;

        for (int mode = 0; mode < 6; ++mode)
        {
            std::complex<double> z =
                interpForceComplex(
                    FnNow,
                    enc.omega,
                    betaQueryDeg,
                    mode);

            if (mirrored)
                z *= static_cast<double>(mirrorSignOfMode(mode));

            const double amp = std::abs(z);
            const double phase = std::atan2(z.imag(), z.real());

            ret.F6(mode) =
                ampScale * amp * std::cos(theta + phase);
        }

        for (int i = 0; i < skCfg_.DOF; ++i)
        {
            const int mode = skCfg_.modes[i];
            if (mode >= 0 && mode < 6)
                ret.Fmode(i) = ret.F6(mode);
        }

        ret.betaRelDeg = betaRelDeg;
        ret.betaQueryDeg = betaQueryDeg;
        ret.mirrored = mirrored ? 1.0 : 0.0;
        ret.theta = theta;
        ret.omegaEncounter = enc.kin.we;
        ret.waveAmp = enc.waveAmp;
        ret.etaCG = enc.waveAmp * std::cos(enc.phaseAtCG);

        return ret;
    }

    void loadHydrostaticsByOldMethod();

    static void trimVelocityHistory(
        std::vector<Eigen::VectorXd>& hist,
        int maxSize,
        int D);
};

inline void FittedWaveForceWindowSolver::loadHydrostaticsByOldMethod()
{
    const std::string elementFile = casePath_ + ship_.Name + ".element";

    if (!std::filesystem::is_regular_file(elementFile))
    {
        throw std::runtime_error(
            "FittedWaveForceWindowSolver: element file not found: " + elementFile);
    }

    auto elementData = CaseLoader::loadelement(
        elementFile,
        skCfg_.Panel.NEType,
        skCfg_.Panel.NE);

    auto element = std::make_shared<Element>(
        skCfg_.Solver,
        skCfg_.Panel.NE,
        std::move(elementData));

    Eigen::Vector3d cg;
    cg << ship_.Mass.CG.at(0),
        ship_.Mass.CG.at(1),
        ship_.Mass.CG.at(2);

    // ���ϴ��Ͳ� processElement() ˼·һ��
    element->Geometry(cg);
    element->RankineSource2();

    // ���ϴ��Ͳ�һ�£���ˮ���ֳ����� hs
    hs_ = SeakeepingDOF::hydrostaticsFromWaterline(
        ship_,
        skCfg_,
        *element);
}


inline void FittedWaveForceWindowSolver::trimVelocityHistory(
    std::vector<Eigen::VectorXd>& hist,
    int maxSize,
    int D)
{
    if (maxSize < 1)
        maxSize = 1;

    if (static_cast<int>(hist.size()) > maxSize)
    {
        hist.erase(
            hist.begin(),
            hist.begin() + (static_cast<int>(hist.size()) - maxSize));
    }

    for (auto& v : hist)
    {
        if (v.size() != D)
            v = Eigen::VectorXd::Zero(D);
    }


}