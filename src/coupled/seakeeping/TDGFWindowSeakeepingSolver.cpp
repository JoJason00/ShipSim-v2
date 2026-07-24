#include "TDGFWindowSeakeepingSolver.h"

#include "../common/CoupledMath.h"
#include "../../const/Const.h"
#include "../../seakeeping/RollDamping.h"
#include "../../seakeeping/SeakeepingDOF.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <seakeeping/RadiationKernelTailSplit.h>
#include <filesystem>

namespace
{
    // Allocate every history field of a CoupledWindowResult to its final shape.
    void allocateWindowOutput(CoupledWindowResult& out, int N, int D)
    {
        out.loads = CoupledWaveLoadsWindow{};
        out.tHist              = Eigen::VectorXd::Zero(N);
        out.qHist              = Eigen::MatrixXd::Zero(N, D);
        out.vHist              = Eigen::MatrixXd::Zero(N, D);
        out.forceHist          = Eigen::MatrixXd::Zero(N, D);
        out.force6Hist         = Eigen::MatrixXd::Zero(N, 6);
        out.betaRelDegHist     = Eigen::VectorXd::Zero(N);
        out.betaQueryDegHist   = Eigen::VectorXd::Zero(N);
        out.mirroredHist       = Eigen::VectorXd::Zero(N);
        out.thetaHist          = Eigen::VectorXd::Zero(N);
        out.omegaEncounterHist = Eigen::VectorXd::Zero(N);
        out.waveAmpHist        = Eigen::VectorXd::Zero(N);
        out.etaCGHist          = Eigen::VectorXd::Zero(N);
    }

    // Build the stepper's initial state from a CoupledFastState carry, or zero.
    cummins::StepperState makeStepperInit(const CoupledFastState& fast, int D)
    {
        cummins::StepperState s;
        if (fast.initialized)
        {
            s.q     = fast.q;
            s.v     = fast.v;
            s.a     = fast.a;
            s.vHist = fast.vMemHist;
        }
        else
        {
            s.q = Eigen::VectorXd::Zero(D);
            s.v = Eigen::VectorXd::Zero(D);
            s.a = Eigen::VectorXd::Zero(D);
        }
        return s;
    }
}

TDGFWindowSeakeepingSolver::TDGFWindowSeakeepingSolver(
    const ShipConfig& ship,
    const SeakeepingConfig& skCfg,
    const std::string& casePath,
    const CoupledRadiationKernelRepo& radRepo,
    const HydroRefreshConfig& refreshCfg,
    IWaveForceProvider& forceProvider,
    std::shared_ptr<Element> elementForHydrostatics,
    bool enableRadiation)
    : ship_(ship),
      skCfg_(skCfg),
      casePath_(casePath),
      radRepo_(radRepo),
      refreshCfg_(refreshCfg),
      forceProvider_(forceProvider),
      element_(std::move(elementForHydrostatics)),
      enableRadiation_(enableRadiation)
{
    if (!element_)
        throw std::runtime_error("TDGFWindowSeakeepingSolver: element is null (needed for hydrostatics).");

    hs_ = SeakeepingDOF::hydrostaticsFromWaterline(ship_, skCfg_, *element_);
    SeakeepingDOF::buildSystemMatrices(ship_, skCfg_, hs_, Mphys_, Bphys_, Cphys_);

    rollIndex_ = SeakeepingDOF::findModeIndex(skCfg_, MODE_ROLL);
    rollVisc_  = RollDampingBuilder::build(
        casePath_, skCfg_.RollDamping, ship_.Mass.Mass, ship_.Mass.GM);

    std::cout << "[TDGFWindowSeakeepingSolver] init\n"
              << "  hs c33=" << hs_.c33
              << " c44=" << hs_.c44
              << " c55=" << hs_.c55
              << " c35=" << hs_.c35
              << "\n  rollIndex=" << rollIndex_
              << "\n  speedTol=" << refreshCfg_.speedTolRatio
              << " betaTolDeg=" << refreshCfg_.betaTolDeg
              << " omegaTol=" << refreshCfg_.omegaTolRatio
              << "\n  enableRadiation=" << (enableRadiation_ ? "true" : "false")
              << "\n  rollVisc (applied): B44_lin=" << rollVisc_.B44_lin
              << " B44_quad=" << rollVisc_.B44_quad
              << " B44_cube=" << rollVisc_.B44_cube
              << "\n";
}

bool TDGFWindowSeakeepingSolver::needRadRefresh(double Fn, double dtFast) const
{
    if (!radValid_ || !stepper_)                    return true;
    if (std::abs(dtFast - dtCached_) > 1e-12)        return true;

    const double FnRef = std::max(std::abs(FnRadCached_), 1e-3);
    double fnTolAbs = refreshCfg_.speedTolRatio * FnRef;

    // Adaptive policy: derive the Fn refresh tolerance from the radiation
    // library's actual Fn grid spacing, so a denser library refreshes more
    // often. Falls back to the fixed ratio when <2 distinct Fn nodes.
    if (refreshCfg_.adaptiveRefresh)
    {
        const double f = (refreshCfg_.adaptiveRefreshFactor > 1e-6)
                              ? refreshCfg_.adaptiveRefreshFactor
                              : 0.5;
        bool valid = false;
        const double dFn = radRepo_.fnSpacingMedian(valid);
        if (valid)
            fnTolAbs = std::max(1e-4, f * dFn);
    }

    if (std::abs(Fn - FnRadCached_) > fnTolAbs)
        return true;

    return false;
}

void TDGFWindowSeakeepingSolver::postprocessRadiationKernel(
    double Fn,
    RadiationKernelData& kernel) const
{
    KernelPost::KernelScaleInfo scale;
    scale.L = ship_.Geometry.Length;
    scale.displacement = ship_.Geometry.Displacement;
    scale.rho = rho;
    scale.g = G;

    KernelPost::TailSplitOptions splitOpt;

    splitOpt.enabled = true;

    // 方案 A：
    // Klag -= Kinf
    // C_prime += Kinf
    splitOpt.addTailToCPrime = true;

    // 你的编号：2 = heave，4 = pitch
    splitOpt.targetModes = { 2, 4 };
    splitOpt.onlyTargetBlock = true;

    splitOpt.tailFraction = 0.20;
    splitOpt.minTailSamples = 30;

    splitOpt.minMeanToPeak = 0.05;
    splitOpt.maxOscToMean = 0.70;

    // K(2,4)/K(4,2) 实际上趋于 0，不需要强制 tail-split；只 K(4,4) 需要。
    // 让 tail-split 自动根据 meanToPeak 决定。

    splitOpt.verbose = true;

    KernelPost::TailSplitWorkflowOptions workflowOpt;

    workflowOpt.enabled = true;

    // 调试阶段建议 true；正式批量计算可改 false。
    //workflowOpt.writeBeforeHistory = true;
    //workflowOpt.writeAfterHistory = true;

    //workflowOpt.writeDimensional = true;
    //workflowOpt.writeNondimensional = true;

    //workflowOpt.outputDir =
    //    (std::filesystem::path(filePath) / "kernel_cache").string();

    //workflowOpt.tag = "Fn" + keyDouble(Fn);

    KernelPost::applyTailSplitWorkflowInPlace(
        kernel,
        scale,
        splitOpt,
        workflowOpt);

    //if (kRadiationKernelSgSmooth && !kernel.Klag.empty())
    //{
    //    impulse_sg::Options sgOpt;
    //    std::vector<Eigen::MatrixXd> kBefore;
    //    if (kRadiationKernelSgCompareCsv)
    //        kBefore = kernel.Klag;
    //    impulse_sg::smoothMatrixLagHistoryInPlace(kernel.Klag, sgOpt);
    //    if (!kernel.Klag.empty())
    //        kernel.K0 = kernel.Klag[0];
    //    if (kRadiationKernelSgCompareCsv && !kBefore.empty())
    //    {
    //        const std::string csvDir =
    //            (std::filesystem::path(workflowOpt.outputDir) / "kernel_sg_compare").string();
    //        const std::string csvPath =
    //            (std::filesystem::path(csvDir) / ("Klag_compare_" + workflowOpt.tag + ".csv")).string();
    //        const std::vector<std::pair<int, int>> pairs = {
    //            {2, 2}, {2, 4}, {4, 4}, {4, 2} }; // 0-based: heave=2, pitch=4
    //        if (impulse_sg::writeMatrixLagCompareCsv(csvPath, kernel.dt, pairs, kBefore, kernel.Klag))
    //            std::cout << "[ImpulseKernelSG] wrote " << csvPath << "\n";
    //    }
    //}

}

void TDGFWindowSeakeepingSolver::refreshRadiation(double Fn, double dtFast)
{
    const int D = skCfg_.DOF;

    cummins::StepperConfig cfg;

    if (enableRadiation_)
    {
        radCached_   = radRepo_.getKernel(Fn, dtFast);
        FnRadCached_ = Fn;
        dtCached_    = dtFast;
        radValid_    = true;

		postprocessRadiationKernel(Fn, radCached_);

        cfg.Mtotal  = Mphys_ + radCached_.A_inf;
        cfg.Bconst  = Bphys_ + radCached_.B + 0.5 * dtFast * radCached_.K0;
        cfg.Ctotal  = Cphys_ + radCached_.C_prime;
        cfg.Konline = radCached_.Klag;

        std::cout << "[TDGF] refresh radiation kernel at Fn=" << Fn
                  << "  coupled_dtFast=" << dtFast
                  << "  TG_resampled=" << radCached_.TG
                  << "  (impulse memory uses this grid; library dt/TG printed at scan)\n";
    }
    else
    {
        // Radiation memory disabled — bring-up path for direct-pressure FK only.
        // TODO: re-enable once per-Fn kernel database is complete.
        FnRadCached_ = Fn;
        dtCached_    = dtFast;
        radValid_    = true;

        cfg.Mtotal  = Mphys_;
        cfg.Bconst  = Bphys_;
        cfg.Ctotal  = Cphys_;
        cfg.Konline = { Eigen::MatrixXd::Zero(D, D) };  // stepper requires at least K[0]

        std::cout << "[TDGF] radiation disabled — Newmark on M_phys+C_phys+rollVisc only"
                  << " (Fn=" << Fn << " dt=" << dtFast << ")\n";
    }

    cfg.dt        = dtFast;
    cfg.rollIndex = rollIndex_;
    cfg.rollVisc  = rollVisc_;

    stepper_ = std::make_unique<cummins::CumminsTimeStepper>();
    stepper_->init(std::move(cfg));
}

CoupledSlowState3DOF TDGFWindowSeakeepingSolver::interpSlow(
    const CoupledSlowState3DOF& a,
    const CoupledSlowState3DOF& b,
    double alpha)
{
    CoupledSlowState3DOF s{};
    const double om = 1.0 - alpha;

    s.t          = om * a.t          + alpha * b.t;
    s.u          = om * a.u          + alpha * b.u;
    s.v          = om * a.v          + alpha * b.v;
    s.r          = om * a.r          + alpha * b.r;
    s.xe         = om * a.xe         + alpha * b.xe;
    s.ye         = om * a.ye         + alpha * b.ye;
    s.psi        = om * a.psi        + alpha * b.psi;
    s.delta      = om * a.delta      + alpha * b.delta;
    s.U          = om * a.U          + alpha * b.U;
    s.betaDrift  = om * a.betaDrift  + alpha * b.betaDrift;
    s.Fn         = om * a.Fn         + alpha * b.Fn;
    s.rDotRadPerS2 = om * a.rDotRadPerS2 + alpha * b.rDotRadPerS2;
    return s;
}

CoupledWindowResult TDGFWindowSeakeepingSolver::solveWindow(
    const CoupledWindowRequest& req,
    const CoupledWaveEnvironment& env)
{
    const int D = skCfg_.DOF;
    const int N = std::max(1, req.nFast);
    const double dt = req.dtFast;

    CoupledWindowResult out;
    allocateWindowOutput(out, N, D);

    // Radiation kernel + stepper: rebuild only when |ΔFn| crosses the threshold.
    //if (needRadRefresh(req.slow0.Fn, dt))
    //{
    //    refreshRadiation(req.slow0.Fn, dt);
    //}

    if (needRadRefresh(Fn_now, dt))
    {
        refreshRadiation(Fn_now, dt);
    }
    stepper_->setInitialState(makeStepperInit(req.fast0, D));

    // Excitation provider continuation (impulse-kernel mode uses η history;
    // direct-FK provider ignores both calls).
    forceProvider_.setEtaHistory(req.fast0.etaHist);
    const auto enc0 = env.evaluateEncounter(req.slow0, req.t0);
    CoupledFastWindowInfo win{};
    win.slowEnd = &req.slow1_pred;
    win.t0      = req.t0;
    win.t1      = (N > 1) ? (req.t0 + static_cast<double>(N - 1) * dt) : req.t0;


    forceProvider_.onWindowBegin(req.slow0, enc0, dt, win);

    // Step through the window.
    for (int n = 0; n < N; ++n)
    {
        const double alpha = (N > 1) ? (static_cast<double>(n) / static_cast<double>(N - 1)) : 0.0;
        const CoupledSlowState3DOF sSlow = interpSlow(req.slow0, req.slow1_pred, alpha);

        const double t = req.t0 + n * dt;
        const auto encNow = env.evaluateEncounter(sSlow, t);

        Eigen::RowVectorXd F6;
        const Eigen::VectorXd Fexc = forceProvider_.evaluate(sSlow, t, encNow, &F6);

        stepper_->step(Fexc);

        out.tHist(n)              = t;
        out.qHist.row(n)          = stepper_->state().q.transpose();
        out.vHist.row(n)          = stepper_->state().v.transpose();
        out.forceHist.row(n)      = Fexc.transpose();
        out.force6Hist.row(n)     = F6;

        {
            double lookupDir = 0.0;
            bool mirrorPS = false;
            double betaRelDeg = 0.0;
            double betaQueryDeg = 0.0;
            coupled::waveHeadingLookup(
                encNow.betaRel, lookupDir, mirrorPS, betaRelDeg, betaQueryDeg);
            out.betaRelDegHist(n)   = betaRelDeg;
            out.betaQueryDegHist(n) = betaQueryDeg;
            out.mirroredHist(n)     = mirrorPS ? 1.0 : 0.0;
        }

        out.thetaHist(n)          = encNow.phaseAtCG;
        out.omegaEncounterHist(n) = encNow.kin.we;
        out.waveAmpHist(n)        = encNow.waveAmp;
        out.etaCGHist(n)          = encNow.waveAmp * std::cos(encNow.phaseAtCG);
    }

    // Window statistics.
    out.etaMean = out.qHist.colwise().mean().transpose();
    out.etaRms  = Eigen::VectorXd::Zero(D);
    for (int j = 0; j < D; ++j)
    {
        const Eigen::ArrayXd x = out.qHist.col(j).array();
        out.etaRms(j) = std::sqrt((x * x).mean());
    }

    // Hand off state for the next window.
    const auto& st = stepper_->state();
    out.fastLast.q           = st.q;
    out.fastLast.v           = st.v;
    out.fastLast.a           = st.a;
    out.fastLast.initialized = true;
    out.fastLast.vMemHist    = st.vHist;
    out.fastLast.etaHist     = forceProvider_.etaHistory();

    return out;
}
