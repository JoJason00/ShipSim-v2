#include "TwoTimeScaleCoordinator.h"

#include "../maneuver/CoupledLowFreqManeuverFactory.h"
#include "../maneuver/CoupledTurningScenario.h"
#include "../maneuver/CoupledZigzagScenario.h"
#include "../../wave/CrossWave.h"
#include "../../wave/Wave.h"

#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../const/Const.h"

namespace
{
    std::string lowerCopy(std::string s)
    {
        for (char& c : s)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    }

    void validateCoupledCase(const CaseConfig& cfg)
    {
        if (!cfg.enable_coupling)
            throw std::runtime_error("Coupling mode is not enabled.");

        if (!cfg.Mmg.has_value())
            throw std::runtime_error("Coupling requires Mmg config.");

        if (!cfg.Coupled.has_value())
            throw std::runtime_error("Coupling requires Coupled config.");

        if (cfg.Coupled->physics.enableFastSeakeeping && !cfg.Seakeeping.has_value())
            throw std::runtime_error("Fast seakeeping is enabled but Seakeeping config is missing.");

        const std::string lm = lowerCopy(cfg.Mmg->lowFreqManeuverModel);
        if (lm == "sjtu_s175_nomoto" || lm == "sjtu" || lm == "sjtu_mmg")
        {
            if (!cfg.Mmg->sjtuMmg.enabled)
                throw std::runtime_error(
                    "Manoeuvring.lowFreqManeuverModel selects the SJTU MMG, but "
                    "Manoeuvring.CoupledSjtuMmg3DOF.enabled is false. Set enabled=true and supply coefficients.");
        }
    }

    CoupledSlowState3DOF makeInitialState(
        const CaseConfig& cfg,
        const ICoupledLowFreqManeuverSolver& slowSolver,
        double Fn,
        int iWave)
    {
        CoupledSlowState3DOF s{};

        s.t = 0.0;
        s.u = slowSolver.referenceSpeedFromFn(Fn);
        s.v = 0.0;
        s.r = 0.0;

        s.xe = 0.0;
        s.ye = 0.0;
        s.psi = 0.0;

        // A crossing-sea start position seeds the initial Earth-frame position.
        // The whole coupled FK pipeline already traces phase from (xe, ye), so
        // no other coupled change is needed. Single waves carry no start
        // position -> xe/ye stay 0 (unchanged trajectory).
        if (cfg.Seakeeping.has_value()
            && iWave >= 0
            && iWave < static_cast<int>(cfg.Seakeeping->waves.size()))
        {
            if (auto cw = std::dynamic_pointer_cast<CrossWave>(
                    cfg.Seakeeping->waves[static_cast<std::size_t>(iWave)]))
            {
                s.xe = cw->startX();
                s.ye = cw->startY();
            }
        }

        s.delta = 0.0;
        s.U = std::abs(s.u);
        s.betaDrift = 0.0;
        s.Fn = Fn;

        return s;
    }

    std::string coupledOutputFolderTag(const CaseConfig& cfg, const std::string& scenarioTag)
    {
        const bool zig = scenarioTag.find("zigzag") != std::string::npos;
        const std::string base = zig ? "coupled_zigzag" : "coupled_turning";
        if (!cfg.Mmg.has_value())
            return base;
        return base + lowFreqManeuverOutputFolderSuffix(*cfg.Mmg);
    }
}

TwoTimeScaleCoordinator::TwoTimeScaleCoordinator(
    const CaseConfig& cfg,
    const std::string& casePath,
    IWindowSeakeepingSolver& fastSolver,
    ICoupledLowFreqManeuverSolver& slowSolver,
    ISecondOrderLoadProvider& secondOrderProvider,
    CoupledWaveEnvironment& waveEnv)
    : cfg_(cfg),
    casePath_(casePath),
    fastSolver_(fastSolver),
    slowSolver_(slowSolver),
    secondOrderProvider_(secondOrderProvider),
    waveEnv_(waveEnv),
    writer_(casePath, cfg.Ship.Geometry.Length)
{
}

CoupledWindowRequest TwoTimeScaleCoordinator::buildWindowRequest(
    const CoupledSlowState3DOF& s0,
    const CoupledSlowState3DOF& s1pred) const
{
    CoupledWindowRequest req;
    req.slow0 = s0;
    req.slow1_pred = s1pred;
    req.t0 = s0.t;
    req.dtSlow = cfg_.Coupled->time.dtSlow;
    req.dtFast = cfg_.Coupled->time.dtFast;
    req.nFast = static_cast<int>(std::ceil(req.dtSlow / req.dtFast));
    return req;
}

void TwoTimeScaleCoordinator::runScenario(
    ICoupledManeuverScenario& scenario,
    double Fn,
    int iWave)
{
    const double dtSlow = cfg_.Coupled->time.dtSlow;

    waveEnv_.setActiveWaveIndex(iWave);

    CoupledSlowState3DOF s = makeInitialState(cfg_, slowSolver_, Fn, iWave);
    scenario.reset(s);

    const std::string scenarioFolder = coupledOutputFolderTag(cfg_, scenario.tag());

    const double U0 = slowSolver_.referenceSpeedFromFn(Fn);
    const auto waveSpec = waveEnv_.activeWaveSpec();

    // Per-wave-case folder name.
    //   Regular wave -> Fn<Fn>_nada<λ/L>_fai<heading_deg>  (LEGACY, do not change:
    //     the turning sweep + postprocess_turning_indices.py parse this exact form)
    //       Fn  : 2 decimals (e.g. 0.15)
    //       nada: λ/L_pp, 1 decimal (e.g. 1.0)
    //       fai : absolute wave direction in degrees, integer (180 = head sea)
    //   Cross / irregular -> Fn<Fn>_<Wave::conditionTag>  — a single (λ/L, heading)
    //     can't describe a crossing or spectral sea, so name both crossing systems
    //     / the spectrum instead (otherwise such runs collapse to one folder).
    // Files inside use just the scenario tag (port / starboard / zigzag).
    auto formatFixed = [](double x, int decimals) -> std::string
    {
        std::ostringstream os;
        os << std::fixed << std::setprecision(decimals) << x;
        return os.str();
    };

    const std::shared_ptr<WaveBase> wave =
        (cfg_.Seakeeping.has_value()
            && iWave >= 0
            && iWave < static_cast<int>(cfg_.Seakeeping->waves.size()))
        ? cfg_.Seakeeping->waves[static_cast<std::size_t>(iWave)]
        : nullptr;

    std::string caseDir;
    if (wave && !std::dynamic_pointer_cast<RegularWave>(wave))
    {
        caseDir = "Fn" + formatFixed(Fn, 2) + "_" + Wave::conditionTag(wave);
    }
    else
    {
        const double Lpp = cfg_.Ship.Geometry.Length;
        const double lambdaOverL =
            (waveSpec.k > 1.0e-12 && Lpp > 1.0e-12)
            ? (2.0 * PI / waveSpec.k) / Lpp
            : 0.0;
        const double headingDeg = waveSpec.directionRad * 180.0 / PI;

        std::ostringstream caseDirOs;
        caseDirOs << "Fn" << formatFixed(Fn, 2)
                  << "_nada" << formatFixed(lambdaOverL, 1)
                  << "_fai" << static_cast<int>(std::lround(headingDeg));
        caseDir = caseDirOs.str();
    }

    const std::string folderTag =
        (std::filesystem::path(scenarioFolder) / caseDir).string();
    const std::string fileTag = scenario.tag();
    if (cfg_.Seakeeping.has_value())
    {
        writer_.setScenarioReference(
            U0,
            waveSpec.omegaIncident,
            waveSpec.zetaA,
            cfg_.Seakeeping->modes);
    }

    writer_.openCase(folderTag, fileTag);
    writer_.writeSlowState(s);

    // Optional debug dump: incident wave-elevation time history at the ship
    // reference point (CG along the trajectory) — the encounter-frame η that
    // drives the coupled wave force. Enabled per wave entry via the same
    // "output_history" switch as pure seakeeping. Slow-step resolution.
    std::ofstream waveHist;
    bool dumpWaveHist = false;
    if (cfg_.Seakeeping.has_value()
        && iWave >= 0
        && iWave < static_cast<int>(cfg_.Seakeeping->waves.size())
        && cfg_.Seakeeping->waves[static_cast<std::size_t>(iWave)]->outputHistory)
    {
        const std::filesystem::path whDir =
            std::filesystem::path(casePath_) / folderTag;
        std::filesystem::create_directories(whDir);
        waveHist.open((whDir / (fileTag + "_wave_history.csv")).string());
        if (waveHist.is_open())
        {
            dumpWaveHist = true;
            waveHist << "t,eta_cg,xe,ye,psi\n";
            waveHist << s.t << "," << waveEnv_.evaluateWaveElevationAtCG(s, s.t)
                     << "," << s.xe << "," << s.ye << "," << s.psi << "\n";
        }
    }

    CoupledFastState fastCarry;

    while (true)
    {
        const CoupledManeuverCommand cmd = scenario.evaluate(s, s.t);
        if (cmd.finished)
            break;

        // The slow MMG solver never advances s.Fn (it carries the approach value
        // forever), so refresh it here from the true speed BEFORE building the
        // window request. Both req.slow0 and req.slow1_pred then carry the
        // current Fn; otherwise interpSlow ramps Fn back toward the approach
        // value inside every window and the excitation provider's mid-window
        // kernel refresh keeps re-selecting wrong-Fn kernels during the turn.
        const double U = std::sqrt(s.u * s.u + s.v * s.v);
        fastSolver_.Fn_now = U / std::sqrt(G * cfg_.Ship.Geometry.Length);
        s.Fn = fastSolver_.Fn_now;

        CoupledSlowState3DOF sPred = s;
        sPred.t += dtSlow;
        sPred.delta = cmd.deltaCmd;

        CoupledWindowRequest req = buildWindowRequest(s, sPred);
        req.fast0 = fastCarry;


        CoupledWindowResult win = fastSolver_.solveWindow(req, waveEnv_);
        fastCarry = win.fastLast;

        CoupledExternalLoads3DOF ext{};

        if (cfg_.Coupled->physics.enableSecondOrderLoads)
        {
            const auto enc = waveEnv_.evaluateEncounter(s, s.t);
            const auto ext2 = secondOrderProvider_.sample(s, enc);

            ext.X += ext2.X;
            ext.Y += ext2.Y;
            ext.N += ext2.N;

            ext.X2 = ext2.X2;
            ext.Y2 = ext2.Y2;
            ext.N2 = ext2.N2;
            // 二阶平均横倾矩 K2：之前漏拷，导致 4-DOF 横摇方程里 ext.K2 恒为 0。
            // 接通后，斜浪/横浪的平均横倾矩才会进入 integrateRoll。迎/随浪 K2≈0，无影响。
            // 注：若斜浪首测 heel 方向反了，是 K2 与 MMG 横摇符号约定不一致——
            //     在 CoupledMmg3DOFCore::integrateRoll 把 “+Kext” 改成 “-Kext” 即可。
            ext.K2 = ext2.K2;
            ext.hasSecondOrder = ext2.hasSecondOrder;
        }

        // Advance the slow MMG step (first-order wave loads are not fed back
        // to the manoeuvring solver in this version — only the second-order
        // drift loads from `ext` are).
        s = slowSolver_.step(s, ext, cmd.deltaCmd, dtSlow);

        writer_.writeSlowState(s);
        writer_.writeExternalLoads(s.t, ext);
        writer_.writeFastSummary(win, s.phi);   // s.phi = slow heel for total-roll column

        if (dumpWaveHist)
            waveHist << s.t << "," << waveEnv_.evaluateWaveElevationAtCG(s, s.t)
                     << "," << s.xe << "," << s.ye << "," << s.psi << "\n";
    }
}

void TwoTimeScaleCoordinator::run()
{
    validateCoupledCase(cfg_);

    std::cout << "\n[Coupling] ---------- seakeeping–manoeuvre run ----------\n";
    if (cfg_.Mmg.has_value())
        std::cout << "[Coupling] lowFreqManeuverModel: " << lowFreqManeuverModelTag(*cfg_.Mmg)
                  << " (JSON: \"" << cfg_.Mmg->lowFreqManeuverModel << "\")\n";
    std::cout << "[Coupling] dtSlow=" << cfg_.Coupled->time.dtSlow << " s, dtFast=" << cfg_.Coupled->time.dtFast
              << " s\n";
    std::cout << "[Coupling] waveForceMode=" << cfg_.Coupled->physics.waveForceMode
              << ", enableSecondOrderLoads=" << (cfg_.Coupled->physics.enableSecondOrderLoads ? "true" : "false")
              << ", enableFastSeakeeping=" << (cfg_.Coupled->physics.enableFastSeakeeping ? "true" : "false")
              << "\n";
    if (cfg_.Mmg.has_value())
    {
        const std::string suf = lowFreqManeuverOutputFolderSuffix(*cfg_.Mmg);
        std::cout << "[Coupling] coupled result folders: coupled_turning" << suf << ", coupled_zigzag" << suf
                  << " (suffix from lowFreqManeuverModel; default MMG uses no suffix)\n";
    }
    std::cout << "[Coupling] ------------------------------------------------\n\n";

    std::vector<double> fns;
    if (cfg_.Mmg.has_value() && !cfg_.Mmg->Fn.empty())
        fns = cfg_.Mmg->Fn;
    else if (cfg_.Seakeeping.has_value() && !cfg_.Seakeeping->Fn.empty())
        fns = cfg_.Seakeeping->Fn;
    else
        fns = { 0.15 };

    const int nWave =
        (cfg_.Seakeeping.has_value() && !cfg_.Seakeeping->waves.empty())
        ? static_cast<int>(cfg_.Seakeeping->waves.size())
        : 1;

    // Progress accounting: how many individual coupled runs in total, so each
    // start/done line can show "[k/total]" and a remaining count.
    const int scnPerWave = (cfg_.Mmg->turn_on_turningcase ? 2 : 0)
                         + (cfg_.Mmg->turn_on_zigzagcase ? 1 : 0);
    const int totalRuns = static_cast<int>(fns.size()) * nWave * scnPerWave;
    int runIdx = 0;
    std::cout << "[Coupling] batch: " << fns.size() << " Fn x " << nWave
              << " wave x " << scnPerWave << " scenario = " << totalRuns
              << " coupled run(s) total.\n\n";

    for (const double Fn : fns)
    {
        for (int iWave = 0; iWave < nWave; ++iWave)
        {
            std::string waveTag = "still-water";
            if (cfg_.Seakeeping.has_value()
                && iWave < static_cast<int>(cfg_.Seakeeping->waves.size()))
                waveTag = Wave::conditionTag(cfg_.Seakeeping->waves[iWave]);

            std::vector<std::unique_ptr<ICoupledManeuverScenario>> scenarios;

            if (cfg_.Mmg->turn_on_turningcase)
            {
                // stop_by == "t_star" 时把无量纲上限换算成有量纲时间 tMax = t*·Lpp/U0；
                // 其余情况 tMax=0，回退到按圈数停。
                double turnTMax = 0.0;
                if (cfg_.Mmg->Turning.stopBy == "t_star" && cfg_.Mmg->Turning.t_star_max > 0.0)
                {
                    const double Lpp = cfg_.Ship.Geometry.Length;
                    const double U0 = slowSolver_.referenceSpeedFromFn(Fn);
                    if (U0 > 1.0e-9)
                        turnTMax = cfg_.Mmg->Turning.t_star_max * Lpp / U0;
                }
                scenarios.push_back(std::make_unique<CoupledTurningScenario>(cfg_.Mmg->Turning, true,  0.0, turnTMax));
                scenarios.push_back(std::make_unique<CoupledTurningScenario>(cfg_.Mmg->Turning, false, 0.0, turnTMax));
            }

            if (cfg_.Mmg->turn_on_zigzagcase)
            {
                scenarios.push_back(std::make_unique<CoupledZigzagScenario>(cfg_.Mmg->Zigzag));
            }

            for (auto& sc : scenarios)
            {
                ++runIdx;
                std::cout << "\n[Coupling][" << runIdx << "/" << totalRuns
                          << "] START  Fn=" << Fn
                          << "  iWave=" << iWave << "/" << (nWave - 1)
                          << "  wave=" << waveTag
                          << "  scenario=" << sc->tag() << std::endl;

                runScenario(*sc, Fn, iWave);

                std::cout << "[Coupling][" << runIdx << "/" << totalRuns
                          << "] DONE   Fn=" << Fn
                          << "  wave=" << waveTag
                          << "  scenario=" << sc->tag()
                          << "  (" << (totalRuns - runIdx) << " remaining)"
                          << std::endl;
            }
        }
    }
}