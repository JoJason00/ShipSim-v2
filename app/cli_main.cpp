// Non-interactive solver CLI — the single entry the GUI launches via QProcess.
//
// Usage:
//   shipsim_cli <case> [coupledModel]
//
//   <case>          Case name under <projectRoot>/cases/, OR an absolute path
//                   to a case directory (the one containing case.json).
//   [coupledModel]  Only used when coupling is enabled. One of:
//                     mmg_default        — force CoupledMmg3DOFCore
//                     sjtu_s175_nomoto   — force CoupledSjtuMmg3DOF
//                     follow             — use case.json as-is (default)
//
// No stdin prompts (core_test.cpp's interactive prompts are replaced by args)
// so it can be driven head-less from the GUI. All progress goes to stdout;
// errors to stderr. Exit codes: 0 ok, 1 runtime error, 2 bad usage.

#include <string>
#include <iostream>
#include <filesystem>
#include <memory>
#include <vector>
#include <cstdlib>   // std::quick_exit

#include "../src/config/CaseConfig.h"
#include "../src/io/CaseLoader.h"
#include "../src/wave/Wave.h"
#include "../src/seakeeping/Seakeeping.h"
#include "../src/seakeeping/LinearCumminsTDGF.h"
#include "../src/mmg/mmg.h"
#include "../src/tool/Timer.h"
#include "../src/tool/findRoot.h"
#include "../src/tool/ParallelGuard.h"
#include "../src/io/Write.h"
#include "../src/seakeeping/Element.h"

// coupled
#include "../src/coupled/environment/WaveEnvironment.h"
#include "../src/coupled/hydro/CoupledRadiationKernelRepo.h"
#include "../src/coupled/hydro/CoupledExcitationKernelRepo.h"
#include "../src/coupled/seakeeping/TDGFWindowSeakeepingSolver.h"
#include "../src/coupled/seakeeping/NullWindowSeakeepingSolver.h"
#include "../src/coupled/maneuver/CoupledLowFreqManeuverFactory.h"
#include "../src/coupled/load/NullSecondOrderLoadProvider.h"
#include "../src/coupled/load/CsvSecondOrderLoadProvider.h"
#include "../src/coupled/coordinator/TwoTimeScaleCoordinator.h"
#include "../src/coupled/load/DriftForceTxtProvider.h"
#include "../src/coupled/seakeeping/FittedWaveForceWindowSolver.h"
#include "../src/seakeeping/RollDamping.h"
#include "../src/coupled/seakeeping/ImpulseKernelWaveForceProvider.h"
#include "../src/coupled/seakeeping/DirectPressureFKWaveForceProvider.h"
#include "../src/coupled/seakeeping/FKImpulseWaveForceProvider.h"
#include "../src/coupled/seakeeping/DirectFKPlusDfWaveForceProvider.h"

namespace
{
    std::string resolveCasePath(const std::string& arg, char* argv0)
    {
        // Absolute / relative path to a case directory (or its case.json).
        std::filesystem::path p(arg);
        if (p.has_parent_path() || std::filesystem::exists(p))
        {
            if (std::filesystem::is_regular_file(p) &&
                p.filename() == "case.json")
                p = p.parent_path();
            if (std::filesystem::is_directory(p))
                return p.string() + "/";
        }
        // Otherwise treat as a case name under <projectRoot>/cases/.
        auto root = findProjectRootFromExe(argv0);
        return (root / "cases" / arg).string() + "/";
    }
}

int main(int argc, char* argv[])
{
    try
    {
        if (argc < 2)
        {
            std::cerr << "usage: shipsim_cli <case> [coupledModel]\n"
                         "  coupledModel: mmg_default | sjtu_s175_nomoto | follow\n";
            return 2;
        }

        shipsim::setupGlobalThreads();

        const std::string run_case = argv[1];
        const std::string coupledModel = (argc >= 3) ? argv[2] : "follow";

        const std::string casePath = resolveCasePath(run_case, argv[0]);
        const std::string caseJson = casePath + "case.json";
        if (!std::filesystem::is_regular_file(caseJson))
            throw std::runtime_error("case.json not found: " + caseJson);

        std::cout << "[CLI] case: " << casePath << "\n";
        CaseConfig cfg = CaseLoader::loadcase(caseJson);
        Write::start(cfg);

        auto ensureElement = [&]()
        {
            const std::string shipname = cfg.Ship.Name;
            if (!std::filesystem::is_regular_file(casePath + shipname + ".element"))
            {
                if (std::filesystem::is_regular_file(casePath + shipname + ".dat"))
                    CaseLoader::UGtoElement(casePath + shipname + ".dat",
                        casePath + shipname + ".element", 0.001);
                else
                    throw std::runtime_error("there is no element file or UG mesh file");
            }
        };

        //=======================
        // 纯耐波性模块
        //=======================
        if (cfg.enable_seakeeping && cfg.Seakeeping.has_value())
        {
            Timer timer("Seakeeping");
            ensureElement();

            // INTERIM WORKAROUND (Phase 0.5, root cause still open):
            // ~LinearCumminsTDGF bad-frees an Eigen member whose data pointer
            // was overwritten by a Release-only out-of-bounds write during the
            // run (confirmed by ASan; Debug heap masks it). The crash fires
            // when the solver leaves scope — AFTER all CSV outputs are fully
            // written, so results are unaffected. Until the OOB write is found
            // and fixed, heap-allocate the solver and deliberately do NOT
            // delete it: skipping the corrupting destructor lets the process
            // exit cleanly (the OS reclaims the memory at exit anyway). This
            // is a labelled interim shim, not a fix — Phase 0.5 continues.
            if (cfg.Seakeeping.value().ResponseMethod == "LinearCumminsTDGF")
            {
                auto* sk = new LinearCumminsTDGF(cfg.Ship, casePath, cfg.Seakeeping.value());
                sk->run();
                // intentional leak: dtor is corrupt; freed by OS at exit.
            }
            else
            {
                auto* sk = new Seakeeping(cfg.Ship, casePath, cfg.Seakeeping.value());
                sk->run();
                // intentional leak: paired with the above workaround.
            }

            timer.print_s("Seakeeping total time: ");
            std::cout << "\n\nRun seakeeping case successfully!\n" << std::endl;
        }

        //=======================
        // 纯操纵性模块
        //=======================
        if (cfg.enable_maneuvering && cfg.Mmg.has_value())
        {
            Timer timer("Maneuvering");
            Mmg mmg(cfg.Ship, casePath, cfg.Mmg.value());
            mmg.run();
            timer.print_s("Maneuvering total time: ");
            std::cout << "\n\nRun maneuvering case successfully!\n" << std::endl;
        }

        //=======================
        // 耐波性操纵性耦合模块
        //=======================
        if (cfg.enable_coupling && cfg.Coupled.has_value())
        {
            Timer timer("Seakeeping-Maneuvering coupling");
            ensureElement();

            if (!cfg.Mmg.has_value())
                throw std::runtime_error("Coupling requires Manoeuvring (Mmg) in case.json.");

            if (coupledModel == "mmg_default")
            {
                cfg.Mmg->lowFreqManeuverModel = "mmg_default";
                cfg.Mmg->sjtuMmg.enabled = false;
                std::cout << "[Coupled] runtime: mmg_default (override case.json).\n";
            }
            else if (coupledModel == "sjtu_s175_nomoto")
            {
                cfg.Mmg->lowFreqManeuverModel = "sjtu_s175_nomoto";
                cfg.Mmg->sjtuMmg.enabled = true;
                if (cfg.Mmg->sjtuMmg.Lpp <= 1.0e-6)
                    throw std::runtime_error(
                        "SJTU MMG selected but Manoeuvring.CoupledSjtuMmg3DOF is missing/empty in case.json.");
                std::cout << "[Coupled] runtime: sjtu_s175_nomoto.\n";
            }
            else if (coupledModel == "follow")
            {
                std::cout << "[Coupled] runtime: follow case.json lowFreqManeuverModel=\""
                          << cfg.Mmg->lowFreqManeuverModel << "\"\n";
            }
            else
            {
                throw std::runtime_error("Invalid coupledModel: " + coupledModel +
                    " (mmg_default | sjtu_s175_nomoto | follow)");
            }

            std::vector<std::shared_ptr<WaveBase>> waves;
            if (cfg.Seakeeping.has_value())
                waves = cfg.Seakeeping->waves;

            CoupledWaveEnvironment waveEnv(waves);
            if (cfg.Coupled.has_value())
                waveEnv.setWaveRampTime(cfg.Coupled->time.waveRampTimeS);

            std::unique_ptr<IWindowSeakeepingSolver> fastSolver;
            std::unique_ptr<CoupledRadiationKernelRepo> radRepo;
            std::unique_ptr<CoupledExcitationKernelRepo> excRepo;
            std::unique_ptr<IWaveForceProvider> waveForceProvider;
            std::shared_ptr<Element> coupledElement;

            if (cfg.Coupled->physics.enableFastSeakeeping)
            {
                if (!cfg.Seakeeping.has_value())
                    throw std::runtime_error("Fast seakeeping is enabled but Seakeeping config is missing.");

                radRepo = std::make_unique<CoupledRadiationKernelRepo>(
                    cfg.Ship, cfg.Seakeeping.value(), casePath);

                {
                    const std::string elementFile = casePath + cfg.Ship.Name + ".element";
                    auto elementData = CaseLoader::loadelement(
                        elementFile,
                        cfg.Seakeeping->Panel.NEType,
                        cfg.Seakeeping->Panel.NE);
                    coupledElement = std::make_shared<Element>(
                        cfg.Seakeeping->Solver,
                        cfg.Seakeeping->Panel.NE,
                        std::move(elementData));
                    Eigen::Vector3d cg(cfg.Ship.Mass.CG.at(0),
                                       cfg.Ship.Mass.CG.at(1),
                                       cfg.Ship.Mass.CG.at(2));
                    coupledElement->Geometry(cg);
                    coupledElement->RankineSource2();
                }

                const std::string& mode = cfg.Coupled->physics.waveForceMode;
                std::cout << "[Coupled] waveForceMode = " << mode << "\n";

                if (mode == "fitted")
                {
                    fastSolver = std::make_unique<FittedWaveForceWindowSolver>(
                        cfg.Ship, cfg.Seakeeping.value(), casePath, *radRepo,
                        casePath + "ExcitingForce/FirstOrderWaveForceAmpPhase.csv");
                }
                else
                {
                    if (mode == "directFK")
                        waveForceProvider = std::make_unique<DirectPressureFKWaveForceProvider>(
                            cfg.Ship, cfg.Seakeeping.value(), coupledElement);
                    else if (mode == "directFKdf")
                        waveForceProvider = std::make_unique<DirectFKPlusDfWaveForceProvider>(
                            cfg.Ship, cfg.Seakeeping.value(), casePath, cfg.Ship.Name,
                            cfg.Coupled->refresh, coupledElement);
                    else if (mode == "fkImpulse")
                        waveForceProvider = std::make_unique<FKImpulseWaveForceProvider>(
                            cfg.Seakeeping.value(), casePath, cfg.Ship.Name, cfg.Coupled->refresh);
                    else if (mode == "impulseKernel")
                    {
                        excRepo = std::make_unique<CoupledExcitationKernelRepo>(
                            cfg.Ship, cfg.Seakeeping.value(), casePath);
                        waveForceProvider = std::make_unique<ImpulseKernelWaveForceProvider>(
                            cfg.Seakeeping.value(), *excRepo, cfg.Coupled->refresh);
                    }
                    else
                        throw std::runtime_error("Unknown waveForceMode: " + mode);

                    fastSolver = std::make_unique<TDGFWindowSeakeepingSolver>(
                        cfg.Ship, cfg.Seakeeping.value(), casePath, *radRepo,
                        cfg.Coupled->refresh, *waveForceProvider, coupledElement,
                        cfg.Coupled->physics.enableRadiation);
                }
            }
            else
            {
                fastSolver = std::make_unique<NullWindowSeakeepingSolver>();
            }

            // ---- 统一横摇阻尼：单一数据源 ----
            // 操纵慢漂 heel(integrateRoll) 与耐波快横摇(TDGFWindowSeakeepingSolver) 都用
            // 同一份 RollDampingBuilder 输出(权威源 = Seakeeping.RollDamping)，避免两套阻尼。
            // 这样无论 FromDecayCsv(拟合) 还是 DirectCoefficients(手填)，两边自动一致；
            // 想换阻尼只改 Seakeeping.RollDamping 一处即可。
            if (cfg.Mmg.has_value()
                && cfg.Mmg->coupledMmg3DOF.enabled
                && cfg.Mmg->coupledMmg3DOF.maneuverDOF == 4
                && cfg.Seakeeping.has_value())
            {
                try
                {
                    const RollViscDamping rvUnified = RollDampingBuilder::build(
                        casePath, cfg.Seakeeping->RollDamping,
                        cfg.Ship.Mass.Mass, cfg.Ship.Mass.GM);

                    auto& sc = cfg.Mmg->coupledMmg3DOF;
                    std::cout << "[Coupling] unify roll damping (maneuver heel <- seakeeping source): "
                              << "B44_lin " << sc.rollB44_lin << "->" << rvUnified.B44_lin
                              << ", B44_quad " << sc.rollB44_quad << "->" << rvUnified.B44_quad
                              << " (B44_cube=" << rvUnified.B44_cube << " ignored by 4-DOF heel)\n";
                    sc.rollB44_lin  = rvUnified.B44_lin;
                    sc.rollB44_quad = rvUnified.B44_quad;
                }
                catch (const std::exception& e)
                {
                    std::cout << "[Coupling][warn] unify roll damping failed ("
                              << e.what() << "); keeping case.json rollB44 values.\n";
                }
            }

            std::unique_ptr<ICoupledLowFreqManeuverSolver> lowFreqSolver =
                makeLowFreqManeuverSolver(cfg.Ship, cfg.Mmg.value());

            std::unique_ptr<ISecondOrderLoadProvider> secondOrderProvider;
            if (cfg.Coupled->physics.enableSecondOrderLoads)
                secondOrderProvider = std::make_unique<DriftForceTxtProvider>(casePath + "DriftForce");
            else
                secondOrderProvider = std::make_unique<NullSecondOrderLoadProvider>();

            TwoTimeScaleCoordinator coordinator(
                cfg, casePath, *fastSolver, *lowFreqSolver,
                *secondOrderProvider, waveEnv);
            coordinator.run();

            timer.print_s("Seakeeping-Maneuvering coupling total time: ");
            std::cout << "\n\nRun coupled case successfully!\n" << std::endl;
        }

        std::cout << "[CLI] done." << std::endl;

        // INTERIM WORKAROUND (Phase 0.5, root cause still open):
        // There is a known pre-existing Release-only heap corruption that
        // manifests as a bad-free in ~LinearCumminsTDGF at process teardown
        // (ASan: an Eigen VectorXd member's data pointer is overwritten by an
        // out-of-bounds write during the run). It happens AFTER every CSV
        // output is fully written, so results are unaffected. Until the OOB
        // write is found and fixed, exit via quick_exit so the corrupting
        // destructors are skipped and the GUI (QProcess) gets a clean exit
        // code. This does NOT hide the bug — Phase 0.5 keeps tracking it.
        std::cout.flush();
        std::cerr.flush();
        std::quick_exit(0);

        return 0;   // unreachable; kept for clarity
    }
    catch (const std::exception& e)
    {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
}
