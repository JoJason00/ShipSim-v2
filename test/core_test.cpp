#include <string>
#include <iostream>
#include "../src/config/CaseConfig.h"
#include "../src/io/CaseLoader.h"
#include "../src/wave/Wave.h"
//#include "../jsoncpp.cpp"
#include "../src/seakeeping/Seakeeping.h"
#include "../src/mmg/mmg.h"
#include "../src/tool/Timer.h"
#include "../src/tool/findRoot.h"
#include "../src/io/Write.h"
#include <filesystem>
#include <memory>
#include <seakeeping/LinearCumminsTDGF.h>
#include "seakeeping/WigleyIWritePanels.h"

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
#include "../src/coupled/seakeeping/ImpulseKernelWaveForceProvider.h"
#include "../src/coupled/seakeeping/DirectPressureFKWaveForceProvider.h"
#include "../src/coupled/seakeeping/FKImpulseWaveForceProvider.h"
#include "../src/coupled/seakeeping/DirectFKPlusDfWaveForceProvider.h"
#include "../src/seakeeping/Element.h"


int main(int argc, char* argv[])
{
    try
    {
        //=======================
        //  加载工况
        //=======================
        std::vector<std::string> cases = { "wigleyI","wigleyI_Lu","S175_1424","S175_1262","S175_924" ,
            "S175_848", "S175_1722",  "S175_1476", "S175_1874", "S175_3418", "S175_18016" };

		std::cout << "chose a case to run:\n";
		for (size_t i = 0; i < cases.size(); ++i)
			std::cout << "  " << i << ":\t " << cases[i] << "\n";
		std::cout << "input case number:\t";
        int caseNum;
        std::cin >> caseNum;

        std::string run_case = (argc >= 2) ? argv[1] : cases[caseNum];
        auto root = findProjectRootFromExe(argv[0]);
        std::string casePath = (root / "cases" / run_case).string() + "/";
        CaseConfig cfg = CaseLoader::loadcase(casePath + "case.json");

        //NE = 2(NS−1)[(NP−2)+WL_NS * WL_NP],水线边缘网格加密
        //wigleyi::WriteWigleyISeakeepingPanels(
        //    30.0,
        //    3.0,
        //    1.875,
        //    casePath + "wigleyI.element",
        //    21,   // NS
        //    7,    // NP
        //    3,    // WL_NS：水线附近每个原始面元沿船长方向切 4 份
        //    2,    // WL_NP：水线附近每个原始面元沿垂向切 2 份
        //    true,
        //    3
        //);

        Write::start(cfg);


        //=======================
        // 纯耐波性模块
        //=======================
        if (cfg.enable_seakeeping && cfg.Seakeeping.has_value())
        {
            Timer timer("Seakeeping");

            //优先读取面元.element文件，没有则读取UG网格文件并转换成面元文件
            std::string shipname = cfg.Ship.Name;
            //std::string shipname = run_case;
            if (!std::filesystem::is_regular_file(casePath + shipname + ".element"))
            {
                if (std::filesystem::is_regular_file(casePath + shipname + ".dat"))
                    CaseLoader::UGtoElement(casePath + shipname + ".dat",
                        casePath + shipname + ".element", 0.001);
                else
                    throw std::runtime_error("there is no element file or UG mesh file");
            }

            //Seakeeping seakeepingCase(cfg.Ship, casePath, cfg.Seakeeping.value());
            //seakeepingCase.run();

            if (cfg.Seakeeping.value().ResponseMethod == "LinearCumminsTDGF")
            {
                LinearCumminsTDGF sk(cfg.Ship, casePath, cfg.Seakeeping.value());
                sk.run();
            }
            else
            {
                Seakeeping sk(cfg.Ship, casePath, cfg.Seakeeping.value());
                sk.run();
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
            /*RAO4DProvider rao;
            auto info = rao.loadIfExists(casePath + "RAO.csv");
            std::cout << info.message << "\n";
            if (!info.ok)
                throw std::runtime_error("load RAO data file error");*/
                //auto Hk = rao.evalDeg(Fn, betaDeg, omega);
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

            std::string shipname = cfg.Ship.Name;
            if (!std::filesystem::is_regular_file(casePath + shipname + ".element"))
            {
                if (std::filesystem::is_regular_file(casePath + shipname + ".dat"))
                    CaseLoader::UGtoElement(casePath + shipname + ".dat",
                        casePath + shipname + ".element", 0.001);
                else
                    throw std::runtime_error("there is no element file or UG mesh file");
            }

            if (!cfg.Mmg.has_value())
                throw std::runtime_error("Coupling requires Manoeuvring (Mmg) in case.json.");

            std::cout << "\nChoose low-frequency manoeuvre model for coupling (like case selection):\n";
            std::cout << "  0: mmg_default — CoupledMmg3DOFCore, results under coupled_turning / coupled_zigzag\n";
            std::cout << "  1: sjtu_s175_nomoto — CoupledSjtuMmg3DOFCore, results under coupled_turning_sjtu / ...\n";
            std::cout << "  2: follow case.json — use Manoeuvring.lowFreqManeuverModel as loaded (no override)\n";
            std::cout << "input model number:\t";
            int modelChoice = 0;
            std::cin >> modelChoice;
            if (modelChoice == 0)
            {
                cfg.Mmg->lowFreqManeuverModel = "mmg_default";
                cfg.Mmg->sjtuMmg.enabled = false;
                std::cout << "[Coupled] runtime: mmg_default (override case.json manoeuvre model).\n";
            }
            else if (modelChoice == 1)
            {
                cfg.Mmg->lowFreqManeuverModel = "sjtu_s175_nomoto";
                cfg.Mmg->sjtuMmg.enabled = true;
                if (cfg.Mmg->sjtuMmg.Lpp <= 1.0e-6)
                    throw std::runtime_error(
                        "SJTU MMG selected but Manoeuvring.CoupledSjtuMmg3DOF is missing or empty in case.json.");
                std::cout << "[Coupled] runtime: sjtu_s175_nomoto (CoupledSjtuMmg3DOF enabled for this run).\n";
            }
            else if (modelChoice == 2)
            {
                std::cout << "[Coupled] runtime: follow case.json lowFreqManeuverModel=\""
                          << cfg.Mmg->lowFreqManeuverModel << "\"\n";
            }
            else
            {
                throw std::runtime_error("Invalid manoeuvre model choice; enter 0, 1, or 2.");
            }

            std::vector<std::shared_ptr<WaveBase>> waves;
            if (cfg.Seakeeping.has_value())
                waves = cfg.Seakeeping->waves;

            CoupledWaveEnvironment waveEnv(waves);
            if (cfg.Coupled.has_value())
                waveEnv.setWaveRampTime(cfg.Coupled->time.waveRampTimeS);

            std::unique_ptr<IWindowSeakeepingSolver> fastSolver;

            // These must outlive the solver — kept in the outer scope.
            std::unique_ptr<CoupledRadiationKernelRepo> radRepo;
            std::unique_ptr<CoupledExcitationKernelRepo> excRepo;
            std::unique_ptr<IWaveForceProvider> waveForceProvider;
            std::shared_ptr<Element> coupledElement;

            if (cfg.Coupled->physics.enableFastSeakeeping)
            {
                if (!cfg.Seakeeping.has_value())
                    throw std::runtime_error("Fast seakeeping is enabled but Seakeeping config is missing.");

                radRepo = std::make_unique<CoupledRadiationKernelRepo>(
                    cfg.Ship,
                    cfg.Seakeeping.value(),
                    casePath);

                // Always build the panel element — needed both for hydrostatics in the
                // TDGF solver and for direct-pressure FK wave force.
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

                // 1) Legacy fitted-table solver — owns its own time integrator.
                if (mode == "fitted")
                {
                    fastSolver = std::make_unique<FittedWaveForceWindowSolver>(
                        cfg.Ship,
                        cfg.Seakeeping.value(),
                        casePath,
                        *radRepo,
                        casePath + "ExcitingForce/FirstOrderWaveForceAmpPhase.csv");
                }
                else
                {
                    // 2) All TDGF-based modes share one solver — only the provider differs.
                    if (mode == "directFK")
                    {
                        waveForceProvider = std::make_unique<DirectPressureFKWaveForceProvider>(
                            cfg.Ship, cfg.Seakeeping.value(), coupledElement);
                    }
                    else if (mode == "directFKdf")
                    {
                        waveForceProvider = std::make_unique<DirectFKPlusDfWaveForceProvider>(
                            cfg.Ship,
                            cfg.Seakeeping.value(),
                            casePath,
                            cfg.Ship.Name,
                            cfg.Coupled->refresh,
                            coupledElement);
                    }
                    else if (mode == "fkImpulse")
                    {
                        waveForceProvider = std::make_unique<FKImpulseWaveForceProvider>(
                            cfg.Seakeeping.value(), casePath, cfg.Ship.Name, cfg.Coupled->refresh);
                    }
                    else if (mode == "impulseKernel")
                    {
                        excRepo = std::make_unique<CoupledExcitationKernelRepo>(
                            cfg.Ship, cfg.Seakeeping.value(), casePath);
                        waveForceProvider = std::make_unique<ImpulseKernelWaveForceProvider>(
                            cfg.Seakeeping.value(), *excRepo, cfg.Coupled->refresh);
                    }
                    else
                    {
                        throw std::runtime_error("Unknown waveForceMode: " + mode);
                    }

                    fastSolver = std::make_unique<TDGFWindowSeakeepingSolver>(
                        cfg.Ship,
                        cfg.Seakeeping.value(),
                        casePath,
                        *radRepo,
                        cfg.Coupled->refresh,
                        *waveForceProvider,
                        coupledElement,
                        cfg.Coupled->physics.enableRadiation);
                }
            }
            else
            {
                fastSolver = std::make_unique<NullWindowSeakeepingSolver>();
            }

            std::unique_ptr<ICoupledLowFreqManeuverSolver> lowFreqSolver =
                makeLowFreqManeuverSolver(cfg.Ship, cfg.Mmg.value());

            std::unique_ptr<ISecondOrderLoadProvider> secondOrderProvider;
            if (cfg.Coupled->physics.enableSecondOrderLoads)
            {
                const std::string driftFolder = casePath + "DriftForce";
                secondOrderProvider = std::make_unique<DriftForceTxtProvider>(driftFolder);
            }
            else
            {
                secondOrderProvider = std::make_unique<NullSecondOrderLoadProvider>();
            }

            TwoTimeScaleCoordinator coordinator(
                cfg,
                casePath,
                *fastSolver,
                *lowFreqSolver,
                *secondOrderProvider,
                waveEnv);

            coordinator.run();

            timer.print_s("Seakeeping-Maneuvering coupling total time: ");
            std::cout << "\n\nRun coupled case successfully!\n" << std::endl;
        }
        return 0;
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << "error: " << e.what() << std::endl;
        return -1;
    }
}