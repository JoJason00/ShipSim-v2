#pragma once

#include <vector>
#include <string>
#include <memory>
#include <Eigen/Dense>
#include "../config/SeakeepingConfig.h"
#include "../config/CaseConfig.h"

class Element;

constexpr int MODE_SURGE = 0;
constexpr int MODE_SWAY = 1;
constexpr int MODE_HEAVE = 2;
constexpr int MODE_ROLL = 3;
constexpr int MODE_PITCH = 4;
constexpr int MODE_YAW = 5;

struct HydrostaticsData {
    double Awp{ 0.0 };
    double xF{ 0.0 };
    double yF{ 0.0 };
    double IT{ 0.0 };
    double IL{ 0.0 };

    double c33{ 0.0 }, c44{ 0.0 }, c55{ 0.0 };
    double c35{ 0.0 }, c53{ 0.0 };
    double c34{ 0.0 }, c43{ 0.0 };
    double c54{ 0.0 }, c45{ 0.0 };

    double Vdisp = 0.0;
    double xB = 0.0;
    double yB = 0.0;
    double zB = 0.0;
};

struct CaseContextLite {
    double Amp{ 0 };
    double W{ 0 };
    double U{ 0 };
    double dt{ 0 };
    double we{ 0 };
};

class SeakeepingDOF {
public:
    static int findModeIndex(const SeakeepingConfig& cfg, int modeId);
    static std::vector<std::string> modeNames(const SeakeepingConfig& cfg);

    static HydrostaticsData defaultHydrostatics(const SeakeepingConfig& cfg);

    static HydrostaticsData hydrostaticsFromWaterline(
        const ShipConfig& Ship, const SeakeepingConfig& Seakeeping,
        const Element& element, double tolRel = 5e-4
    );

    static void buildSystemMatrices(
        const ShipConfig& Ship, const SeakeepingConfig& Seakeeping,
        const HydrostaticsData& hs,
        Eigen::MatrixXd& M,
        Eigen::MatrixXd& B,
        Eigen::MatrixXd& C);

    static void buildRadiationVn(const SeakeepingConfig& cfg,
        const std::shared_ptr<Element> element,
        const std::vector<double>& y,
        double U,
        Eigen::VectorXd& rVn);

    static void scaleMotionsForOutput(const SeakeepingConfig& cfg,
        const CaseContextLite& ctx,
        Eigen::MatrixXd& motionsInOut);

    static double rollDampingFromZeta(const ShipConfig& cfg,
        const HydrostaticsData& hs,
        double zeta);

    static void DumpWaterlineDebugCSV(
        const ShipConfig& Ship, const SeakeepingConfig& Seakeeping,
        const Element& element,
        const std::string& prefix,
        double tolRel = 1e-5);

    static void computeBuoyancyFromSurfacePanels(
        const Element& element,
        HydrostaticsData& hs,
        double zWL = 0.0);
};
