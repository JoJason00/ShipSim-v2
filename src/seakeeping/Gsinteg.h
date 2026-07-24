#pragma once

#include "Greenf.h"
#include "Element.h"
#include "GreenTable.h"
#include "TDGFProvider.h"

struct GreenData
{
    double sG, xdG, ydG, zdG, tdG;
};

class Gsinteg
{
private:
    double U;
    double sG, xdG, ydG, zdG, tdG;

    const std::shared_ptr<Element> element;
    const std::vector<ElementMatrix> element_test;
    const std::vector<Vector3d> point_test;

    void GpMultiply(double&, double&, double&);
    void Gderivative(double&, double&, double&, double&, double&);

public:
    Greenf GF1;

    Gsinteg(const std::shared_ptr<Element> element_s, double u);
    Gsinteg(const std::vector<ElementMatrix>& element,
        const std::vector<Vector3d>& point_test,
        double u);

    // 旧接口：保留，避免破坏其它代码
    GreenData GreenCal(
        const int j,
        const int i,
        const double tn,
        Greenf& GF,
        GreenTable& gGreenTable);

    GreenData GreenCalPanelGauss(
        const int j,
        const int i,
        const double tn,
        Greenf& GF,
        GreenTable& gGreenTable,
        int order = 2);

    std::vector<GreenData> GreenCalPanelGauss_test(
        const double tn,
        Greenf& GF,
        int order = 2);

    GreenData GreenCal_WL_dl(
        int& j,
        int& i,
        double& tn,
        Greenf& GF,
        GreenTable& gGreenTable);

    GreenData GreenCal_WL_Gauss12(
        int& j,
        int& i,
        double& tn,
        Greenf& GF,
        GreenTable& gGreenTable);

    GreenData GreenCal_WL(
        int& j,
        int& i,
        double& tn,
        Greenf& GF,
        GreenTable& gGreenTable);

    // 新接口：由 TDGFProvider 控制 Exact / Interp / ODETable
    GreenData GreenCal(
        const int j,
        const int i,
        const double tn,
        const TDGFProvider& tdgf);

    GreenData GreenCalPanelGauss(
        const int j,
        const int i,
        const double tn,
        const TDGFProvider& tdgf,
        int order = 2);

    GreenData GreenCal_WL_Gauss12(
        int& j,
        int& i,
        double& tn,
        const TDGFProvider& tdgf);

    GreenData GreenCal_WL_Gauss16(
        int& j,
        int& i,
        double& tn,
        const TDGFProvider& tdgf);
};