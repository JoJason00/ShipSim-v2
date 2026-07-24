#pragma once

#include <string>
#include <vector>
#include <array>

#include "Greenf.h"
#include "GreenTable.h"

struct TDGFValue
{
    double F = 0.0;   // Greenf.Gf
    double Ft = 0.0;  // Greenf.Gbd = dF / d tau
    double Fm = 0.0;  // Greenf.Gmd = dF / d mu
};

struct TDGFNode
{
    double F = 0.0;
    double Ft = 0.0;
    double Ftt = 0.0;
    double Fttt = 0.0;
};

class TDGFProvider
{
public:
    void initExact();

    void initInterp(GreenTable* table);

    void initODETable(
        const std::string& tablePath,
        double rk4Step);

    TDGFValue eval(double tau, double mu) const;

private:
    enum class Method
    {
        Exact,
        Interp,
        ODETable
    };

    Method method_ = Method::Exact;

    GreenTable* interpTable_ = nullptr;

    double rk4Step_ = 0.001;

    std::vector<double> tauNodes_;
    std::vector<double> muNodes_;
    std::vector<TDGFNode> table_;

private:
    TDGFValue evalExact(double tau, double mu) const;
    TDGFValue evalInterp(double tau, double mu) const;
    TDGFValue evalODETable(double tau, double mu) const;

    bool loadODETable(const std::string& path);

    int findLeft(const std::vector<double>& xs, double x) const;

    TDGFNode nodeAt(int it, int im) const;

    TDGFNode interpMu(int it, double mu) const;

    TDGFNode advanceODE(
        TDGFNode node,
        double tau0,
        double tau1,
        double mu) const;

    static std::array<double, 4> rhs(
        double tau,
        double mu,
        const std::array<double, 4>& y);

    void checkODETableRange(double tau, double mu) const;
};