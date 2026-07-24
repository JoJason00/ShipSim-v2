#include "TDGFProvider.h"

#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>

void TDGFProvider::initExact()
{
    method_ = Method::Exact;
}

void TDGFProvider::initInterp(GreenTable* table)
{
    if (!table)
        throw std::runtime_error("TDGFProvider::initInterp: null GreenTable.");

    method_ = Method::Interp;
    interpTable_ = table;
}

void TDGFProvider::initODETable(
    const std::string& tablePath,
    double rk4Step)
{
    if (rk4Step <= 0.0)
        throw std::runtime_error("TDGFProvider::initODETable: RK4Step must be positive.");

    if (!loadODETable(tablePath))
        throw std::runtime_error("TDGFProvider: cannot load ODE table: " + tablePath);

    rk4Step_ = rk4Step;
    method_ = Method::ODETable;
}

TDGFValue TDGFProvider::eval(double tau, double mu) const
{
    if (tau < 0.0)
    {
        throw std::runtime_error(
            "TDGFProvider::eval: tau/beta must be positive. tau = "
            + std::to_string(tau));
    }

    if (mu < 0.0 || mu > 1.0)
    {
        throw std::runtime_error(
            "TDGFProvider::eval: mu must be in [0, 1]. mu = "
            + std::to_string(mu));
    }

    switch (method_)
    {
    case Method::Exact:
    {
        // Exact 方法不依赖离线表。
        // 为避免 mu=0 处数值奇异，内部轻微抬高。
        //const double tauSafe = std::max(tau, 1.0e-8);
        //const double muSafe = std::max(mu, 1.0e-8);

        return evalExact(tau, mu);
    }

    case Method::Interp:
    {
        // Interp 使用现有 GreenTable，是否越界交给 GreenTable 自己判断。
        return evalInterp(tau, mu);
    }

    case Method::ODETable:
    {
        // ODETable 必须严格检查表格范围，不允许自动截断。
        checkODETableRange(tau, mu);
        return evalODETable(tau, mu);
    }
    }

    throw std::runtime_error("TDGFProvider::eval: unknown method.");
}

TDGFValue TDGFProvider::evalExact(double tau, double mu) const
{
    Greenf gf;
    gf.GreenFunctionCal(tau, mu);

    TDGFValue v;
    v.F = gf.Gf;
    v.Ft = gf.Gbd;
    v.Fm = gf.Gmd;
    return v;
}

TDGFValue TDGFProvider::evalInterp(double tau, double mu) const
{
    if (!interpTable_)
        throw std::runtime_error("TDGFProvider::evalInterp: table not initialized.");

    TDGFValue v;
    interpTable_->eval(tau, mu, v.F, v.Ft, v.Fm);
    return v;
}

int TDGFProvider::findLeft(
    const std::vector<double>& xs,
    double x) const
{
    if (xs.size() < 2)
        throw std::runtime_error("TDGFProvider::findLeft: bad node vector.");

    if (x <= xs.front())
        return 0;

    if (x >= xs.back())
        return static_cast<int>(xs.size()) - 2;

    auto it = std::upper_bound(xs.begin(), xs.end(), x);
    return static_cast<int>(std::distance(xs.begin(), it)) - 1;
}

TDGFNode TDGFProvider::nodeAt(int it, int im) const
{
    const int nMu = static_cast<int>(muNodes_.size());
    return table_[static_cast<std::size_t>(it * nMu + im)];
}

TDGFNode TDGFProvider::interpMu(
    int it,
    double mu) const
{
    const int im = findLeft(muNodes_, mu);

    const double m0 = muNodes_[im];
    const double m1 = muNodes_[im + 1];

    const double w = (std::abs(m1 - m0) < 1.0e-14)
        ? 0.0
        : (mu - m0) / (m1 - m0);

    TDGFNode a = nodeAt(it, im);
    TDGFNode b = nodeAt(it, im + 1);

    TDGFNode r;
    r.F = a.F * (1.0 - w) + b.F * w;
    r.Ft = a.Ft * (1.0 - w) + b.Ft * w;
    r.Ftt = a.Ftt * (1.0 - w) + b.Ftt * w;
    r.Fttt = a.Fttt * (1.0 - w) + b.Fttt * w;

    return r;
}

std::array<double, 4> TDGFProvider::rhs(
    double tau,
    double mu,
    const std::array<double, 4>& y)
{
    std::array<double, 4> dy{};

    dy[0] = y[1];
    dy[1] = y[2];
    dy[2] = y[3];

    dy[3] =
        -9.0 / 4.0 * y[0]
        - 7.0 / 4.0 * tau * y[1]
        - (4.0 * mu + tau * tau / 4.0) * y[2]
        - mu * tau * y[3];

    return dy;
}

TDGFNode TDGFProvider::advanceODE(
    TDGFNode node,
    double tau0,
    double tau1,
    double mu) const
{
    if (tau1 <= tau0)
        return node;

    std::array<double, 4> y{
        node.F,
        node.Ft,
        node.Ftt,
        node.Fttt
    };

    double tau = tau0;

    while (tau < tau1)
    {
        const double h = std::min(rk4Step_, tau1 - tau);

        auto k1 = rhs(tau, mu, y);

        std::array<double, 4> yt{};
        for (int i = 0; i < 4; ++i)
            yt[i] = y[i] + 0.5 * h * k1[i];

        auto k2 = rhs(tau + 0.5 * h, mu, yt);

        for (int i = 0; i < 4; ++i)
            yt[i] = y[i] + 0.5 * h * k2[i];

        auto k3 = rhs(tau + 0.5 * h, mu, yt);

        for (int i = 0; i < 4; ++i)
            yt[i] = y[i] + h * k3[i];

        auto k4 = rhs(tau + h, mu, yt);

        for (int i = 0; i < 4; ++i)
        {
            y[i] += h * (
                k1[i]
                + 2.0 * k2[i]
                + 2.0 * k3[i]
                + k4[i]
                ) / 6.0;
        }

        tau += h;
    }

    node.F = y[0];
    node.Ft = y[1];
    node.Ftt = y[2];
    node.Fttt = y[3];

    return node;
}


TDGFValue TDGFProvider::evalODETable(
    double tau,
    double mu) const
{
    if (tauNodes_.empty() || muNodes_.empty() || table_.empty())
    {
        throw std::runtime_error(
            "TDGFProvider::evalODETable: table is empty.");
    }

    // 这里再次检查，防止有人绕过 eval() 直接调用 evalODETable()
    checkODETableRange(tau, mu);

    const int it = findLeft(tauNodes_, tau);
    const double tau0 = tauNodes_[it];

    TDGFNode base = interpMu(it, mu);
    TDGFNode y = advanceODE(base, tau0, tau, mu);

    TDGFValue v;
    v.F = y.F;
    v.Ft = y.Ft;

    // 第一版简化：Fmu 用有限差分。
    // 注意有限差分点也不能越出表格范围。
    const double dm = 1.0e-4;

    const double muMin = muNodes_.front();
    const double muMax = muNodes_.back();

    double m0 = mu - dm;
    double m1 = mu + dm;

    if (m0 < muMin)
    {
        m0 = mu;
        m1 = mu + dm;

        if (m1 > muMax)
        {
            throw std::runtime_error(
                "TDGFProvider::evalODETable: cannot compute Fmu near lower mu boundary. "
                "mu = " + std::to_string(mu)
                + ", table range = ["
                + std::to_string(muMin) + ", "
                + std::to_string(muMax) + "].");
        }
    }
    else if (m1 > muMax)
    {
        m1 = mu;
        m0 = mu - dm;

        if (m0 < muMin)
        {
            throw std::runtime_error(
                "TDGFProvider::evalODETable: cannot compute Fmu near upper mu boundary. "
                "mu = " + std::to_string(mu)
                + ", table range = ["
                + std::to_string(muMin) + ", "
                + std::to_string(muMax) + "].");
        }
    }

    TDGFNode y1 = advanceODE(interpMu(it, m1), tau0, tau, m1);
    TDGFNode y0 = advanceODE(interpMu(it, m0), tau0, tau, m0);

    v.Fm = (y1.F - y0.F) / (m1 - m0);

    return v;
}



bool TDGFProvider::loadODETable(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in)
        return false;

    uint32_t magic = 0;
    uint32_t version = 0;
    uint64_t nt = 0;
    uint64_t nm = 0;

    in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    in.read(reinterpret_cast<char*>(&nt), sizeof(nt));
    in.read(reinterpret_cast<char*>(&nm), sizeof(nm));

    if (magic != 0x54444746 || version != 1 || nt == 0 || nm == 0)
        return false;

    tauNodes_.resize(static_cast<std::size_t>(nt));
    muNodes_.resize(static_cast<std::size_t>(nm));
    table_.resize(static_cast<std::size_t>(nt * nm));

    in.read(reinterpret_cast<char*>(tauNodes_.data()),
        static_cast<std::streamsize>(tauNodes_.size() * sizeof(double)));

    in.read(reinterpret_cast<char*>(muNodes_.data()),
        static_cast<std::streamsize>(muNodes_.size() * sizeof(double)));

    in.read(reinterpret_cast<char*>(table_.data()),
        static_cast<std::streamsize>(table_.size() * sizeof(TDGFNode)));

    return static_cast<bool>(in);

    if (in && !tauNodes_.empty() && !muNodes_.empty())
    {
        std::cout << "[TDGFProvider] ODE table loaded.\n"
            << "  tau/beta range: ["
            << tauNodes_.front() << ", "
            << tauNodes_.back() << "]\n"
            << "  mu range: ["
            << muNodes_.front() << ", "
            << muNodes_.back() << "]\n"
            << "  tau nodes: " << tauNodes_.size() << "\n"
            << "  mu nodes: " << muNodes_.size() << "\n";
    }
}


void TDGFProvider::checkODETableRange(double tau, double mu) const
{
    if (tauNodes_.empty() || muNodes_.empty())
    {
        throw std::runtime_error(
            "TDGFProvider::checkODETableRange: ODE table is empty.");
    }

    const double tauMin = tauNodes_.front();
    const double tauMax = tauNodes_.back();
    const double muMin = muNodes_.front();
    const double muMax = muNodes_.back();

    if (tau < tauMin || tau > tauMax)
    {
        throw std::runtime_error(
            "TDGFProvider::ODETable tau/beta out of range. "
            "tau = " + std::to_string(tau)
            + ", table range = ["
            + std::to_string(tauMin) + ", "
            + std::to_string(tauMax) + "]. "
            "Please rebuild the TDGF ODE table with a larger tauMax.");
    }

    if (mu < muMin || mu > muMax)
    {
        throw std::runtime_error(
            "TDGFProvider::ODETable mu out of range. "
            "mu = " + std::to_string(mu)
            + ", table range = ["
            + std::to_string(muMin) + ", "
            + std::to_string(muMax) + "]. "
            "Please rebuild the TDGF ODE table with wider mu range.");
    }
}