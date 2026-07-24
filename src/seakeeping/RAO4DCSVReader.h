#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <complex>
#include "RAO4DTable.h"

struct RAO4DLoadInfo
{
    bool ok = false;
    std::string message;
    size_t nanCount = 0;
    size_t total = 0;
};

class RAO4DProvider
{
public:
    // 读取：path = ".../RAO.csv"
    RAO4DLoadInfo loadIfExists(const std::string& path)
    {
        RAO4DLoadInfo info;

        if (!RAO4DTable::fileExists(path)) {
            info.ok = false;
            info.message = "RAO.csv not found: " + path;
            return info;
        }

        try {
            table_ = RAO4DTable::readCSV(path);
            loaded_ = true;

            info.ok = true;
            info.message = "RAO.csv loaded: " + path;

            info.nanCount = table_.countNaN();
            info.total = (size_t)table_.nFn * (size_t)table_.nBe * (size_t)table_.nOm * (size_t)table_.nDof;

            if (info.nanCount > 0) {
                info.message += " (warning: grid has NaNs, likely missing cases)";
            }
            return info;
        }
        catch (const std::exception& e) {
            loaded_ = false;
            info.ok = false;
            info.message = std::string("load failed: ") + e.what();
            return info;
        }
    }

    bool loaded() const { return loaded_; }
    const RAO4DTable& table() const { return table_; }

    // 插值：beta 用 DEG
    std::vector<std::complex<double>> evalDeg(double Fn, double betaDeg, double omegaRadPerSec) const
    {
        if (!loaded_) throw std::runtime_error("RAO4DProvider: not loaded");
        return table_.interp(Fn, betaDeg, omegaRadPerSec);
    }

    // 插值：beta 用 RAD
    std::vector<std::complex<double>> evalRad(double Fn, double betaRad, double omegaRadPerSec) const
    {
        constexpr double R2D = 180.0 / 3.141592653589793238462643383279502884;
        return evalDeg(Fn, betaRad * R2D, omegaRadPerSec);
    }

private:
    bool loaded_ = false;
    RAO4DTable table_;
};
