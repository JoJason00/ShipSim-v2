#include "Wave.h"
#include "WaveBase.h"
#include "../config/WaveConfig.h"
#include "../config/SpectrumParams.h"
#include "../const/Const.h"
#include <Eigen/Dense>
#include <cmath>
#include <sstream>
#include <iomanip>

std::shared_ptr<WaveBase> Wave::createWave(const WaveConfig& cfg)
{
    switch (cfg.type)
    {
        case WaveType::regularwave:
        {
            const RegularWaveConfig& regular = std::get<RegularWaveConfig>(cfg.config);
            return std::make_shared<RegularWave>(regular);
        }

        case WaveType::irregularwave:
        {
            const IrregularWaveConfig& irregular = std::get<IrregularWaveConfig>(cfg.config);
            return std::make_shared<IrregularWave>(irregular);
        }

        case WaveType::crosswave:
        {
            const CrossWaveConfig& cross = std::get<CrossWaveConfig>(cfg.config);
            return std::make_shared<CrossWave>(cross);
        }
    }
    return nullptr;
}

namespace
{
    // Whole-degree heading (waves store direction in radians internally).
    std::string dirDeg(double dirRad)
    {
        return std::to_string(std::lround(dirRad * 180.0 / PI));
    }

    // Fixed-precision number with trailing zeros trimmed: 0.70 -> "0.7".
    std::string shortNum(double x, int prec = 2)
    {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(prec) << x;
        std::string s = ss.str();
        if (s.find('.') != std::string::npos)
        {
            s.erase(s.find_last_not_of('0') + 1);
            if (!s.empty() && s.back() == '.') s.pop_back();
        }
        return s;
    }
}

std::string Wave::systemDescriptor(const std::shared_ptr<WaveBase>& w)
{
    if (auto rw = std::dynamic_pointer_cast<RegularWave>(w))
        return "Reg-d" + dirDeg(rw->direction()) + "-w" + shortNum(rw->getFreq());

    if (auto ir = std::dynamic_pointer_cast<IrregularWave>(w))
    {
        const SpectrumParams& s = ir->spectrum();
        const std::string d = "-d" + dirDeg(ir->direction());
        // Seed in the tag so repeated runs of the SAME spectrum under different
        // random phase realisations (different wave patterns) get distinct
        // output folders instead of overwriting each other.
        const std::string seedTag = "-s" + std::to_string(s.seed);
        switch (s.type)
        {
        case SpectrumType::JONSWAP:
            return "JON" + d + "-Hs" + shortNum(s.Hs) + "-Tp" + shortNum(s.Tp)
                 + (std::abs(s.gamma - 3.3) > 1e-6 ? "-g" + shortNum(s.gamma, 1)
                                                   : std::string())
                 + seedTag;
        case SpectrumType::PM:
            return "PM" + d + "-Uw" + shortNum(s.Uwind, 1) + seedTag;
        case SpectrumType::ITTC:
            return "ITTC" + d + "-Hs" + shortNum(s.Hs) + "-T1" + shortNum(s.T1) + seedTag;
        case SpectrumType::OH:
        {
            // bimodal: "x" joins the two systems' Hs; lo/hi flags ohMode.
            const std::string m = (s.ohMode == 1) ? "-lo"
                                 : (s.ohMode == 2) ? "-hi"
                                                   : std::string();
            return "OH" + d + "-Hs" + shortNum(s.Hs1) + "x" + shortNum(s.Hs2) + m + seedTag;
        }
        }
    }
    return "Wave";
}

std::string Wave::conditionTag(const std::shared_ptr<WaveBase>& w)
{
    if (auto cw = std::dynamic_pointer_cast<CrossWave>(w))
        return "X_" + systemDescriptor(cw->wave1())
             + "__" + systemDescriptor(cw->wave2());
    return systemDescriptor(w);
}
