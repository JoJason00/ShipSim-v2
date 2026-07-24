#include "CoupledLowFreqManeuverFactory.h"
#include "CoupledMmg3DOFCore.h"
#include "CoupledSjtuMmg3DOFCore.h"

#include <cctype>
#include <stdexcept>

namespace
{
    std::string toLower(std::string s)
    {
        for (char& c : s)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    }
}

std::string lowFreqManeuverModelTag(const MmgConfig& mmg)
{
    const std::string m = toLower(mmg.lowFreqManeuverModel);
    if (m.empty() || m == "mmg_default" || m == "default" || m == "coupled_mmg")
        return "mmg_default";
    if (m == "sjtu_s175_nomoto" || m == "sjtu" || m == "sjtu_mmg")
        return "sjtu_s175_nomoto";
    return mmg.lowFreqManeuverModel;
}

std::string lowFreqManeuverOutputFolderSuffix(const MmgConfig& mmg)
{
    const std::string t = lowFreqManeuverModelTag(mmg);
    if (t.empty() || t == "mmg_default")
        return "";
    if (t == "sjtu_s175_nomoto")
        return "_sjtu";
    std::string s = "_";
    for (char c : t)
    {
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_')
            s += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        else if (c == '-' || c == '.')
            s += '_';
    }
    return s.size() > 1 ? s : "";
}

std::unique_ptr<ICoupledLowFreqManeuverSolver> makeLowFreqManeuverSolver(
    const ShipConfig& ship,
    const MmgConfig& mmg)
{
    const std::string m = toLower(mmg.lowFreqManeuverModel);
    if (m.empty() || m == "mmg_default" || m == "default" || m == "coupled_mmg")
        return std::make_unique<CoupledMmg3DOFCore>(ship, mmg);

    if (m == "sjtu_s175_nomoto" || m == "sjtu" || m == "sjtu_mmg")
    {
        if (!mmg.sjtuMmg.enabled)
            throw std::runtime_error(
                "lowFreqManeuverModel selects SJTU MMG but Manoeuvring.CoupledSjtuMmg3DOF.enabled "
                "is false or block missing. Set enabled=true and fill coefficients.");
        return std::make_unique<CoupledSjtuMmg3DOFCore>(ship, mmg);
    }

    throw std::runtime_error("Unknown Manoeuvring.lowFreqManeuverModel: \"" + mmg.lowFreqManeuverModel + "\"");
}
