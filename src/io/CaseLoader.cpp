#include "CaseLoader.h"
#include "../const/Const.h"
#include <../json/json.h>  
#include <string>
#include <fstream>
#include <stdexcept>
#include <cctype>
#include <iostream>
#include <iomanip>
#include <wave/WaveBase.h>
#include <wave/RegularWave.h>
#include <wave/IrregularWave.h>
#include <wave/CrossWave.h>

namespace {
    WaveType parseWaveType(const std::string& s)
    {
        if (s == "RegularWave")   return WaveType::regularwave;
        if (s == "IrregularWave") return WaveType::irregularwave;
        if (s == "CrossWave")     return WaveType::crosswave;

        throw std::runtime_error("Unknown wave type: " + s);
    }

    SpectrumType parseSpectrumType(const std::string& s)
    {
        if (s == "JONSWAP") return SpectrumType::JONSWAP;
        if (s == "PM")      return SpectrumType::PM;
        if (s == "ITTC")    return SpectrumType::ITTC;
        if (s == "OH")      return SpectrumType::OH;

        throw std::runtime_error("Unknown spectrum type: " + s);
    }

    static FKForceMethod parseFKForceMethod(const std::string& s)
    {
        if (s == "PrescribedAmp")     return FKForceMethod::PrescribedAmp;
        if (s == "TDGFImpulse")       return FKForceMethod::TDGFImpulse;
        if (s == "DirectPressure")    return FKForceMethod::DirectPressure;
        // Back-compat: merged into DirectPressure (numerically identical).
        if (s == "DirectFKImpulseDf") return FKForceMethod::DirectPressure;

        throw std::runtime_error("Unknown FKForceMethod: " + s);
    }

    template <typename T>
    T read(const Json::Value& node, const std::string& name);

    template <>
    std::string read<std::string>(const Json::Value& node, const std::string& name)
    {
        if (!node.isMember(name) || !node[name].isString())
        {
            throw std::runtime_error("there is no:" + name + "in json file or " + name + "is not string");
        }
        return node[name].asString();
    }

    template <>
    int read<int>(const Json::Value& node, const std::string& name)
    {
        if (!node.isMember(name) || !node[name].isInt())
        {
            throw std::runtime_error("there is no " + name + " in json file or " + name + " is not int");
        }
        return node[name].asInt();
    }

    template <>
    double read<double>(const Json::Value& node, const std::string& name)
    {
        if (!node.isMember(name) || !node[name].isNumeric())
        {
            throw std::runtime_error("there is no " + name + " in json file or " + name + " is not numeric");
        }
        return node[name].asDouble();
    }

    template <>
    bool read<bool>(const Json::Value& node, const std::string& name)
    {
        if (!node.isMember(name) || !node[name].isBool())
        {
            throw std::runtime_error("there is no " + name + " in json file or " + name + " is not bool");
        }
        return node[name].asBool();
    }

    template <>
    std::vector<int> read<std::vector<int>>(const Json::Value& node, const std::string& name)
    {
        if (!node.isMember(name) || !node[name].isArray())
        {
            throw std::runtime_error("there is no " + name + " in json file or " + name + " is not array");
        }
        std::vector<int> vec;
        const Json::Value& arr = node[name];
        for (int i = 0; i < arr.size(); ++i)
        {
            if (!arr[i].isInt())
            {
                throw std::runtime_error("element " + std::to_string(i) + " of " + name + " is not int");
            }
            vec.push_back(arr[i].asInt());
        }
        return vec;
    }

    template <>
    std::vector<double> read<std::vector<double>>(const Json::Value& node, const std::string& name)
    {
        if (!node.isMember(name) || !node[name].isArray())
        {
            throw std::runtime_error("there is no " + name + " in json file or " + name + " is not array");
        }
        std::vector<double> vec;
        const Json::Value& arr = node[name];
        for (int i = 0; i < arr.size(); ++i)
        {
            if (!arr[i].isNumeric())
            {
                throw std::runtime_error("element " + std::to_string(i) + " of " + name + " is not numeric");
            }
            vec.push_back(arr[i].asDouble());
        }
        return vec;
    }


    template <>
    Json::Value read<Json::Value>(const Json::Value& node, const std::string& name)
    {
        if (!node.isMember(name) || !node[name].isObject())
        {
            throw std::runtime_error("there is no " + name + " in json file or " + name + " is not Object");
        }
        return node[name];
    }



    std::vector<std::shared_ptr<WaveBase>> loadWaves(const Json::Value& node)
    {
        std::vector<std::shared_ptr<WaveBase>> waves;

        std::string wave_type_str = read<std::string>(node, "type");
        WaveType type = parseWaveType(wave_type_str);

        switch (type)
        {
        case WaveType::regularwave:
        {
            Json::Value regular_node = read<Json::Value>(node, "RegularWave");

            double H = read<double>(regular_node, "H");

            // Optional initial phase ε [rad]; default 0. Per (sub-)wave, so a
            // CrossWave of two RegularWave blocks can set each independently.
            double phase0 = 0.0;
            if (regular_node.isMember("phase0"))
                phase0 = read<double>(regular_node, "phase0");

            auto directions = read<std::vector<double>>(regular_node, "direction");
            for (double& d : directions)
            {
                d *= (2.0 * PI / 360.0);
            }

            std::string set = read<std::string>(regular_node, "set");
            std::vector<double> omegas;

            if (set == "nada")
            {
                auto nada_list = read<std::vector<double>>(regular_node, "nada");
                for (double nada : nada_list)
                    omegas.push_back(sqrt(2.0 * PI * G / nada));
            }
            else if (set == "nada/L")
            {
                auto nada_list = read<std::vector<double>>(regular_node, "nada/L");
                double nada_i;
				double L = read<double>(regular_node, "L");
                for (double nada : nada_list)
                {
					nada_i = nada * L;
                    omegas.push_back(sqrt(2.0 * PI * G / nada_i));
                }
            }
            else if (set == "omiga")
                omegas = read<std::vector<double>>(regular_node, "omiga");
            else
                throw std::runtime_error("unknow set of regular wave: " + set);

            for (double dir : directions)
            {
                for (double w : omegas)
                {
                    RegularWaveConfig reg(dir, H, w, phase0);
                    waves.push_back(std::make_shared<RegularWave>(reg));

                    /*std::cout << "std::make_shared<RegularWave>(reg): " << std::make_shared<RegularWave>(reg)->direction() << std::endl;
                    std::cout << "std::make_shared<RegularWave>(reg): " << std::make_shared<RegularWave>(reg)->getFreq() << std::endl;

                    std::cout << " waves: " << waves.front()->getFreq() << std::endl;*/
                    //std::cout << " waves: " << std::make_shared<RegularWave>(reg)->getFreq() << std::endl;
                }
            }
            break;
        }
        case WaveType::irregularwave:
        {
            // ----------------------------------------------------------------
            // "set"-style batch expansion (back-compatible with the old scalar
            // form, which is just the N=1 special case):
            //   * Per-CASE spectrum parameters (ITTC: Hs,T1 / JONSWAP: Hs,Tp,..
            //     / PM: Uwind / OH: Hs1,Hs2,wm1,wm2,lam1,lam2) are PAIRED equal-
            //     length arrays — element i defines case i. omegaLo/omegaHi/
            //     nComponents may be scalar (shared) or per-case arrays.
            //   * "direction" (deg) and "seed" are CROSS-multiplied.
            // Total waves = nCase * nDir * nSeed. A plain scalar field is just a
            // length-1 array, so a fully scalar block yields exactly one wave —
            // identical to the previous behaviour.
            // ----------------------------------------------------------------
            Json::Value irr = read<Json::Value>(node, "IrregularWave");

            // scalar-or-array readers (missing key -> empty vector)
            auto vecD = [](const Json::Value& nd, const std::string& key)
                -> std::vector<double>
            {
                std::vector<double> v;
                if (!nd.isMember(key)) return v;
                const Json::Value& x = nd[key];
                if (x.isArray())
                    for (Json::ArrayIndex i = 0; i < x.size(); ++i) v.push_back(x[i].asDouble());
                else if (x.isNumeric())
                    v.push_back(x.asDouble());
                return v;
            };
            auto vecI = [](const Json::Value& nd, const std::string& key)
                -> std::vector<int>
            {
                std::vector<int> v;
                if (!nd.isMember(key)) return v;
                const Json::Value& x = nd[key];
                if (x.isArray())
                    for (Json::ArrayIndex i = 0; i < x.size(); ++i) v.push_back(x[i].asInt());
                else if (x.isInt())
                    v.push_back(x.asInt());
                return v;
            };

            SpectrumType stype = parseSpectrumType(read<std::string>(irr, "SpectrumType"));

            // discretisation controls (scalar or per-case)
            std::vector<double> omLo = vecD(irr, "omegaLo");
            std::vector<double> omHi = vecD(irr, "omegaHi");
            std::vector<double> nCmp = vecD(irr, "nComponents");
            const bool randVal = irr.isMember("randomFreqInBand")
                ? read<bool>(irr, "randomFreqInBand") : true;

            // cross-multiplied axes
            std::vector<double> dirsDeg = vecD(irr, "direction");
            if (dirsDeg.empty())
                throw std::runtime_error("IrregularWave: 'direction' is required");
            std::vector<int> seeds = vecI(irr, "seed");
            if (seeds.empty()) seeds.push_back(1);

            // per-case spectrum parameter arrays (paired)
            std::vector<double> Hs   = vecD(irr, "Hs");
            std::vector<double> T1   = vecD(irr, "T1");
            std::vector<double> Tp   = vecD(irr, "Tp");
            std::vector<double> gam  = vecD(irr, "gamma");
            std::vector<double> Uw   = vecD(irr, "Uwind");
            std::vector<double> Hs1  = vecD(irr, "Hs1");
            std::vector<double> Hs2  = vecD(irr, "Hs2");
            std::vector<double> wm1  = vecD(irr, "wm1");
            std::vector<double> wm2  = vecD(irr, "wm2");
            std::vector<double> lam1 = vecD(irr, "lam1");
            std::vector<double> lam2 = vecD(irr, "lam2");
            std::vector<int>    ohM  = vecI(irr, "ohMode");

            // number of paired cases = longest defining array for this type
            std::size_t nCase = 1;
            auto bump = [&nCase](const std::vector<double>& v)
            { if (v.size() > nCase) nCase = v.size(); };
            switch (stype)
            {
            case SpectrumType::ITTC:    bump(Hs); bump(T1); break;
            case SpectrumType::JONSWAP: bump(Hs); bump(Tp); break;
            case SpectrumType::PM:      bump(Uw); break;
            case SpectrumType::OH:
                bump(Hs1); bump(Hs2); bump(wm1); bump(wm2); bump(lam1); bump(lam2);
                break;
            }

            auto checkLen = [nCase](const std::vector<double>& v, const char* nm)
            {
                if (!v.empty() && v.size() != 1 && v.size() != nCase)
                    throw std::runtime_error(std::string("IrregularWave: array '")
                        + nm + "' length must be 1 or equal to the number of cases ("
                        + std::to_string(nCase) + ")");
            };
            checkLen(Hs, "Hs");   checkLen(T1, "T1");   checkLen(Tp, "Tp");
            checkLen(Uw, "Uwind"); checkLen(Hs1, "Hs1"); checkLen(Hs2, "Hs2");
            checkLen(wm1, "wm1"); checkLen(wm2, "wm2"); checkLen(lam1, "lam1");
            checkLen(lam2, "lam2"); checkLen(omLo, "omegaLo");
            checkLen(omHi, "omegaHi"); checkLen(nCmp, "nComponents");

            auto at = [](const std::vector<double>& v, std::size_t i, double dflt)
                -> double
            { if (v.empty()) return dflt; return v.size() == 1 ? v.front() : v.at(i); };
            auto atI = [](const std::vector<int>& v, std::size_t i, int dflt)
                -> int
            { if (v.empty()) return dflt; return v.size() == 1 ? v.front() : v.at(i); };

            for (std::size_t c = 0; c < nCase; ++c)
            {
                SpectrumParams sp;
                sp.type = stype;
                sp.omegaLo = at(omLo, c, sp.omegaLo);
                sp.omegaHi = at(omHi, c, sp.omegaHi);
                sp.nComponents = static_cast<int>(at(nCmp, c, sp.nComponents));
                sp.randomFreqInBand = randVal;

                switch (stype)
                {
                case SpectrumType::PM:
                    sp.Uwind = at(Uw, c, 0.0);
                    break;
                case SpectrumType::ITTC:
                    sp.Hs = at(Hs, c, 0.0);
                    sp.T1 = at(T1, c, 0.0);
                    break;
                case SpectrumType::JONSWAP:
                    sp.Hs = at(Hs, c, 0.0);
                    sp.Tp = at(Tp, c, 0.0);
                    sp.gamma = at(gam, c, 3.3);
                    break;
                case SpectrumType::OH:
                    sp.Hs1 = at(Hs1, c, 0.0);
                    sp.Hs2 = at(Hs2, c, 0.0);
                    sp.wm1 = at(wm1, c, 0.0);
                    sp.wm2 = at(wm2, c, 0.0);
                    sp.lam1 = at(lam1, c, 0.0);
                    sp.lam2 = at(lam2, c, 0.0);
                    sp.ohMode = atI(ohM, c, 0);
                    break;
                }

                for (double dDeg : dirsDeg)
                {
                    for (int sd : seeds)
                    {
                        IrregularWaveConfig irregular{};
                        irregular.direction = dDeg * (2.0 * PI / 360.0);
                        irregular.spectrum = sp;
                        irregular.spectrum.seed = static_cast<std::uint64_t>(sd);
                        waves.push_back(std::make_shared<IrregularWave>(irregular));
                    }
                }
            }
            break;
        }
        case WaveType::crosswave:
        {
            // "CrossWave": { "wave1": {...}, "wave2": {...}}
            Json::Value cross_node = read<Json::Value>(node, "CrossWave");

            auto w1_list = loadWaves(read<Json::Value>(cross_node, "wave1"));
            auto w2_list = loadWaves(read<Json::Value>(cross_node, "wave2"));
            if (w1_list.empty() || w2_list.empty())
                throw std::runtime_error("CrossWave needs wave1 and wave2");

            CrossWaveConfig cross{};

            // Optional ship start position [x0, y0] in the Earth frame [m].
            // Only meaningful for a crossing sea: it shifts the two sub-waves'
            // phases by different amounts, changing their relative phase at the
            // ship. Absent / [0,0] -> bit-identical to before.
            if (cross_node.isMember("start_position"))
            {
                auto sp = read<std::vector<double>>(cross_node, "start_position");
                if (sp.size() != 2)
                    throw std::runtime_error(
                        "CrossWave.start_position must be [x0, y0] (2 values)");
                cross.startX = sp[0];
                cross.startY = sp[1];
            }

            for (auto& w1 : w1_list)
            {
                for (auto& w2 : w2_list)
                {
                    cross.wave1 = w1;
                    cross.wave2 = w2;
                    waves.push_back(std::make_shared<CrossWave>(cross));
                }
            }
            break;
        }
        default:
            throw std::runtime_error("the wave type is not support\n");
        }

        // Per-wave-entry debug switch (uniform across all wave types): dump the
        // incident wave-elevation time history at the fixed body reference
        // point. Read from this entry's node so it applies to the wave the
        // ship actually sees (the CrossWave object, not its sub-blocks).
        if (node.isMember("output_history"))
        {
            const bool oh = read<bool>(node, "output_history");
            for (auto& w : waves)
                w->outputHistory = oh;
        }

        if (!waves.empty())
            std::cout << "waves:  " << waves.front()->getFreq() << "\t" << waves.front()->getAmp() << std::endl;
        return waves;
    }


    ShipConfig loadShipConfig(const Json::Value& shipNode)
    {
        ShipConfig ship;

        //name
        ship.Name = read<std::string>(shipNode, "Name");

        // Geometry
        Json::Value geo = read<Json::Value>(shipNode, "Geometry");
        ship.Geometry.Length = read<double>(geo, "Length");
        ship.Geometry.Draft = read<double>(geo, "Draft");
        ship.Geometry.Displacement = read<double>(geo, "Displacement");
        ship.Geometry.Breadth = read<double>(geo, "Breadth");
        ship.Geometry.CB = read<double>(geo, "CB");
        ship.Geometry.Trim = read<double>(geo, "Trim");

        // Mass
        Json::Value mass = read<Json::Value>(shipNode, "Mass");
        ship.Mass.Mass = read<double>(mass, "Mass");
        ship.Mass.Ixx = read<double>(mass, "Ixx");
        ship.Mass.Iyy = read<double>(mass, "Iyy");
        ship.Mass.Izz = read<double>(mass, "Izz");
        ship.Mass.CG = read<std::vector<double>>(mass, "CG");
        ship.Mass.GM = read<double>(mass, "GM");

        return ship;
    }


    static RollDampingMode parseRollDampingMode(const std::string& s)
    {
        if (s == "DirectCoefficients") return RollDampingMode::DirectCoefficients;
        if (s == "FromDecayCsv")       return RollDampingMode::FromDecayCsv;
        throw std::runtime_error("Unknown RollDamping mode: " + s);
    }

    static RollDampingConfig loadRollDampingConfig(const Json::Value& node)
    {
        RollDampingConfig cfg;

        // 不写 RollDamping 就默认关闭
        if (!node.isMember("RollDamping") || !node["RollDamping"].isObject())
        {
            cfg.enabled = false;
            return cfg;
        }

        Json::Value rollNode = read<Json::Value>(node, "RollDamping");

        // enabled 可选，默认 true
        if (rollNode.isMember("enabled"))
            cfg.enabled = read<bool>(rollNode, "enabled");
        else
            cfg.enabled = true;

        // mode 必填
        cfg.mode = parseRollDampingMode(read<std::string>(rollNode, "mode"));

        if (cfg.mode == RollDampingMode::DirectCoefficients)
        {
            Json::Value directNode = read<Json::Value>(rollNode, "direct");

            cfg.direct.B44_lin = read<double>(directNode, "B44_lin");
            cfg.direct.B44_quad = read<double>(directNode, "B44_quad");
            cfg.direct.B44_cube = read<double>(directNode, "B44_cube");
        }
        else if (cfg.mode == RollDampingMode::FromDecayCsv)
        {
            Json::Value decayNode = read<Json::Value>(rollNode, "decay");

            cfg.decay.csvPath = read<std::string>(decayNode, "csvPath");

            if (decayNode.isMember("polyOrder"))
                cfg.decay.polyOrder = read<int>(decayNode, "polyOrder");
            else
                cfg.decay.polyOrder = 2;

            if (decayNode.isMember("minPeakGap"))
                cfg.decay.minPeakGap = read<int>(decayNode, "minPeakGap");
            else
                cfg.decay.minPeakGap = 3;

            if (decayNode.isMember("IeffOverride"))
                cfg.decay.IeffOverride = read<double>(decayNode, "IeffOverride");
            else
                cfg.decay.IeffOverride = -1.0;

            if (decayNode.isMember("angleInDeg"))
                cfg.decay.angleInDeg = read<bool>(decayNode, "angleInDeg");
            else
                cfg.decay.angleInDeg = false;

            if (decayNode.isMember("refine"))
                cfg.decay.refine = read<bool>(decayNode, "refine");
            if (decayNode.isMember("refineSkipFirstPeaks"))
                cfg.decay.refineSkipFirstPeaks =
                    read<int>(decayNode, "refineSkipFirstPeaks");
        }
        else
        {
            throw std::runtime_error("unsupported RollDamping mode");
        }

        return cfg;
    }

    FreeRollDecayConfig loadFreeRollDecayConfig(const Json::Value& node)
    {
        FreeRollDecayConfig cfg;

        if (!node.isMember("FreeRollDecay") || !node["FreeRollDecay"].isObject())
        {
            cfg.enabled = false;
            return cfg;
        }

        Json::Value frd = read<Json::Value>(node, "FreeRollDecay");

        cfg.enabled = read<bool>(frd, "enabled");
        cfg.phi0_deg = read<double>(frd, "phi0_deg");
        cfg.phidot0_deg_s = read<double>(frd, "phidot0_deg_s");
        cfg.dt = read<double>(frd, "dt");
        cfg.duration = read<double>(frd, "duration");

        return cfg;
    }

    SeakeepingConfig loadSeakeepingConfig(const Json::Value& node)
    {
        SeakeepingConfig seakeeping;

        //if (node.isMember("waves_seakeeping") && node["waves_seakeeping"].isArray())
        //{
        //    Json::Value waves_array = node["waves_seakeeping"];
        //    Json::Value wave_node;
        //    for (int i = 0; i < waves_array.size(); ++i)
        //    {
        //        wave_node = waves_array[i];
        //        auto wave_cfgs = loadWaves(wave_node);
        //        seakeeping.waves.insert(seakeeping.waves.end(), wave_cfgs.begin(), wave_cfgs.end());
        //    }
        //}
        if (node.isMember("waves_seakeeping") && node["waves_seakeeping"].isArray())
        {
            const Json::Value& waves_array = node["waves_seakeeping"];
            seakeeping.waves.reserve(waves_array.size());
            for (unsigned int i = 0; i < waves_array.size(); ++i)
            {
                auto wave_cfgs = loadWaves(waves_array[i]);
                seakeeping.waves.reserve(seakeeping.waves.size() + wave_cfgs.size());
                std::move(wave_cfgs.begin(), wave_cfgs.end(), std::back_inserter(seakeeping.waves));
            }
        }
        else
        {
            throw std::runtime_error("there is no waves in json file or waves is not array");
        }

      
        seakeeping.Fn = read<std::vector<double>>(node, "Fn");
        seakeeping.DOF = read<int>(node, "DOF");
        seakeeping.modes = read<std::vector<int>>(node, "modes");

        // PanelSettings
        Json::Value panel = read<Json::Value>(node, "PanelSettings");
        seakeeping.Panel.NEType = read<std::string>(panel, "NEType");
        seakeeping.Panel.NE = read<int>(panel, "NE");

        // TimeSettings
        Json::Value time = read<Json::Value>(node, "TimeSettings");
        seakeeping.Time.dt = read<double>(time, "dt");
        seakeeping.Time.TimeCircle = read<double>(time, "TimeCircle");
        seakeeping.Time.PreCircle = read<double>(time, "PreCircle");
        seakeeping.Time.GreenCircle = read<double>(time, "GreenCircle");

        Json::Value solver = read<Json::Value>(node, "Solver");
        seakeeping.Solver = read<std::string>(solver, "Method");

        if (solver.isMember("ResponseMethod"))
            seakeeping.ResponseMethod = read<std::string>(solver, "ResponseMethod");
        else
            seakeeping.ResponseMethod = "LegacyAB";

        if (solver.isMember("FKForceMethod"))
            seakeeping.FKMethod =
                parseFKForceMethod(read<std::string>(solver, "FKForceMethod"));
        else
            seakeeping.FKMethod = FKForceMethod::PrescribedAmp;


        // 新增：格林函数计算方式
        if (node.isMember("GreenFunction") && node["GreenFunction"].isObject())
        {
            Json::Value gfNode = read<Json::Value>(node, "GreenFunction");

            if (gfNode.isMember("Method"))
                seakeeping.GreenFunction.Method = read<std::string>(gfNode, "Method");

            if (gfNode.isMember("TablePath"))
                seakeeping.GreenFunction.TablePath = read<std::string>(gfNode, "TablePath");

            if (gfNode.isMember("RK4Step"))
                seakeeping.GreenFunction.RK4Step = read<double>(gfNode, "RK4Step");
        }
        else
        {
            seakeeping.GreenFunction.Method = "Exact";
            seakeeping.GreenFunction.TablePath = "green_table.bin";
            seakeeping.GreenFunction.RK4Step = 0.001;
        }


        // Optional temporary K(t) override (default off).
        if (node.isMember("KradOverride") && node["KradOverride"].isObject())
        {
            Json::Value ko = read<Json::Value>(node, "KradOverride");
            if (ko.isMember("enabled"))
                seakeeping.KradOverride.enabled = read<bool>(ko, "enabled");
            if (ko.isMember("file"))
                seakeeping.KradOverride.file = read<std::string>(ko, "file");
        }

        // Optional non-uniform chi-time grid for K(t) offline build (default off → 走老路径).
        // 示例：
        //   "ChiTimeGrid": {
        //     "enabled": true,
        //     "dt_fine": 0.05,
        //     "blocks": [
        //       { "stride": 1,  "count": 100 },
        //       { "stride": 4,  "count": 125 },
        //       { "stride": 20, "count": 170 }
        //     ],
        //     "memory_cutoff_lag": 800
        //   }
        if (node.isMember("ChiTimeGrid") && node["ChiTimeGrid"].isObject())
        {
            Json::Value cg = read<Json::Value>(node, "ChiTimeGrid");
            auto& dst = seakeeping.ChiTimeGrid;

            if (cg.isMember("enabled"))            dst.enabled = read<bool>(cg, "enabled");
            if (cg.isMember("dt_fine"))            dst.dt_fine = read<double>(cg, "dt_fine");
            if (cg.isMember("memory_cutoff_lag"))  dst.memory_cutoff_lag = read<int>(cg, "memory_cutoff_lag");

            if (cg.isMember("blocks") && cg["blocks"].isArray())
            {
                const Json::Value& arr = cg["blocks"];
                dst.blocks.clear();
                dst.blocks.reserve(arr.size());
                for (Json::ArrayIndex i = 0; i < arr.size(); ++i)
                {
                    if (!arr[i].isObject()) continue;
                    ChiTimeGridBlock b;
                    if (arr[i].isMember("stride")) b.stride = read<int>(arr[i], "stride");
                    if (arr[i].isMember("count"))  b.count  = read<int>(arr[i], "count");
                    if (b.stride >= 1 && b.count >= 1)
                        dst.blocks.push_back(b);
                }
            }

            // 完整性 check：开启了但参数不合法就降级回等距路径，给个 warn
            if (dst.enabled)
            {
                if (dst.dt_fine <= 0.0 || dst.blocks.empty())
                {
                    std::cerr << "[CaseLoader] ChiTimeGrid.enabled=true but config "
                                 "incomplete (dt_fine="
                              << dst.dt_fine << ", blocks=" << dst.blocks.size()
                              << "); falling back to uniform chi grid.\n";
                    dst.enabled = false;
                }
                else
                {
                    int totalNodes = 1;
                    int totalSteps = 0;
                    for (const auto& b : dst.blocks) {
                        totalNodes += b.count;
                        totalSteps += b.stride * b.count;
                    }
                    std::cout << "[CaseLoader] ChiTimeGrid enabled: dt_fine="
                              << dst.dt_fine
                              << ", blocks=" << dst.blocks.size()
                              << ", chi nodes=" << totalNodes
                              << ", tMax=" << totalSteps * dst.dt_fine << "s"
                              << ", memory_cutoff_lag=" << dst.memory_cutoff_lag
                              << " (= " << dst.memory_cutoff_lag * dst.dt_fine << "s)\n";
                }
            }
        }


        if (node.isMember("FKImpulseKernel") && node["FKImpulseKernel"].isObject())
        {
            auto& fk = seakeeping.FKImpulseKernel;
            Json::Value FKNode = read<Json::Value>(node, "FKImpulseKernel");

            if (FKNode.isMember("Method"))            fk.Method                  = read<std::string>(FKNode, "Method");
            if (FKNode.isMember("IncludeShipName"))   fk.IncludeShipName         = read<bool>       (FKNode, "IncludeShipName");

            // Physics-driven (dt, tMot) — all optional, defaults give sensible behaviour.
            if (FKNode.isMember("SamplesPerWavePeriod"))   fk.SamplesPerWavePeriod   = read<double>(FKNode, "SamplesPerWavePeriod");
            if (FKNode.isMember("DtFloorNd"))              fk.DtFloorNd              = read<double>(FKNode, "DtFloorNd");
            if (FKNode.isMember("DtCapNd"))                fk.DtCapNd                = read<double>(FKNode, "DtCapNd");
            if (FKNode.isMember("MemoryEncounterPeriods")) fk.MemoryEncounterPeriods = read<double>(FKNode, "MemoryEncounterPeriods");
            if (FKNode.isMember("MemoryFloorNd"))          fk.MemoryFloorNd          = read<double>(FKNode, "MemoryFloorNd");
            if (FKNode.isMember("MemoryCapNd"))            fk.MemoryCapNd            = read<double>(FKNode, "MemoryCapNd");
            if (FKNode.isMember("TailTolRel"))             fk.TailTolRel             = read<double>(FKNode, "TailTolRel");

            if (FKNode.isMember("RegionAdaptiveKernelGrid")) fk.RegionAdaptiveKernelGrid = read<bool>(FKNode, "RegionAdaptiveKernelGrid");
            if (FKNode.isMember("SamplesPerWavePeriodHead")) fk.SamplesPerWavePeriodHead = read<double>(FKNode, "SamplesPerWavePeriodHead");
            if (FKNode.isMember("DtCapNdHead"))             fk.DtCapNdHead             = read<double>(FKNode, "DtCapNdHead");
            if (FKNode.isMember("MemoryEncounterPeriodsHead")) fk.MemoryEncounterPeriodsHead = read<double>(FKNode, "MemoryEncounterPeriodsHead");
            if (FKNode.isMember("MemoryFloorNdHead"))       fk.MemoryFloorNdHead       = read<double>(FKNode, "MemoryFloorNdHead");
            if (FKNode.isMember("MemoryCapNdHead"))         fk.MemoryCapNdHead         = read<double>(FKNode, "MemoryCapNdHead");
            if (FKNode.isMember("SamplesPerWavePeriodFollowing")) fk.SamplesPerWavePeriodFollowing = read<double>(FKNode, "SamplesPerWavePeriodFollowing");
            if (FKNode.isMember("DtCapNdFollowing"))        fk.DtCapNdFollowing        = read<double>(FKNode, "DtCapNdFollowing");
            if (FKNode.isMember("MemoryEncounterPeriodsFollowing")) fk.MemoryEncounterPeriodsFollowing = read<double>(FKNode, "MemoryEncounterPeriodsFollowing");
            if (FKNode.isMember("MemoryFloorNdFollowing"))  fk.MemoryFloorNdFollowing  = read<double>(FKNode, "MemoryFloorNdFollowing");
            if (FKNode.isMember("MemoryCapNdFollowing"))    fk.MemoryCapNdFollowing    = read<double>(FKNode, "MemoryCapNdFollowing");
            if (FKNode.isMember("MemoryCapNdF1"))           fk.MemoryCapNdF1           = read<double>(FKNode, "MemoryCapNdF1");
            if (FKNode.isMember("MemoryEncounterPeriodsF2")) fk.MemoryEncounterPeriodsF2 = read<double>(FKNode, "MemoryEncounterPeriodsF2");
            if (FKNode.isMember("MemoryCapNdF2"))           fk.MemoryCapNdF2           = read<double>(FKNode, "MemoryCapNdF2");
            if (FKNode.isMember("MemoryCapNdF3"))           fk.MemoryCapNdF3           = read<double>(FKNode, "MemoryCapNdF3");
            if (FKNode.isMember("FollowingMemoryMinIncidentPeriods")) fk.FollowingMemoryMinIncidentPeriods = read<double>(FKNode, "FollowingMemoryMinIncidentPeriods");
            if (FKNode.isMember("MaxFKImpulseRows")) fk.MaxFKImpulseRows = read<int>(FKNode, "MaxFKImpulseRows");
        }


        seakeeping.RollDamping = loadRollDampingConfig(node);
        seakeeping.FreeRollDecay = loadFreeRollDecayConfig(node);
        
        if (seakeeping.RollDamping.enabled)
        {
            std::cout << "RollDamping enabled\n";
            if (seakeeping.RollDamping.mode == RollDampingMode::DirectCoefficients)
            {
                std::cout << "RollDamping mode: DirectCoefficients\n";
                std::cout << "B44_lin:  " << seakeeping.RollDamping.direct.B44_lin << "\n";
                std::cout << "B44_quad: " << seakeeping.RollDamping.direct.B44_quad << "\n";
                std::cout << "B44_cube: " << seakeeping.RollDamping.direct.B44_cube << "\n";
            }
            else
            {
                std::cout << "RollDamping mode: FromDecayCsv\n";
                std::cout << "csvPath:      " << seakeeping.RollDamping.decay.csvPath << "\n";
                std::cout << "polyOrder:    " << seakeeping.RollDamping.decay.polyOrder << "\n";
                std::cout << "minPeakGap:   " << seakeeping.RollDamping.decay.minPeakGap << "\n";
                std::cout << "IeffOverride: " << seakeeping.RollDamping.decay.IeffOverride << "\n";
                std::cout << "angleInDeg:   " << seakeeping.RollDamping.decay.angleInDeg << "\n";
            }
        }

        return seakeeping;
    }

    MmgConfig loadMmgConfig(const Json::Value& node)
    {
        MmgConfig mmg;

        if (node.isMember("waves_manoeuvring") && node["waves_manoeuvring"].isArray())
        {
            const Json::Value& waves_array = node["waves_manoeuvring"];
            mmg.waves.reserve(waves_array.size());
            for (unsigned int i = 0; i < waves_array.size(); ++i)
            {
                const Json::Value& wave_node = waves_array[i];
                auto wave_cfgs = loadWaves(wave_node);
                mmg.waves.reserve(mmg.waves.size() + wave_cfgs.size());
                std::move(wave_cfgs.begin(), wave_cfgs.end(), std::back_inserter(mmg.waves));
            }
        }

        else
        {
            // 静水允许没有 waves_manoeuvring
            mmg.waves.clear();
            //throw std::runtime_error("there is no waves in json file or waves is not array");
        }

        // MmgTimeConfig
        Json::Value time = read<Json::Value>(node, "MmgTimeConfig");
        mmg.Time.dt = read<double>(time, "dt");
        mmg.Time.SeakeepingStep = read<int>(time, "SeakeepingStep");
        mmg.Time.v_interavl = read<double>(time, "v_interavl");
        mmg.Time.angle_interval = read<double>(time, "angle_interval");

        // hullCoeff
        Json::Value hullCoeff = read<Json::Value>(node, "hullCoeff");
        mmg.Hull.Xvv = read<double>(hullCoeff, "Xvv");
        mmg.Hull.Xvr = read<double>(hullCoeff, "Xvr");
        mmg.Hull.Xrr = read<double>(hullCoeff, "Xrr");
        mmg.Hull.Xvvvv = read<double>(hullCoeff, "Xvvvv");
        mmg.Hull.Xuu = read<double>(hullCoeff, "Xuu");

        mmg.Hull.Yv = read<double>(hullCoeff, "Yv");
        mmg.Hull.Yr = read<double>(hullCoeff, "Yr");
        mmg.Hull.Yvvv = read<double>(hullCoeff, "Yvvv");
        mmg.Hull.Yvvr = read<double>(hullCoeff, "Yvvr");
        mmg.Hull.Yvrr = read<double>(hullCoeff, "Yvrr");
        mmg.Hull.Yrrr = read<double>(hullCoeff, "Yrrr");

        mmg.Hull.Nv = read<double>(hullCoeff, "Nv");
        mmg.Hull.Nr = read<double>(hullCoeff, "Nr");
        mmg.Hull.Nvvv = read<double>(hullCoeff, "Nvvv");
        mmg.Hull.Nvvr = read<double>(hullCoeff, "Nvvr");
        mmg.Hull.Nvrr = read<double>(hullCoeff, "Nvrr");
        mmg.Hull.Nrrr = read<double>(hullCoeff, "Nrrr");

        // propellerCoeff
        Json::Value propellerCoeff = read<Json::Value>(node, "propellerCoeff");
        mmg.Propeller.tp = read<double>(propellerCoeff, "tp");
        mmg.Propeller.np = read<double>(propellerCoeff, "np");
        mmg.Propeller.Dp = read<double>(propellerCoeff, "Dp");
        mmg.Propeller.k0 = read<double>(propellerCoeff, "k0");
        mmg.Propeller.k1 = read<double>(propellerCoeff, "k1");
        mmg.Propeller.k2 = read<double>(propellerCoeff, "k2");
        mmg.Propeller.C1 = read<double>(propellerCoeff, "C1");
        mmg.Propeller.C2 = read<double>(propellerCoeff, "C2");
        mmg.Propeller.xp = read<double>(propellerCoeff, "xp");
        mmg.Propeller.wP0 = read<double>(propellerCoeff, "wP0");

        // rudderCoeff
        Json::Value rudderCoeff = read<Json::Value>(node, "rudderCoeff");
        mmg.Rudder.omigaR = read<double>(rudderCoeff, "omigaR");
        mmg.Rudder.HR = read<double>(rudderCoeff, "HR");
        mmg.Rudder.kappa = read<double>(rudderCoeff, "kappa");
        mmg.Rudder.lR = read<double>(rudderCoeff, "lR");
        mmg.Rudder.AR = read<double>(rudderCoeff, "AR");
        mmg.Rudder.f_alpha = read<double>(rudderCoeff, "f_alpha");
        mmg.Rudder.tR = read<double>(rudderCoeff, "tR");
        mmg.Rudder.aH = read<double>(rudderCoeff, "aH");
        mmg.Rudder.xR = read<double>(rudderCoeff, "xR");
        mmg.Rudder.xH = read<double>(rudderCoeff, "xH");
        mmg.Rudder.epsilonP = read<double>(rudderCoeff, "epsilonP");
        mmg.Rudder.kP = read<double>(rudderCoeff, "kP");
        mmg.Rudder.gammaR_plus = read<double>(rudderCoeff, "gammaR_plus");
        mmg.Rudder.gammaR_minus = read<double>(rudderCoeff, "gammaR_minus");

        // addedCoeff
        Json::Value addedCoeff = read<Json::Value>(node, "addedCoeff");
        mmg.Added.mx = read<double>(addedCoeff, "mx");
        mmg.Added.my = read<double>(addedCoeff, "my");
        mmg.Added.Jz = read<double>(addedCoeff, "Jz");

        mmg.Fn = read<std::vector<double>>(node, "Fn");

        Json::Value TurningCfg = read<Json::Value>(node, "TurningCfg");
        mmg.Turning.TurningCaseCicle = read<double>(TurningCfg, "TurningCaseCicle");
        mmg.Turning.rudder_angle = read<double>(TurningCfg, "rudder_angle") * PI / 180.0;
        mmg.Turning.rudder_rate = read<double>(TurningCfg, "rudder_rate") * PI / 180.0;
        // 可选停止判据（缺省 "circle" + 0，等同原按圈数停的行为）
        if (TurningCfg.isMember("stop_by") && TurningCfg["stop_by"].isString())
            mmg.Turning.stopBy = TurningCfg["stop_by"].asString();
        if (TurningCfg.isMember("t_star_max") && TurningCfg["t_star_max"].isNumeric())
            mmg.Turning.t_star_max = TurningCfg["t_star_max"].asDouble();

        Json::Value ZigzagCfg = read<Json::Value>(node, "ZigzagCfg");
        mmg.Zigzag.ZigzagCaseCicle = read<double>(ZigzagCfg, "ZigzagCaseCicle");
        mmg.Zigzag.rudder_angle = read<double>(ZigzagCfg, "rudder_angle") * PI / 180.0;
        mmg.Zigzag.rudder_rate = read<double>(ZigzagCfg, "rudder_rate") * PI / 180.0;

        mmg.hullCoeff_defined = read<bool>(node, "use_hullCoeff_defined");
        mmg.propellerCoeff_defined = read<bool>(node, "use_propellerCoeff_defined");
        mmg.rudderCoeff_defined = read<bool>(node, "use_rudderCoeff_defined");
        mmg.addedCoeff_defined = read<bool>(node, "addedCoeff_defined");
        mmg.turn_on_turningcase = read<bool>(node, "turn_on_turningcase");
        mmg.turn_on_zigzagcase = read<bool>(node, "turn_on_zigzagcase");

        if (node.isMember("CoupledMmg3DOF") && node["CoupledMmg3DOF"].isObject())
        {
            const Json::Value& c = node["CoupledMmg3DOF"];
            mmg.coupledMmg3DOF.enabled = read<bool>(c, "enabled");
            if (mmg.coupledMmg3DOF.enabled)
            {
                if (c.isMember("reuseManoeuvringNondimHullProp"))
                    mmg.coupledMmg3DOF.reuseManoeuvringNondimHullProp =
                        read<bool>(c, "reuseManoeuvringNondimHullProp");
                mmg.coupledMmg3DOF.rho = read<double>(c, "rho");
                mmg.coupledMmg3DOF.g = read<double>(c, "g");
                mmg.coupledMmg3DOF.Lpp = read<double>(c, "Lpp");
                mmg.coupledMmg3DOF.Lwl = read<double>(c, "Lwl");
                mmg.coupledMmg3DOF.B = read<double>(c, "B");
                mmg.coupledMmg3DOF.d = read<double>(c, "d");
                mmg.coupledMmg3DOF.Cb = read<double>(c, "Cb");
                mmg.coupledMmg3DOF.massKg = read<double>(c, "massKg");
                mmg.coupledMmg3DOF.xG_m = read<double>(c, "xG_m");
                mmg.coupledMmg3DOF.kzz_nd = read<double>(c, "kzz_nd");
                mmg.coupledMmg3DOF.Dp = read<double>(c, "Dp");
                mmg.coupledMmg3DOF.HR = read<double>(c, "HR");
                mmg.coupledMmg3DOF.BR = read<double>(c, "BR");
                mmg.coupledMmg3DOF.AR = read<double>(c, "AR");
                mmg.coupledMmg3DOF.xP_nd = read<double>(c, "xP_nd");
                mmg.coupledMmg3DOF.xR_nd = read<double>(c, "xR_nd");
                mmg.coupledMmg3DOF.xH_nd = read<double>(c, "xH_nd");
                mmg.coupledMmg3DOF.lR_nd = read<double>(c, "lR_nd");
                mmg.coupledMmg3DOF.epsilon = read<double>(c, "epsilon");
                mmg.coupledMmg3DOF.kappa = read<double>(c, "kappa");
                mmg.coupledMmg3DOF.tR = read<double>(c, "tR");
                mmg.coupledMmg3DOF.aH = read<double>(c, "aH");
                mmg.coupledMmg3DOF.gammaR_plus = read<double>(c, "gammaR_plus");
                mmg.coupledMmg3DOF.gammaR_minus = read<double>(c, "gammaR_minus");
                mmg.coupledMmg3DOF.deltaRateDegPerS = read<double>(c, "deltaRateDegPerS");
                mmg.coupledMmg3DOF.maxDeltaDeg = read<double>(c, "maxDeltaDeg");
                mmg.coupledMmg3DOF.Fn = read<double>(c, "Fn");
                mmg.coupledMmg3DOF.npRps = read<double>(c, "npRps");
                mmg.coupledMmg3DOF.mx_nd = read<double>(c, "mx_nd");
                mmg.coupledMmg3DOF.my_nd = read<double>(c, "my_nd");
                mmg.coupledMmg3DOF.Jz_nd = read<double>(c, "Jz_nd");

                // Optional 4-DOF roll block (all fields optional -> defaults keep 3-DOF).
                if (c.isMember("maneuverDOF"))   mmg.coupledMmg3DOF.maneuverDOF   = read<int>(c, "maneuverDOF");
                if (c.isMember("rollJxxAdd"))    mmg.coupledMmg3DOF.rollJxxAdd    = read<double>(c, "rollJxxAdd");
                if (c.isMember("rollB44_lin"))   mmg.coupledMmg3DOF.rollB44_lin   = read<double>(c, "rollB44_lin");
                if (c.isMember("rollB44_quad"))  mmg.coupledMmg3DOF.rollB44_quad  = read<double>(c, "rollB44_quad");
                if (c.isMember("zH_over_d"))     mmg.coupledMmg3DOF.zH_over_d     = read<double>(c, "zH_over_d");
                if (c.isMember("zR_over_d"))     mmg.coupledMmg3DOF.zR_over_d     = read<double>(c, "zR_over_d");
                if (c.isMember("zG_over_d"))     mmg.coupledMmg3DOF.zG_over_d     = read<double>(c, "zG_over_d");
            }
        }

        if (node.isMember("lowFreqManeuverModel") && !node["lowFreqManeuverModel"].isNull())
            mmg.lowFreqManeuverModel = node["lowFreqManeuverModel"].asString();

        if (node.isMember("CoupledSjtuMmg3DOF") && node["CoupledSjtuMmg3DOF"].isObject())
        {
            const Json::Value& sj = node["CoupledSjtuMmg3DOF"];
            mmg.sjtuMmg.enabled = read<bool>(sj, "enabled");
            mmg.sjtuMmg.rho = read<double>(sj, "rho");
            mmg.sjtuMmg.g = read<double>(sj, "g");
            mmg.sjtuMmg.Lpp = read<double>(sj, "Lpp");
            mmg.sjtuMmg.Lwl = read<double>(sj, "Lwl");
            mmg.sjtuMmg.B = read<double>(sj, "B");
            mmg.sjtuMmg.d = read<double>(sj, "d");
            mmg.sjtuMmg.massKg = read<double>(sj, "massKg");
            mmg.sjtuMmg.kzz_nd = read<double>(sj, "kzz_nd");
            mmg.sjtuMmg.S_wetted = read<double>(sj, "S_wetted");
            mmg.sjtuMmg.Dp = read<double>(sj, "Dp");
            mmg.sjtuMmg.HR = read<double>(sj, "HR");
            mmg.sjtuMmg.BR = read<double>(sj, "BR");
            mmg.sjtuMmg.AR = read<double>(sj, "AR");
            mmg.sjtuMmg.xP_nd = read<double>(sj, "xP_nd");
            mmg.sjtuMmg.xR_nd = read<double>(sj, "xR_nd");
            mmg.sjtuMmg.xH_nd = read<double>(sj, "xH_nd");
            mmg.sjtuMmg.lR_nd = read<double>(sj, "lR_nd");
            mmg.sjtuMmg.Xvv = read<double>(sj, "Xvv");
            mmg.sjtuMmg.Xvr = read<double>(sj, "Xvr");
            mmg.sjtuMmg.Xrr = read<double>(sj, "Xrr");
            mmg.sjtuMmg.Yv = read<double>(sj, "Yv");
            mmg.sjtuMmg.Yr = read<double>(sj, "Yr");
            mmg.sjtuMmg.Yvvv = read<double>(sj, "Yvvv");
            mmg.sjtuMmg.Yvvr = read<double>(sj, "Yvvr");
            mmg.sjtuMmg.Yvrr = read<double>(sj, "Yvrr");
            mmg.sjtuMmg.Yrrr = read<double>(sj, "Yrrr");
            mmg.sjtuMmg.Nv = read<double>(sj, "Nv");
            mmg.sjtuMmg.Nr = read<double>(sj, "Nr");
            mmg.sjtuMmg.Nvvv = read<double>(sj, "Nvvv");
            mmg.sjtuMmg.Nvvr = read<double>(sj, "Nvvr");
            mmg.sjtuMmg.Nvrr = read<double>(sj, "Nvrr");
            mmg.sjtuMmg.Nrrr = read<double>(sj, "Nrrr");
            mmg.sjtuMmg.mx_nd = read<double>(sj, "mx_nd");
            mmg.sjtuMmg.my_nd = read<double>(sj, "my_nd");
            mmg.sjtuMmg.Jz_nd = read<double>(sj, "Jz_nd");
            mmg.sjtuMmg.Cr = read<double>(sj, "Cr");
            mmg.sjtuMmg.kf = read<double>(sj, "kf");
            if (sj.isMember("nuWater"))
                mmg.sjtuMmg.nuWater = read<double>(sj, "nuWater");
            mmg.sjtuMmg.tp = read<double>(sj, "tp");
            mmg.sjtuMmg.wP0 = read<double>(sj, "wP0");
            mmg.sjtuMmg.j0 = read<double>(sj, "j0");
            mmg.sjtuMmg.j1 = read<double>(sj, "j1");
            mmg.sjtuMmg.j2 = read<double>(sj, "j2");
            if (sj.isMember("wakeExpBetaP"))
                mmg.sjtuMmg.wakeExpBetaP = read<double>(sj, "wakeExpBetaP");
            mmg.sjtuMmg.epsilon = read<double>(sj, "epsilon");
            mmg.sjtuMmg.kappa = read<double>(sj, "kappa");
            mmg.sjtuMmg.tR = read<double>(sj, "tR");
            mmg.sjtuMmg.aH = read<double>(sj, "aH");
            mmg.sjtuMmg.gammaR_plus = read<double>(sj, "gammaR_plus");
            mmg.sjtuMmg.gammaR_minus = read<double>(sj, "gammaR_minus");
            if (sj.isMember("tauPropInflow"))
                mmg.sjtuMmg.tauPropInflow = read<double>(sj, "tauPropInflow");
            mmg.sjtuMmg.C_Rv = read<double>(sj, "C_Rv");
            mmg.sjtuMmg.C_Rrv = read<double>(sj, "C_Rrv");
            mmg.sjtuMmg.C_Rrrv = read<double>(sj, "C_Rrrv");
            mmg.sjtuMmg.deltaRateDegPerS = read<double>(sj, "deltaRateDegPerS");
            mmg.sjtuMmg.maxDeltaDeg = read<double>(sj, "maxDeltaDeg");
            mmg.sjtuMmg.Fn = read<double>(sj, "Fn");
            mmg.sjtuMmg.npRps = read<double>(sj, "npRps");
        }

        return mmg;
    }


    CoupledConfig loadCouplingConfig(const Json::Value& root)
    {
        CoupledConfig cfg{};

        const Json::Value& c = root;

        if (c.isMember("Time") && c["Time"].isObject())
        {
            const Json::Value& t = c["Time"];
            cfg.time.dtSlow = read<double>(t, "dtSlow");
            cfg.time.dtFast = read<double>(t, "dtFast");
            cfg.time.totalTime = read<double>(t, "totalTime");
            cfg.time.predictorCorrectorIters =
                read<int>(t, "predictorCorrectorIters");
            if (t.isMember("waveRampTimeS"))
                cfg.time.waveRampTimeS = read<double>(t, "waveRampTimeS");
        }

        if (c.isMember("Refresh") && c["Refresh"].isObject())
        {
            const Json::Value& r = c["Refresh"];
            cfg.refresh.betaTolDeg = read<double>(r, "betaTolDeg");
            cfg.refresh.omegaTolRatio = read<double>(r, "omegaTolRatio");
            cfg.refresh.speedTolRatio = read<double>(r, "speedTolRatio");
            cfg.refresh.enableExcitationInterpolation =
                read<bool>(r, "enableExcitationInterpolation");
            cfg.refresh.enableRadiationInterpolation =
                read<bool>(r, "enableRadiationInterpolation");
            if (r.isMember("adaptiveRefresh"))
                cfg.refresh.adaptiveRefresh = read<bool>(r, "adaptiveRefresh");
            if (r.isMember("adaptiveRefreshFactor"))
                cfg.refresh.adaptiveRefreshFactor =
                    read<double>(r, "adaptiveRefreshFactor");
        }

        if (c.isMember("Output") && c["Output"].isObject())
        {
            const Json::Value& o = c["Output"];
            cfg.output.saveSlowStates = read<bool>(o, "saveSlowStates");
            cfg.output.saveWindowLoads = read<bool>(o, "saveWindowLoads");
            cfg.output.saveFastSummary = read<bool>(o, "saveFastSummary");
            if (o.isMember("referenceModelLppM"))
                cfg.output.referenceModelLppM = read<double>(o, "referenceModelLppM");
        }

        if (c.isMember("physics"))
        {
            const Json::Value& phy = c["physics"];

            if (phy.isMember("enableFastSeakeeping"))
                cfg.physics.enableFastSeakeeping = phy["enableFastSeakeeping"].asBool();

            if (phy.isMember("enableSecondOrderLoads"))
                cfg.physics.enableSecondOrderLoads = phy["enableSecondOrderLoads"].asBool();

            if (phy.isMember("waveForceMode"))
                cfg.physics.waveForceMode = phy["waveForceMode"].asString();

            if (phy.isMember("enableRadiation"))
                cfg.physics.enableRadiation = phy["enableRadiation"].asBool();
        }

        return cfg;
    }

    struct Vec3 { double x{}, y{}, z{}; };

    static inline bool startsWith(const std::string& s, const char* p) {
        return s.rfind(p, 0) == 0;
    }

    static inline std::string trim(std::string s) {
        auto notSpace = [](unsigned char c) { return !std::isspace(c); };
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
        s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
        return s;
    }

    static inline void sanitizeNum(std::string& s) {
        s = trim(s);
        while (!s.empty() && (s.back() == '+' || s.back() == '*')) s.pop_back();
        for (char& c : s) if (c == 'D') c = 'E';
    }

    static inline std::string field16(const std::string& line, int k) {
        // Nastran large-field  8(    )+4*16  ֶ 
        const int start = 8 + k * 16;
        if (start >= (int)line.size()) return {};
        std::string s = line.substr(start, 16);
        sanitizeNum(s);
        return s;
    }

    static void readNodes(std::ifstream& in,
        std::unordered_map<int, Vec3>& nodes,
        double scale)
    {
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '$') continue;

            if (startsWith(line, "GRID*")) {
                std::string line2;
                if (!std::getline(in, line2)) break;

                int nid = std::stoi(field16(line, 0));
                double x = std::stod(field16(line, 2)) * scale;
                double y = std::stod(field16(line, 3)) * scale;
                double z = std::stod(field16(line2, 0)) * scale;

                nodes[nid] = Vec3{ x, y, z };
            }
            else if (startsWith(line, "GRID")) {
                //    ɸ ʽ  GRID id cp x y z ...
                std::istringstream iss(line);
                std::string card;
                int nid = 0, cp = 0;
                double x = 0, y = 0, z = 0;
                iss >> card >> nid >> cp >> x >> y >> z;
                if (iss && card == "GRID") {
                    nodes[nid] = Vec3{ x * scale, y * scale, z * scale };
                }
            }
        }
    }



    void validateElement(ElementMatrix& Elem, const double EPS)
    {
        // 1) ˮ ¼ 飺   е  z      <= EPS  z=0  ˮ    
        if (Elem(0, 2) > EPS || Elem(1, 2) > EPS || Elem(2, 2) > EPS || Elem(3, 2) > EPS)
            throw std::runtime_error("Element are not under water!\n");

        // 2) ȡ   ĵ 
        Eigen::Vector3d c = 0.25 * (Elem.row(0) + Elem.row(1) + Elem.row(2) + Elem.row(3)).transpose();

        Eigen::Vector3d p[4];
        for (int i = 0; i < 4; ++i) p[i] = Elem.row(i).transpose();

        // 3)     һ   ɿ   ƽ 淨      ͶӰ   򣬲 Ҫ      ȷ  
        Eigen::Vector3d n = (p[1] - p[0]).cross(p[2] - p[0]);
        //if (n.squaredNorm() < EPS * EPS) n = (p[2] - p[0]).cross(p[3] - p[0]);
        //std::cout << std::setprecision(10) << "n.squaredNorm():\t" << n.squaredNorm() << std::endl;
        //if (n.squaredNorm() < EPS * EPS) n = (p[3] - p[0]).cross(p[1] - p[0]);
        //std::cout << std::setprecision(10) << "n.squaredNorm():\t" << n.squaredNorm() << std::endl;
        //if (n.squaredNorm() < EPS * EPS)
        //    throw std::runtime_error("Degenerate element (cannot build plane).\n");
        n.normalize();

        // 4) ƽ   ڻ    (u,v)
        auto proj = [&](const Eigen::Vector3d& v) { return v - n * (v.dot(n)); };

        Eigen::Vector3d u = proj(p[0] - c);
        if (u.squaredNorm() < EPS * EPS) u = proj(p[1] - c);
        if (u.squaredNorm() < EPS * EPS) u = proj(p[2] - c);
        if (u.squaredNorm() < EPS * EPS) u = proj(p[3] - c);
        if (u.squaredNorm() < EPS * EPS)
            throw std::runtime_error("Degenerate element (cannot build in-plane axis).\n");
        u.normalize();

        Eigen::Vector3d v = n.cross(u); // in-plane

        // 5)          ->  õ    Խ  Ļ 
        std::array<int, 4>      idx{ 0, 1, 2, 3 };
        std::array<double, 4>   ang{};

        for (int k = 0; k < 4; ++k)
        {
            Eigen::Vector3d d = p[k] - c;
            ang[k] = std::atan2(d.dot(v), d.dot(u));
        }

        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return ang[a] < ang[b]; });

        //   ת        ԭ   ĵ 0    Ϊ  һ  
        int pos0 = 0;
        for (int i = 0; i < 4; ++i) if (idx[i] == 0) { pos0 = i; break; }
        std::array<int, 4> id2{};
        for (int i = 0; i < 4; ++i) id2[i] = idx[(pos0 + i) % 4];

        ElementMatrix sorted;
        for (int i = 0; i < 4; ++i) sorted.row(i) = Elem.row(id2[i]);
        Elem = sorted;

        //               +            12/34  Խ  ߽   
        auto to2D = [&](const Eigen::Vector3d& P) {
            Eigen::Vector3d d = P - c;
            return Eigen::Vector2d(d.dot(u), d.dot(v));
        };

        Eigen::Vector2d q2d[4];
        for (int i = 0; i < 4; ++i) q2d[i] = to2D(Elem.row(i).transpose());

        auto orient = [&](const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c2) {
            return (b.x() - a.x()) * (c2.y() - a.y()) - (b.y() - a.y()) * (c2.x() - a.x());
        };

        auto properIntersect = [&](const Eigen::Vector2d& a, const Eigen::Vector2d& b,
            const Eigen::Vector2d& c2, const Eigen::Vector2d& d2) {
            double o1 = orient(a, b, c2);
            double o2 = orient(a, b, d2);
            double o3 = orient(c2, d2, a);
            double o4 = orient(c2, d2, b);

            // ֻ   ġ  ϸ  ཻ            / ˵㣩     EPS     
            if (std::abs(o1) < EPS || std::abs(o2) < EPS || std::abs(o3) < EPS || std::abs(o4) < EPS) return false;

            //o1 * o2 < 0.0  c,d  ab   ࣬ o3 * o4 < 0.0  a,b  cd    
            return (o1 * o2 < 0.0) && (o3 * o4 < 0.0);
        };

        //    (0-1)  (2-3) ཻ =>      1    2
        if (properIntersect(q2d[0], q2d[1], q2d[2], q2d[3]))
        {
            Eigen::RowVector3d tmp = Elem.row(1);
            Elem.row(1) = Elem.row(2);
            Elem.row(2) = tmp;
            //     2D  
            for (int i = 0; i < 4; ++i) q2d[i] = to2D(Elem.row(i).transpose());
        }

        //    (1-2)  (3-0) ཻ =>      2    3
        if (properIntersect(q2d[1], q2d[2], q2d[3], q2d[0]))
        {
            Eigen::RowVector3d tmp = Elem.row(2);
            Elem.row(2) = Elem.row(3);
            Elem.row(3) = tmp;
            //     2D  
            for (int i = 0; i < 4; ++i) q2d[i] = to2D(Elem.row(i).transpose());
        }

        // 6)   󣺱 ֤     z > 0        ͷ ת    
        Eigen::Vector3d P0 = Elem.row(0).transpose();
        Eigen::Vector3d P1 = Elem.row(1).transpose();
        Eigen::Vector3d P2 = Elem.row(2).transpose();
        Eigen::Vector3d P3 = Elem.row(3).transpose();

        Eigen::Vector3d nn = (P1 - P0).cross(P2 - P0) + (P2 - P0).cross(P3 - P0);
        //if (nn.squaredNorm() < EPS * EPS)
        //    throw std::runtime_error("Degenerate element (normal too small).\n");

        if (nn.z() < 0.0)
        {
            //   ת         1    3  ȼ  ڰ ˳  ת  0,1,2,3 -> 0,3,2,1  
            Eigen::RowVector3d tmp = Elem.row(1);
            Elem.row(1) = Elem.row(3);
            Elem.row(3) = tmp;
        }
    }

}


CaseConfig CaseLoader::loadcase(const std::string& json_file)
{
    std::ifstream json_ifstream(json_file);
    if (!json_ifstream.is_open())
    {
        throw std::runtime_error("can't open json file: " + json_file);
    }

    Json::Value  root;
    Json::Reader reader;

    if (!reader.parse(json_ifstream, root))
        throw std::runtime_error(reader.getFormattedErrorMessages());

    json_ifstream.close();

    CaseConfig cfg;

    cfg.Ship = loadShipConfig(root["Ship"]);
    cfg.enable_seakeeping = read<bool>(root, "enable_seakeeping");
    cfg.enable_maneuvering = read<bool>(root, "enable_maneuvering");
    cfg.enable_coupling = read<bool>(root, "enable_coupling");

    //if (root.isMember("waves") && root["waves"].isArray())
    //{
    //    Json::Value waves_array = root["waves"];
    //    Json::Value wave_node;
    //    for (int i = 0; i < waves_array.size(); ++i)
    //    {
    //        wave_node = waves_array[i];
    //        auto wave_cfgs = loadWaves(wave_node);
    //        cfg.waves.insert(cfg.waves.end(), wave_cfgs.begin(), wave_cfgs.end());
    //    }
    //}
    //else
    //{
    //    throw std::runtime_error("there is no waves in json file or waves is not array");
    //}

    if (cfg.enable_seakeeping)
    {
        if (root.isMember("Seakeeping") && root["Seakeeping"].isObject())
        {
            cfg.Seakeeping = loadSeakeepingConfig(root["Seakeeping"]);
        }
        else
        {
            throw std::runtime_error("use seakeeping but there is no seakeeping in json file!\n");
        }
    }

    if (cfg.enable_maneuvering)
    {
        if (root.isMember("Manoeuvring") && root["Manoeuvring"].isObject())
        {
            cfg.Mmg = loadMmgConfig(root["Manoeuvring"]);
        }
        else
        {
            throw std::runtime_error("use maneuvering but there is no Manoeuvring in json file!\n");
        }
    }

    if (cfg.enable_coupling)
    {
        if (root.isMember("Seakeeping") && root["Seakeeping"].isObject())
        {
            cfg.Seakeeping = loadSeakeepingConfig(root["Seakeeping"]);
        }
        else
        {
            throw std::runtime_error("use coupling but there is no seakeeping in json file!\n");
        }

        if (root.isMember("Manoeuvring") && root["Manoeuvring"].isObject())
        {
            cfg.Mmg = loadMmgConfig(root["Manoeuvring"]);
        }
        else
        {
            throw std::runtime_error("use coupling but there is no Manoeuvring in json file!\n");
        }


        if (root.isMember("Coupling") && root["Coupling"].isObject())
        {
            cfg.Coupled = loadCouplingConfig(root["Coupling"]);
        }
        else
        {
            throw std::runtime_error("use coupling but there is no coupling in json file!\n");
        }
    }
    return cfg;
}


void CaseLoader::UGtoElement(const std::string& datFile, const std::string& elementFile, const double scale)
{
    std::ifstream in(datFile);
    if (!in.is_open())
        throw std::runtime_error("UGtoElement: cannot open dat file: " + datFile);

    std::unordered_map<int, Vec3> nodes;
    nodes.reserve(100000);

    //        GRID/GRID*
    readNodes(in, nodes, scale);
    if (nodes.empty())
        throw std::runtime_error("UGtoElement: no GRID/GRID* found in: " + datFile);

    //    CQUAD4   д .element
    in.clear();
    in.seekg(0, std::ios::beg);

    std::ofstream out(elementFile);
    if (!out.is_open())
        throw std::runtime_error("UGtoElement: cannot open output: " + elementFile);

    out.setf(std::ios::fixed);
    out << std::setprecision(10);

    auto getNode = [&](int nid) -> const Vec3& {
        auto it = nodes.find(nid);
        if (it == nodes.end())
            throw std::runtime_error("UGtoElement: CQUAD4 references missing GRID id: " + std::to_string(nid));
        return it->second;
    };

    std::string line;
    std::size_t quadCount = 0;

    while (std::getline(in, line)) {
        if (!startsWith(line, "CQUAD4")) continue;

        std::istringstream iss(line);
        std::string card;
        int eid = 0, pid = 0, n1 = 0, n2 = 0, n3 = 0, n4 = 0;
        iss >> card >> eid >> pid >> n1 >> n2 >> n3 >> n4;
        if (!iss) continue;

        const Vec3& p1 = getNode(n1);
        const Vec3& p2 = getNode(n2);
        const Vec3& p3 = getNode(n3);
        const Vec3& p4 = getNode(n4);

        // ˳   ϸ  n1 n2 n3 n4    
        out << eid << "\n";
        out << p1.x << " " << p1.y << " " << p1.z << "\n";
        out << p2.x << " " << p2.y << " " << p2.z << "\n";
        out << p3.x << " " << p3.y << " " << p3.z << "\n";
        out << p4.x << " " << p4.y << " " << p4.z << "\n";

        ++quadCount;
    }

    if (quadCount == 0)
        throw std::runtime_error("UGToElement: no CQUAD4 found in: " + datFile);
}

std::unique_ptr<std::vector<ElementMatrix>> CaseLoader::loadelement
(std::string file, std::string NEType, int NE)
{
    std::ifstream ifile(file);
    if (!ifile.is_open())
        throw std::runtime_error("can't open element file!\n");

    auto Element = std::make_unique<std::vector<ElementMatrix>>(NE);

    int temp;
    int row;
    ElementMatrix       Elem;
    Eigen::RowVector3d  El_middle;

    const double EPS = 1.0e-3;

    if (NEType == "halfship")
    {
        if (NE % 2 != 0)
        {
            throw std::runtime_error("halfship but number of element is singler number\n");
        }
        int NE_half = NE / 2;

        for (int i = 0; i < NE_half; ++i)
        {
            ifile >> temp;
            for (row = 0; row < 4; ++row)
            {
                if (!(ifile >> Elem(row, 0) >> Elem(row, 1) >> Elem(row, 2)))
                {
                    throw std::runtime_error("load Element unsuccessful\n");
                };
            }
            validateElement(Elem, EPS);
            Element->at(i) = Elem;

            Elem(0, 1) = -Elem(0, 1);
            Elem(1, 1) = -Elem(1, 1);
            Elem(2, 1) = -Elem(2, 1);
            Elem(3, 1) = -Elem(3, 1);
            El_middle = Elem.row(1);
            Elem.row(1) = Elem.row(3);
            Elem.row(3) = El_middle;

            Element->at(i + NE_half) = Elem;
        }
    }
    else if (NEType == "fullship")
    {
        for (int i = 0; i < NE; ++i)
        {
            ifile >> temp;
            for (row = 0; row < 4; ++row)
            {
                if (!(ifile >> Elem(row, 0) >> Elem(row, 1) >> Elem(row, 2)))
                    throw std::runtime_error("load Element unsuccessful\n");
            }
            validateElement(Elem, EPS);
            Element->at(i) = Elem;
        }
    }
    else
        throw std::runtime_error("unknow ElementType\n");

    return Element;
}
