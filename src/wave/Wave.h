#pragma once

#include <memory>
#include "WaveBase.h"
#include "RegularWave.h"
#include "IrregularWave.h"
#include "CrossWave.h"
#include "../config/WaveConfig.h" 

#include <string>

namespace Wave
{
    std::shared_ptr<WaveBase> createWave(const WaveConfig& cfg);

    // Compact, file-name-safe descriptors that identify an incident sea from the
    // name alone, for every wave type. Shared by the seakeeping run tag and the
    // coupled output folder so both stay consistent. Direction is in whole
    // degrees (exactly what the case file specifies).
    //   single regular   -> "Reg-d{deg}-w{ω}"
    //   single irregular -> "{SPEC}-d{deg}-{params}"   (JON/PM/ITTC/OH)
    std::string systemDescriptor(const std::shared_ptr<WaveBase>& w);

    // Whole-sea tag: a single system's descriptor, or both crossing systems
    // joined as "X_{sys1}__{sys2}" so a crossing condition is identifiable.
    std::string conditionTag(const std::shared_ptr<WaveBase>& w);
}