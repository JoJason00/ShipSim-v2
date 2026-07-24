#pragma once

#include <string>

namespace shipsim {

// Version string of the build, e.g. "2.0.0". Written into result directories so a
// set of results can always be traced back to the code that produced it.
std::string version();

} // namespace shipsim
