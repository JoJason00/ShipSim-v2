# Modules. Each is its own target with explicit dependencies — the dependency graph
# in ShipSim-Refactor/architecture.md is enforced here, not just documented.
#
# Arrows only point downward. If a link line below violates the graph, that is a
# design error, and the build will say so.
#
# Step 2 lands the reference code here; Phase 3 splits it into the modules below.
#
#   core        math / types / constants          (no deps)
#   geometry    panel mesh                        core
#   wave        regular / irregular / spectrum    core
#   hydro       Green function / BEM / kernels    core geometry
#   seakeeping  Cummins time domain / RAO         hydro wave geometry
#   maneuver    MMG                               core wave
#   coupled     two-time-scale coordinator        seakeeping maneuver wave
#   config      case loading + validation         core
#   io          result writing                    core

add_library(shipsim_core
    core/Version.cpp
    core/Version.h
)

target_include_directories(shipsim_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(shipsim_core PRIVATE SHIPSIM_VERSION="${PROJECT_VERSION}")

target_link_libraries(shipsim_core
    PUBLIC  third_party_eigen
    PRIVATE project_warnings
)

if(ENABLE_COVERAGE)
    enable_coverage(shipsim_core)
endif()
