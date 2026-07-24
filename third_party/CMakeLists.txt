# Vendored dependencies, carried over unchanged from the reference implementation.
#
# These stay vendored on purpose: the migration changes one variable at a time, and
# swapping Eigen/jsoncpp for vcpkg packages happens only AFTER the regression
# baselines are green (see ShipSim-Refactor/README.md, Step 5).
#
# Populated in Step 2. Until then this directory only declares the interface targets
# so the rest of the build can already link against them.

# --- Eigen (header-only) ---
add_library(third_party_eigen INTERFACE)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/eigen-5.0.0/Eigen)
    target_include_directories(third_party_eigen SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/eigen-5.0.0)
endif()

# --- jsoncpp (single translation unit) ---
add_library(third_party_json INTERFACE)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/jsoncpp/jsoncpp.cpp)
    target_sources(third_party_json INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/jsoncpp/jsoncpp.cpp)
    target_include_directories(third_party_json SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/jsoncpp)
endif()
