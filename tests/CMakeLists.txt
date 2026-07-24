# unit/        单个函数/类,硬编码期望值
# validation/  对标试验或解析解,宽容差(2~5%)
# regression/  对标冻结基线,严容差(1e-10);见 ShipSim-Refactor/testing-strategy.md

add_executable(unit_tests
    unit/test_version.cpp
)

target_link_libraries(unit_tests
    PRIVATE shipsim_core Catch2::Catch2WithMain project_warnings
)

include(Catch)
catch_discover_tests(unit_tests)

if(ENABLE_COVERAGE)
    enable_coverage(unit_tests)
endif()
