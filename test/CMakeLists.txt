add_executable(core_test
    core_test.cpp)

 add_executable(test_green
    test_green.cpp)

 add_executable(build_tdgf_table
    build_tdgf_table.cpp)

target_link_libraries(core_test
    PRIVATE
        shipsim_core
)

target_link_libraries(test_green
    PRIVATE
        shipsim_core
)

target_link_libraries(build_tdgf_table
    PRIVATE
        shipsim_core
)
