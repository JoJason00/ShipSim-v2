add_executable(shipsim_cli
    cli_main.cpp)

target_link_libraries(shipsim_cli
    PRIVATE
        shipsim_core
)
