cmake_minimum_required(VERSION 3.21)

project(shipsim VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Sources carry Chinese comments; MSVC needs to be told they are UTF-8.
# Carried over from the reference implementation unchanged.
if(MSVC)
    add_compile_options(/utf-8)
endif()

option(ENABLE_COVERAGE "Build with coverage instrumentation" OFF)
option(ENABLE_ASAN "Build with AddressSanitizer (+UBSan on GCC/Clang)" OFF)

enable_testing()

if(ENABLE_ASAN)
    if(MSVC)
        add_compile_options(/fsanitize=address)
        string(REPLACE "/RTC1" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
        add_compile_definitions(_DISABLE_STRING_ANNOTATION _DISABLE_VECTOR_ANNOTATION)
    else()
        add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
        add_link_options(-fsanitize=address,undefined)
    endif()
endif()

# Ninja picks the compiler off PATH, so a shell without MSVC silently selects MinGW
# while vcpkg still serves MSVC-built libraries. Fail here with a readable message.
if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND VCPKG_TARGET_TRIPLET MATCHES "^x64-windows")
    message(FATAL_ERROR
        "ABI mismatch: MinGW GCC with the '${VCPKG_TARGET_TRIPLET}' triplet (MSVC-built).\n"
        "Delete ${CMAKE_BINARY_DIR} (the failed cache pins the compiler), then configure "
        "again from a Developer PowerShell.")
endif()

find_package(Catch2 3 CONFIG REQUIRED)
include(CodeCoverage)  # provides enable_coverage()

# Warnings live in an opt-in target. shipsim_core deliberately does NOT link it yet:
# the reference build had no /W4, and turning it on across 33k legacy lines would
# bury real problems under thousands of messages. New modules should link it.
add_library(project_warnings INTERFACE)
if(MSVC)
    target_compile_options(project_warnings INTERFACE /W4 /permissive-)
else()
    target_compile_options(project_warnings INTERFACE -Wall -Wextra -Wpedantic)
endif()

add_subdirectory(src)     # shipsim_core — the ported reference code
add_subdirectory(app)     # shipsim_cli
add_subdirectory(test)    # legacy driver programs (core_test, test_green, ...)
add_subdirectory(tests)   # new Catch2 suite: unit / validation / regression
