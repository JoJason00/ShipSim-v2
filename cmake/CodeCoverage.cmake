# enable_coverage(<target>) — adds coverage instrumentation.
#   clang-cl     -> llvm-profdata + llvm-cov  (see scripts/coverage.ps1)
#   GCC / Clang  -> gcov + gcovr
#   MSVC         -> unsupported, warns and does nothing

function(enable_coverage target_name)
    get_target_property(target_type ${target_name} TYPE)

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
        target_compile_options(${target_name} PRIVATE -fprofile-instr-generate -fcoverage-mapping)

        if(NOT target_type STREQUAL "STATIC_LIBRARY")
            # link.exe is the linker here, so the clang driver never gets a chance to
            # add its profile runtime. Locate it next to the compiler and link it.
            get_filename_component(_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
            get_filename_component(_root "${_bin}" DIRECTORY)
            file(GLOB _rt
                "${_root}/lib/clang/*/lib/windows/clang_rt.profile-x86_64.lib"
                "${_root}/lib/clang/*/lib/x86_64-pc-windows-msvc/clang_rt.profile.lib")
            if(NOT _rt)
                message(FATAL_ERROR "clang_rt.profile not found near ${CMAKE_CXX_COMPILER}")
            endif()
            list(GET _rt 0 _rt)
            target_link_libraries(${target_name} PRIVATE "${_rt}")
            # That runtime uses the static CRT while we build /MDd; silences LNK4098.
            target_link_options(${target_name} PRIVATE /NODEFAULTLIB:libcmt)
        endif()

    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target_name} PRIVATE -O0 -g --coverage)
        if(NOT target_type STREQUAL "STATIC_LIBRARY")
            target_link_options(${target_name} PRIVATE --coverage)
        endif()

    else()
        message(WARNING "${CMAKE_CXX_COMPILER_ID} cannot instrument; skipping '${target_name}'.")
    endif()
endfunction()
