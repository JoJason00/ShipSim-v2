# 获取所有 cpp 文件
file(GLOB_RECURSE SHIPSIM_SOURCES CONFIGURE_DEPENDS "*.cpp")

add_library(shipsim_core STATIC ${SHIPSIM_SOURCES})


target_compile_features(shipsim_core PUBLIC cxx_std_20)

target_compile_options(shipsim_core PRIVATE
    $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:
        -O3
        -march=native
        -fno-math-errno
        -ffp-contract=fast
    >

    $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:
        /O2
        /Ob2
        /Ot
        /arch:AVX2
        /fp:fast
        /GL
    >
)

target_link_options(shipsim_core PRIVATE
    $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:
        -flto
    >
    $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:
        /LTCG
    >
)

target_compile_definitions(shipsim_core PRIVATE
    $<$<CONFIG:Release>:EIGEN_NO_DEBUG>
)


target_include_directories(shipsim_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/eigen-5.0.0
        ${CMAKE_SOURCE_DIR}/json
)

# [v2] core/Version.cpp — the only source added on top of the reference code.
target_compile_definitions(shipsim_core PRIVATE SHIPSIM_VERSION="${PROJECT_VERSION}")

# jsoncpp
target_sources(shipsim_core PRIVATE
    ${CMAKE_SOURCE_DIR}/jsoncpp.cpp
)

# OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    target_link_libraries(shipsim_core PRIVATE OpenMP::OpenMP_CXX)
endif()

# ------------------------------------------------------------------
# Optional BLAS / LAPACK backend for Eigen large-matrix operations.
# Off by default (uses Eigen's built-in implementation).
#
#   -DSHIPSIM_USE_MKL=ON       use Intel MKL (requires MKLROOT or env)
#   -DSHIPSIM_USE_OPENBLAS=ON  use OpenBLAS  (requires BLAS::BLAS findable)
#
# Both flags simply set Eigen preprocessor switches; no source changes needed.
# ------------------------------------------------------------------
option(SHIPSIM_USE_MKL "Use Intel MKL as Eigen backend" OFF)
option(SHIPSIM_USE_OPENBLAS "Use OpenBLAS as Eigen backend" OFF)

if(SHIPSIM_USE_MKL)
    find_package(MKL CONFIG)
    if(MKL_FOUND)
        target_compile_definitions(shipsim_core PRIVATE EIGEN_USE_MKL_ALL)
        target_link_libraries(shipsim_core PRIVATE MKL::MKL)
        message(STATUS "Eigen backend: Intel MKL")
    else()
        message(WARNING "SHIPSIM_USE_MKL=ON but MKL not found, falling back to default Eigen")
    endif()
elseif(SHIPSIM_USE_OPENBLAS)
    find_package(BLAS)
    find_package(LAPACK)
    if(BLAS_FOUND AND LAPACK_FOUND)
        target_compile_definitions(shipsim_core PRIVATE EIGEN_USE_BLAS EIGEN_USE_LAPACKE)
        target_link_libraries(shipsim_core PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        message(STATUS "Eigen backend: OpenBLAS / LAPACK")
    else()
        message(WARNING "SHIPSIM_USE_OPENBLAS=ON but BLAS/LAPACK not found, falling back to default Eigen")
    endif()
endif()


