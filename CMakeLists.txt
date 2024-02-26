cmake_minimum_required(VERSION 3.1)
project(SequenceAlignment VERSION 0.1.0)

# Activate c++11 flag
set(CMAKE_CXX_STANDARD 11)
set(COMP_OPT "-O3" "-mavx" "-mavx2" "-fopenmp" "-march=skylake")
#set(CMAKE_CXX_FLAGS "-O0 -Wall")

# set default build type to Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif ()

include(CTest)
enable_testing()

add_subdirectory(src/sequentialNW/sequentialBaseline_00)
add_subdirectory(src/sequentialNW/sequentialOpt_01)
add_subdirectory(src/sequentialNW/sequentialDiag_02)
add_subdirectory(src/sequentialNW/sequentialDiagOpt_03)
add_subdirectory(src/sequentialNW/sequentialDiagSimd_04)
add_subdirectory(src/sequentialNW/sequentialDiagOptNested_05)
add_subdirectory(src/sequentialNW/sequentialDiagSimdPragma_06)
add_subdirectory(src/parallelNW/parallelBaseline_10)
add_subdirectory(src/parallelNW/parallelBaselineOpt_11)
add_subdirectory(src/parallelNW/parallelShared_12)
add_subdirectory(src/parallelNW/parallelOpenMp_13)
add_subdirectory(src/parallelNW/parallelOpenMpOpt_14)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
