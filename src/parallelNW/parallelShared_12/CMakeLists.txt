cmake_minimum_required(VERSION 3.1)

#project(shared_memory_NW)

set(COMP_OPT_SPECIFIC "-Wno-c++11-extensions") # "-fstrict-aliasing")#"-ftree-switch-conversion" "-ftree-vectorize" "-ftree-loop-vectorize" "-ftree-slp-vectorize" "-fvect-cost-model=dynamic")


set(headers
        include/nwHelper.hpp
        include/nwMpi.hpp
)

set(sources
        src/nwHelper.cpp
        src/userFunctions.cpp
        src/nwMpi.cpp
)

# Find MPI
find_package(MPI REQUIRED)

# Allow finding directories
include_directories(${MPI_INCLUDE_PATH})

# Build library
add_library(ParNW2 ${sources} ${headers})
target_link_libraries(ParNW2 ${MPI_LIBRARIES})
target_include_directories(ParNW2 PUBLIC include)
target_compile_options(ParNW2 PUBLIC ${COMP_OPT} ${COMP_OPT_SPECIFIC})

# Create executable
add_executable(alg_12 main.cpp)
target_link_libraries(alg_12 ParNW2)
target_compile_options(alg_12 PUBLIC ${COMP_OPT} ${COMP_OPT_SPECIFIC})

