set(COMP_OPT_SPECIFIC "-Wno-c++11-extensions") #"-fstrict-aliasing")#"-ftree-switch-conversion" "-ftree-vectorize" "-ftree-loop-vectorize" "-ftree-slp-vectorize" "-fvect-cost-model=dynamic")


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
add_library(baseParOptNWLib ${sources} ${headers})
target_link_libraries(baseParOptNWLib ${MPI_LIBRARIES})
target_include_directories(baseParOptNWLib PUBLIC include)
target_compile_options(baseParOptNWLib PUBLIC ${COMP_OPT} ${COMP_OPT_SPECIFIC})

# Create executable
add_executable(alg_11 main.cpp)
target_link_libraries(alg_11 baseParOptNWLib)
target_compile_options(alg_11 PUBLIC ${COMP_OPT} ${COMP_OPT_SPECIFIC})

