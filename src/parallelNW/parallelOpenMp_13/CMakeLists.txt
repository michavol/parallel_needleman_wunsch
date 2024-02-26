set(COMP_OPT_SPECIFIC "")#"-ftree-switch-conversion" "-ftree-vectorize" "-ftree-loop-vectorize" "-ftree-slp-vectorize" "-fvect-cost-model=dynamic")

set(headers
        include/nwHelper.hpp
        include/userFunctions.hpp
)

set(sources
        src/nwHelper.cpp
)

# Find OpenMP
find_package(OpenMP REQUIRED)

if(OPENMP_FOUND)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(parallelOpenMpLib ${sources} ${headers})
target_include_directories(parallelOpenMpLib PUBLIC include)
target_compile_options(parallelOpenMpLib PUBLIC ${COMP_OPT} ${COMP_OPT_SPECIFIC})

add_executable(alg_13 main.cpp)
target_link_libraries(alg_13 parallelOpenMpLib)
target_compile_options(alg_13 PUBLIC ${COMP_OPT})