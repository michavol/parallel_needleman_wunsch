set(COMP_OPT_SPECIFIC "")#"-ftree-switch-conversion" "-ftree-vectorize" "-ftree-loop-vectorize" "-ftree-slp-vectorize" "-fvect-cost-model=dynamic")

set(headers
        include/nwHelper.hpp
        include/nwSimd.hpp
        include/nwSolver.hpp
        include/userFunctions.hpp
)

set(sources
        src/nwHelper.cpp
        src/nwSimd.cpp
        src/nwSolver.cpp 
)

add_library(sequentialDiagSimdLib ${sources} ${headers})
target_include_directories(sequentialDiagSimdLib PUBLIC include)
target_compile_options(sequentialDiagSimdLib PUBLIC ${COMP_OPT} ${COMP_OPT_SPECIFIC})

add_executable(alg_04 main.cpp)
target_link_libraries(alg_04 sequentialDiagSimdLib)
target_compile_options(alg_04 PUBLIC ${COMP_OPT})