set(COMP_OPT_SPECIFIC "")#"-ftree-switch-conversion" "-ftree-vectorize" "-ftree-loop-vectorize" "-ftree-slp-vectorize" "-fvect-cost-model=dynamic")

set(headers
        include/nwHelper.hpp
        include/nwSolver.hpp
        include/userFunctions.hpp
)

set(sources
        src/nwHelper.cpp
        src/nwSolver.cpp 
)

add_library(sequentialDiagOptNestedLib ${sources} ${headers})
target_include_directories(sequentialDiagOptNestedLib PUBLIC include)
target_compile_options(sequentialDiagOptNestedLib PUBLIC ${COMP_OPT} ${COMP_OPT_SPECIFIC})

add_executable(alg_05 main.cpp)
target_link_libraries(alg_05 sequentialDiagOptNestedLib)
target_compile_options(alg_05 PUBLIC ${COMP_OPT})