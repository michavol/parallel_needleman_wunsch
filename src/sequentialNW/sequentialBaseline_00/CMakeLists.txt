set(COMP_OPT_SPECIFIC "") #"-fstrict-aliasing")

set(headers
        include/userFunctions.hpp
        include/nwHelper.hpp
        include/nwSolver.hpp
)

set(sources
        src/userFunctions.cpp
        src/nwHelper.cpp
        src/nwSolver.cpp 
)

add_library(sequentialBaselineLib ${sources} ${headers})
target_include_directories(sequentialBaselineLib PUBLIC include)
target_compile_options(sequentialBaselineLib PUBLIC ${COMP_OPT} ${COMP_OPT_SPECIFIC})

add_executable(alg_00 main.cpp)
target_link_libraries(alg_00 sequentialBaselineLib)
target_compile_options(alg_00 PUBLIC ${COMP_OPT})