set(EXAMPLE_INC
        include/SDARefine.h
        include/BenchDataTypes.h
        )

set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(bench-construction  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(bench-construction PUBLIC include)
target_link_libraries(bench-construction dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
