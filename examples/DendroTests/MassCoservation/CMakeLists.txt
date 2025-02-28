set(EXAMPLE_SRC
        src/main.cpp
        )

set(EXAMPLE_INC
        include/MassComputation.h
        include/Refinement.h
        include/utils.h
        include/MassCalculatorTotal.h

        include/MEquation.h
        include/MNodeData.h
        )


add_executable(mconv  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(mconv PUBLIC include)
target_link_libraries(mconv dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
