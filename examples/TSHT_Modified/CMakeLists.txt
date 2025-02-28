set(EXAMPLE_INC
        include/HTEquation.h
        include/HTNodeData.h
        include/HTAnalytic.h
        include/HTInputData.h
        )

set(EXAMPLE_SRC
        src/TSHT.cpp
        )

add_executable(transient-heat-modified  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(transient-heat-modified PUBLIC include)
target_link_libraries(transient-heat-modified dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
