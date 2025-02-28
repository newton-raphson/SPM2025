set(EXAMPLE_INC
        include/bdRefineInputData.h
        include/SDARefine.h
        )

set(EXAMPLE_SRC
        src/bdRefine.cpp
        include/NeighborSearch.h
)

add_executable(Neighbors ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(Neighbors PUBLIC include)
target_link_libraries(Neighbors dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
