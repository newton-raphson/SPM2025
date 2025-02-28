set(EXAMPLE_INC
        include/SSHTNodeData.h
        include/SSHTNodeData.h
        include/ElementMarker.h
        include/BFS.h
        include/utils.h
        include/IMGALoop.h
        include/CheckSurface.h
        include/DACoarsen.h
        )

set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(sbm-test  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(sbm-test PUBLIC include)
target_link_libraries(sbm-test dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
