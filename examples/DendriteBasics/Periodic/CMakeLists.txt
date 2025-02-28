
set(EXAMPLE_INC
        include/PRefine.h
        include/ADEquation.h
        include/ADNodeData.h

        )


set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(periodic ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(periodic PUBLIC include)
target_link_libraries(periodic dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
