set(EXAMPLE_INC
#        include/HTEquation.h
#        include/HTNodeData.h
#        include/HTAnalytic.h
#        include/HTInputData.h
        )

set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(dendro-periodic  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(dendro-periodic PUBLIC include)
target_link_libraries(dendro-periodic dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
