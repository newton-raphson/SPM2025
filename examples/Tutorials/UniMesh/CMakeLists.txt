set(EXAMPLE_INC
        include/UniMeshInputData.h
        )

set(EXAMPLE_SRC
        src/UniMesh.cpp
        )

add_executable(tut_UniMesh ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(tut_UniMesh PUBLIC include)
target_link_libraries(tut_UniMesh dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
