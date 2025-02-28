set(EXAMPLE_INC
        include/SSHTEquation.h
        include/SSHTNodeData.h
        include/SSHTInputData.h
        include/UniMeshInputData.h
        include/SDARefine.h
        )

set(EXAMPLE_SRC
        src/SSHT.cpp
        )

add_executable(steady-heat-stl  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(steady-heat-stl PUBLIC include)
target_link_libraries(steady-heat-stl dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
