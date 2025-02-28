set(EXAMPLE_INC
        include/bdRefineInputData.h
        include/SDARefine.h
        )

set(EXAMPLE_SRC
        src/bdRefine.cpp
        )

add_executable(tut_bdRefine ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(tut_bdRefine PUBLIC include)
target_link_libraries(tut_bdRefine dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
