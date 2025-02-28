
set(EXAMPLE_INC
    include/DARefine.h
#    include/DARefineValues.h
        )


set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(intergrid-many ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(intergrid-many PUBLIC include)
target_link_libraries(intergrid-many dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
