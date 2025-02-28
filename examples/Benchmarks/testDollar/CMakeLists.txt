set(EXAMPLE_INC

        )

set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(tst-dollar  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(tst-dollar  PUBLIC include)
target_link_libraries(tst-dollar  dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
