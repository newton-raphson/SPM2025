set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(sda  ${EXAMPLE_SRC})
target_include_directories(sda PUBLIC include)
target_link_libraries(sda dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)

set(EXAMPLE_SRC
        src/boundary.cpp
        )

add_executable(sda-bnd  ${EXAMPLE_SRC})
target_include_directories(sda-bnd PUBLIC include)
target_link_libraries(sda-bnd dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)