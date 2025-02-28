set(EXAMPLE_SRC
        src/main.cpp
        )

add_executable(cpoint  ${EXAMPLE_SRC})
target_include_directories(cpoint PUBLIC include)
target_link_libraries(cpoint dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)

set(EXAMPLE_SRC
        src/subDABug.cpp
        )

add_executable(sbug  ${EXAMPLE_SRC})
target_include_directories(sbug PUBLIC include)
target_link_libraries(sbug dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)

set(EXAMPLE_SRC
        src/IntergridBug.cpp
        )

add_executable(intergridCheckpoint  ${EXAMPLE_SRC})
target_include_directories(intergridCheckpoint PUBLIC include)
target_link_libraries(intergridCheckpoint dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
