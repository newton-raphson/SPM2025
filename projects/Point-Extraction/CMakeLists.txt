set(EXAMPLE_INC
        include/ptextractInputData.h
        include/SDARefine.h
        )

set(EXAMPLE_SRC
        src/ptextract.cpp
        )

add_executable(ptextract ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(ptextract PUBLIC include)
target_link_libraries(ptextract dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
