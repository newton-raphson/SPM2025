set(EXAMPLE_INC
        include/NSEquation.h
        include/NSNodeData.h
        include/NSRefine.h
        include/NSUtils.h
        include/NSInputData.h
        include/NSBoundaryConditions.h
        )

set(EXAMPLE_SRC
        src/main.cpp
        include/NSPost.h
)

add_executable(ns-ld-cons  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(ns-ld-cons PUBLIC include)
target_link_libraries(ns-ld-cons dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
