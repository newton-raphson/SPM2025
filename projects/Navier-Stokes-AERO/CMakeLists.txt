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
        include/Post.h
)

add_executable(ns-aero  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(ns-aero PUBLIC include)
target_link_libraries(ns-aero dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
