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
        include/SBMMarker.h
        include/BFS.h
        include/NSDataStructure.h
        include/SBMCalc.h
        include/nanoflann.hpp
)

add_executable(lin-ns-ld-cons-imga  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(lin-ns-ld-cons-imga PUBLIC include)
target_link_libraries(lin-ns-ld-cons-imga dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
