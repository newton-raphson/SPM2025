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
        )

add_executable(ns-vms  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(ns-vms PUBLIC include)
target_link_libraries(ns-vms dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
