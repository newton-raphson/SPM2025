set(NS_INC
        include/NSPNPNodeData.hpp
        include/NSRefine.h
        include/NSUtils.h
        include/NSPNPInputData.h
        include/NSBoundaryConditions.h
        include/PressurePoissonEquation.h
        include/NSParams.h
        include/NSEquation.h
        include/NSPNPExpr.hpp
        include/PNPEquation.h
        )

set(NS_SRC src/main.cpp)

add_executable(ns-pnp-electroosmotic  ${NS_SRC} ${NS_INC})
target_include_directories(ns-pnp-electroosmotic PUBLIC include)
target_link_libraries(ns-pnp-electroosmotic dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
