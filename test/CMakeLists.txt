set(PYRAMID_SRC
        pyramidElement.cpp
        )

add_executable(tstpyramid  ${PYRAMID_SRC} )
target_include_directories(tstpyramid PUBLIC include)
target_link_libraries(tstpyramid dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)


set(HEAT_SRC
        heat.cpp
)



set(REFINEMENT_INC
        Refinement/BoundaryRefine.h
        RefinementUtils.h
)

set(MVEC_INC
        HeatEquation/SSHTNodeData.h
        HeatEquation/SSHTEquation.h)

add_executable(tstheat  ${HEAT_SRC})
target_include_directories(tstheat PUBLIC include)
target_link_libraries(tstheat dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)