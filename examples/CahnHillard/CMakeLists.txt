set(CH_INC
        include/CHEquation.h
        include/CHNodeData.h
        include/CHInputData.h
        include/CHInitialCondition.h
        include/CHRefine.h
        include/CHUtils.h
        include/MassComputation.h
        include/MEquation.h
        include/MergeChildrenMass.h
        )

set(CH_SRC
        src/CH.cpp
        )

add_executable(ch  ${CH_SRC} ${CH_INC})
target_include_directories(ch PUBLIC include)
target_link_libraries(ch dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
