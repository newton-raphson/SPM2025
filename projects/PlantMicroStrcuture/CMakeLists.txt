
set(LE_INC
        include/util.h
        include/LEInputData.h
        include/LENodeData.h
        include/LEBCSetup.h
        include/SSLEEquation.h
        include/LEInputDataStructs.h
        include/MaterialMap2D.h
        include/stb_image.h
)

set(SLE_SRC src/SSLE.cpp
        include/CalcStress.h)

add_executable(plant-sle ${SLE_SRC} ${LE_INC})
target_include_directories(plant-sle PUBLIC include)
target_link_libraries(plant-sle dendriteKT dendroKT ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
