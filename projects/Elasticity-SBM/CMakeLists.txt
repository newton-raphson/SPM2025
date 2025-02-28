include_directories(${ONNX_RUNTIME_LIB_DIR}/../include)

set(LE_INC
        include/util.h
        include/LEInputData.h
        include/LENodeData.h
        include/LEBCSetup.h
        include/TSLEEquation.h
        include/SSLEEquation.h
        include/LERefine.h
        include/LEInputDataStructs.h
)

set(SLE_SRC src/SSLE.cpp include/SBMCalc.h include/GetSurfaceGp.h include/GetTrueSurfaceGP.h include/BFS.h include/ElementMarker.h include/CheckSurface.h
        include/CalcStress.h
        include/DARefine.h
        include/DACoarse.h
        include/SBMCalcDeepTrace.h
        include/DeepTrace.h)

add_executable(sle-sbm ${SLE_SRC} ${LE_INC})
target_include_directories(sle-sbm PUBLIC include)
target_link_libraries(sle-sbm dendriteKT dendroKT ${ONNX_RUNTIME_LIB_DIR}/libonnxruntime.so ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)

#add_executable(sle-sbm  ${EXAMPLE_SRC} ${EXAMPLE_INC})
#target_include_directories(sle-sbm  PUBLIC include)
#target_link_libraries(sle-sbm   dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
