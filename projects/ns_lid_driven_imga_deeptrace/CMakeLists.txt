set(ONNX_RUNTIME_LIB_DIR "/home/pc-14-2/PACKAGES/onnxruntime-linux-x64-1.17.1/lib")
include_directories("/home/pc-14-2/PACKAGES/onnxruntime-linux-x64-1.17.1/include")

set(EXAMPLE_INC
        include/NSNodeData.h
        include/NSRefine.h
        include/NSUtils.h
        include/NSInputData.h
        include/NSBoundaryConditions.h
        )

set(EXAMPLE_SRC
        src/main.cpp
        include/CheckSurface.h
        include/IMGALoop.h
        include/NSPost.h
        include/BoundaryCalc_ExtraInter.h
        include/Extrapolation.h

        include/DeepTrace.h
        include/SBMCalcDeepTrace.h
)

add_executable(ns-ld-cons-imga-deeptrace  ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(ns-ld-cons-imga-deeptrace PUBLIC include)
target_link_libraries(ns-ld-cons-imga-deeptrace  dendriteKT dendroKT ${ONNX_RUNTIME_LIB_DIR}/libonnxruntime.so  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES}  m)
