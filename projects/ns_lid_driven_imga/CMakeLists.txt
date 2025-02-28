if (DEEPTRACE)
        set(ONNX_RUNTIME_LIB_DIR "/work/mech-ai-scratch/samundra/packages/onnxruntime-linux-x64-1.17.1/lib")
        include_directories("/work/mech-ai-scratch/samundra/packages/onnxruntime-linux-x64-1.17.1/include")
endif ()

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
        include/CheckSurface.h
        include/IMGALoop.h
        include/NSPost.h
        include/BoundaryCalc_ExtraInter.h
        include/Extrapolation.h
        include/DARefine.h
        include/DACoarsen.h
#        These are set for deeptrace here
        include/DeepTrace.h
        include/SBMCalcDeepTrace.h
        include/DeepTrace.h
)
if(DEEPTRACE)
        add_executable(ns-ld-cons-imga  ${EXAMPLE_SRC} ${EXAMPLE_INC})
        target_include_directories(ns-ld-cons-imga PUBLIC include)
        target_link_libraries(ns-ld-cons-imga  dendriteKT dendroKT ${ONNX_RUNTIME_LIB_DIR}/libonnxruntime.so  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES}  m)
else ()

        add_executable(ns-ld-cons-imga  ${EXAMPLE_SRC} ${EXAMPLE_INC})
        target_include_directories(ns-ld-cons-imga PUBLIC include)
        target_link_libraries(ns-ld-cons-imga dendriteKT dendroKT  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
endif ()
