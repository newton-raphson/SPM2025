

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_PREFIX_PATH /home/mehdi/Packages/onnxruntime-linux-x64-1.16.2/include)
set(ONNX_RUNTIME_LIB_DIR "/home/mehdi/Packages/onnxruntime-linux-x64-1.16.2/lib")
include_directories("/home/mehdi/Packages/onnxruntime-linux-x64-1.16.2/include")

#set(TORCH_LIB_DIR "/home/mehdi/Packages/libtorch-shared-with-deps-2.1.1+cu118/libtorch/lib")
#include_directories("/home/mehdi/Packages/libtorch-shared-with-deps-2.1.1+cu118/libtorch/include")

#set(CMAKE_PREFIX_PATH /home/mehdi/Packages/libtorch)
#find_package(Torch REQUIRED)

# Link against the ONNX Runtime library
#target_link_libraries(Inference ${ONNX_RUNTIME_LIB_DIR}/libonnxruntime.so) # or .a or .dll



set(EXAMPLE_INC
        include/DeepTraceInputData.h
        include/SDARefine.h
        )

set(EXAMPLE_SRC
        src/deeptrace.cpp
        )

add_executable(deeptrace ${EXAMPLE_SRC} ${EXAMPLE_INC})
target_include_directories(deeptrace PUBLIC include)
target_link_libraries(deeptrace dendriteKT dendroKT ${ONNX_RUNTIME_LIB_DIR}/libonnxruntime.so ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
#