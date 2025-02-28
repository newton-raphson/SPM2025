#!/bin/bash
cd ../../

mkdir build_2d
cd build_2d
( cmake .. -DENABLE_2D=Yes -DCMAKE_BUILD_TYPE=Release
-DCMAKE_C_COMPILEER=mpicc -DCMAKE_CXX_COMPILER=mpicxx )
make tut_UniMesh -j6
cd ..

mkdir build_3d
cd build_3d
( cmake .. -DENABLE_3D=Yes -DCMAKE_BUILD_TYPE=Release
-DCMAKE_C_COMPILEER=mpicc -DCMAKE_CXX_COMPILER=mpicxx )
make tut_UniMesh -j6
cd ..