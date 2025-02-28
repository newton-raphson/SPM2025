//
// Created by chenghau on 2/1/24.
//

#include <PETSc/PetscUtils.h>

#ifndef DENDRITEKT_NSDATASTRUCTURE_H
#define DENDRITEKT_NSDATASTRUCTURE_H

struct GeomInfo{
    std::vector<GEOMETRY::Geometry *> carved_geoms;
#if (DIM == 3)
    std::vector<GEOMETRY::STL *> stls; // 3D
#endif
#if (DIM == 2)
    std::vector<GEOMETRY::MSH *> mshs; // 2D
#endif

};

struct BoundaryDistanceData {
    unsigned int distanceArraySize;
    double* distanceToBdry;
    double* bdryNormal;
    double* surfaceNodeCoords;
    bool* distanceSet;

    /**
     *  Store the distance to boundary vectors and geometry normals
     *  This size is too big, but I'm not sure how to determine the correct size.
     * @param octDA
     * @param dim
     */
    BoundaryDistanceData(DA* octDA, unsigned int dim) {
        distanceArraySize = octDA->getLocalElementSz() * octDA->getNumNodesPerElement();

        distanceToBdry = new double[dim * distanceArraySize];
        bdryNormal = new double[dim * distanceArraySize];
        surfaceNodeCoords = new double[dim * distanceArraySize];
        distanceSet = new bool[distanceArraySize];

        for (unsigned int i = 0; i < distanceArraySize; i++) {
            distanceSet[i] = false;
        }
    }

    // Destructor to free the allocated memory
    ~BoundaryDistanceData() {
        delete[] distanceToBdry;
        delete[] bdryNormal;
        delete[] surfaceNodeCoords;
        delete[] distanceSet;
    }

    // Prevent copying
    BoundaryDistanceData(const BoundaryDistanceData&) = delete;
    BoundaryDistanceData& operator=(const BoundaryDistanceData&) = delete;

    // Enable move semantics
    BoundaryDistanceData(BoundaryDistanceData&& other) noexcept
            : distanceArraySize(other.distanceArraySize),
              distanceToBdry(other.distanceToBdry),
              bdryNormal(other.bdryNormal),
              surfaceNodeCoords(other.surfaceNodeCoords),
              distanceSet(other.distanceSet) {
        other.distanceToBdry = nullptr;
        other.bdryNormal = nullptr;
        other.surfaceNodeCoords = nullptr;
        other.distanceSet = nullptr;
    }

    BoundaryDistanceData& operator=(BoundaryDistanceData&& other) noexcept {
        if (this != &other) {
            delete[] distanceToBdry;
            delete[] bdryNormal;
            delete[] surfaceNodeCoords;
            delete[] distanceSet;

            distanceArraySize = other.distanceArraySize;
            distanceToBdry = other.distanceToBdry;
            bdryNormal = other.bdryNormal;
            surfaceNodeCoords = other.surfaceNodeCoords;
            distanceSet = other.distanceSet;

            other.distanceToBdry = nullptr;
            other.bdryNormal = nullptr;
            other.surfaceNodeCoords = nullptr;
            other.distanceSet = nullptr;
        }
        return *this;
    }
};


#endif //DENDRITEKT_NSDATASTRUCTURE_H
