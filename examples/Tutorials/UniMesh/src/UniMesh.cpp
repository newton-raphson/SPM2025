//
// Created by dhruv on 5/26/23.
//

#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <UniMeshInputData.h>
#include <IO/VTU.h>
#include "SDARefine.h"

void performNoChangeRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                               DENDRITE_UINT maxLevel, SubDomain & subDomain){
    // Use this function to perform no refinement but to remove the physical boundary octants in case of
    // retain inside case
    SDARefine refine(octDA, distTree.getTreePartFiltered(),domainInfo,maxLevel);
    DA * newDA = refine.getForceRefineSubDA(distTree,0.01);
    std::swap(octDA, newDA);
    delete newDA;
    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);
}

template<typename T>
void CreateMesh(UniMeshInputData &idata, T &geoObj, DENDRITE_UINT level, DENDRITE_UINT eleOrder, std::string geoName){
    /// Physical dimensions.
    DomainInfo cubeDomain, physDomain;
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
        cubeDomain.min[dim] = idata.cubeDomainMin[dim];
        cubeDomain.max[dim] = idata.cubeDomainMax[dim];
        physDomain.min[dim] = idata.physDomainMin[dim];
        physDomain.max[dim] = idata.physDomainMax[dim];
    }
    DomainExtents domainExtents(cubeDomain,physDomain);
    // ********************************************************************************************************************
    // Create STL file carved out domain using the same domainExtents
    SubDomain subDomain(domainExtents);

    subDomain.addObject(geoObj);

    // Retain Function
    auto functionToRetain = [&](const double *physCoords,double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };

    /// Create DA
    DistTREE distTree;
    auto octDA = createSubDA(distTree, functionToRetain, level, eleOrder, 0.3);
    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
    performNoChangeRefinement(octDA, distTree, domainExtents, 5, subDomain);
    IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), geoName.c_str(), geoName.c_str(), domainExtents);
    delete octDA;
}

int main(int argc, char *argv[]) {
  dendrite_init(argc, argv);
  /// read parameters from config.txt
  UniMeshInputData idata;

  /// use config.txt first
  std::ifstream configFile("config.txt");
  DENDRITE_UINT eleOrder = 0;
  DENDRITE_UINT level = 0;
  if (configFile.good()) {
    if (!idata.ReadFromFile()) {  /// read from file named "config.txt"
      throw std::runtime_error("[ERR] Error reading input data, check the config file!");
    }
    if (!idata.CheckInputData()) {
      throw std::runtime_error("[ERR] Problem with input data, check the config file!");
    }
    eleOrder = idata.basisFunction;
    level = idata.refine_lvl;
  } else {
    TALYFEMLIB::PrintStatus("Fix the config file!");
    return -1;
  }

  TALYFEMLIB::PrintInfo("Total number of processor = ", TALYFEMLIB::GetMPISize());
  TALYFEMLIB::PrintInfo("size of DendroInt ", sizeof(DendroIntL));
  TALYFEMLIB::PrintInfo("size of PetscInt ", sizeof(PetscInt));

  TALYFEMLIB::PrintStatus("eleOrder ", eleOrder);
  TALYFEMLIB::PrintStatus("Level ", level);
  TALYFEMLIB::PrintStatus("DIM =  ", DIM);

#if(DIM==2)
    {
        VOXEL::Circle circle(idata.circleCenter, idata.circleRadius, idata.circleRetainSide);
        CreateMesh(idata,circle,level,eleOrder,"circle");
    }

    {
        VOXEL::Box box(idata.rectLowerLeft, idata.rectUpperRight, idata.rectRetainSide);
        CreateMesh(idata,box,level,eleOrder,"rectangle");
    }
    {
        GEOMETRY::MSH msh(idata.mshFileName.c_str());
        Point<DIM> temp(0.0);
        GEOMETRY::Geometry *mshGeo = new GEOMETRY::Geometry(&msh, temp, idata.mshRetainSide);
        CreateMesh(idata,mshGeo,level,eleOrder,"msh");
        delete mshGeo;
    }
#endif
#if(DIM==3)
    {
        // Use CreateMesh to create a mesh with VOXEL::Sphere as carve out object
        VOXEL::Sphere sphere(idata.sphereCenter, idata.sphereRadius, idata.sphereRetainSide);
        CreateMesh(idata, sphere, level, eleOrder, "sphere");
    }

    {
        // Use CreateMesh to create a mesh with VOXEL::Box as carve out object
        VOXEL::Box box(idata.boxLowerLeft, idata.boxUpperRight, idata.boxRetainSide);
        CreateMesh(idata, box, level, eleOrder, "box");
    }

    {
        // Use CreateMesh to create a mesh with GEOMETRY::STL as carve out object
        GEOMETRY::STL stl(idata.stlFileName.c_str());
        Point<DIM> temp2(0.0);
        GEOMETRY::Geometry *stlGeo = new GEOMETRY::Geometry(&stl, temp2, idata.stlRetainSide);
        CreateMesh(idata, stlGeo, level, eleOrder, "stl");
        delete stlGeo;
    }
#endif
    dendrite_finalize();
//    dendrite_finalize(octDA);
//    dendrite_finalize(octDA2);
//    dendrite_finalize(octDA3);
}