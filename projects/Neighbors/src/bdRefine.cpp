//
// Created by dhruv on 5/26/23.
//

#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <bdRefineInputData.h>
#include <IO/VTU.h>
#include "SDARefine.h"
#include "Traversal/GenerateMesh.h"
#include "NeighborSearch.h"

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

void performRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                               DENDRITE_UINT maxLevel, SubDomain & subDomain){
    while(true) {
        SDARefine refine(octDA, distTree.getTreePartFiltered(), domainInfo, maxLevel);

        DA *newDA = refine.getRefineSubDA(distTree, 0.01,RefinementStrategy::FULL_DOMAIN,ot::RemeshPartition::SurrogateInByOut);
        if(newDA == nullptr){
            break;
        }
        std::swap(octDA, newDA);
        delete newDA;
        subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);
    }

}

template<typename T>
void CreateMesh(const DeepTraceInputData &idata, T &geoObj, DENDRITE_UINT uniLevel, DENDRITE_UINT bdLevel, DENDRITE_UINT
eleOrder, std::string geoName){
    /// Physical dimensions.
    DomainInfo cubeDomain, physDomain;
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
        cubeDomain.min[dim] = idata.cubeDomainMin[dim];
        cubeDomain.max[dim] = idata.cubeDomainMax[dim];
        physDomain.min[dim] = idata.physDomainMin[dim];
        physDomain.max[dim] = idata.physDomainMax[dim];
    }
    DomainExtents domainExtents(cubeDomain,physDomain);
    SubDomain subDomain(domainExtents);

//    subDomain.addObject(geoObj);

    // Retain Function
    auto functionToRetain = [&](const double *physCoords,double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };

    /// Create DA
    DistTREE distTree;
    DA * octDA = createSubDA(distTree, functionToRetain, uniLevel, eleOrder, 0.3);
    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
//    performRefinement(octDA, distTree, domainExtents, bdLevel, subDomain);
//    performNoChangeRefinement(octDA, distTree, domainExtents, 5, subDomain);

    SubDomainBoundary boundary(&subDomain,octDA,domainExtents);

//    GenerateMeshInfo meshInfo(octDA,distTree.getTreePartFiltered(),domainExtents,&boundary);
//    meshInfo.generate();
//    meshInfo.print("GMSH/nodes.txt","GMSH/elements.txt","GMSH/boundary.txt");
//    meshInfo.printGMSH("mesh.msh");

    IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), "test", "test", domainExtents);
//    delete octDA;
}

int main(int argc, char *argv[]) {
  dendrite_init(argc, argv);
  /// read parameters from config.txt
  DeepTraceInputData idata;

  /// use config.txt first
  std::ifstream configFile("config.txt");
  DENDRITE_UINT eleOrder = 0;
  DENDRITE_UINT uniLevel = 0;
  DENDRITE_UINT bdLevel = 0;
  if (configFile.good()) {
    if (!idata.ReadFromFile()) {  /// read from file named "config.txt"
      throw std::runtime_error("[ERR] Error reading input data, check the config file!");
    }
    if (!idata.CheckInputData()) {
      throw std::runtime_error("[ERR] Problem with input data, check the config file!");
    }
    eleOrder = idata.basisFunction;
    uniLevel = idata.refine_lvl_uni;
    bdLevel = idata.refine_lvl_bd;
  } else {
    TALYFEMLIB::PrintStatus("Fix the config file!");
    return -1;
  }

  TALYFEMLIB::PrintInfo("Total number of processor = ", TALYFEMLIB::GetMPISize());
  TALYFEMLIB::PrintInfo("size of DendroInt ", sizeof(DendroIntL));
  TALYFEMLIB::PrintInfo("size of PetscInt ", sizeof(PetscInt));

  TALYFEMLIB::PrintStatus("eleOrder ", eleOrder);
  TALYFEMLIB::PrintStatus("Level ", uniLevel);
  TALYFEMLIB::PrintStatus("DIM =  ", DIM);

#if(DIM==2)
    {
        VOXEL::Circle circle(idata.circleCenter, idata.circleRadius, idata.circleRetainSide);
//        CreateMesh(idata, circle, uniLevel, bdLevel, eleOrder, "circle");
    }

    {
        VOXEL::Box box(idata.rectLowerLeft, idata.rectUpperRight, idata.rectRetainSide);
//        CreateMesh(idata,box,uniLevel, bdLevel, eleOrder,"rectangle");
    }
    {
        GEOMETRY::MSH msh(idata.mshFileName.c_str());
        Point<DIM> temp(0.0);
        GEOMETRY::Geometry *mshGeo = new GEOMETRY::Geometry(&msh, temp, idata.mshRetainSide);
//      CreateMesh(idata,mshGeo,uniLevel, bdLevel, eleOrder,"msh");



        DomainInfo cubeDomain, physDomain;
        for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
            cubeDomain.min[dim] = idata.cubeDomainMin[dim];
            cubeDomain.max[dim] = idata.cubeDomainMax[dim];
            physDomain.min[dim] = idata.physDomainMin[dim];
            physDomain.max[dim] = idata.physDomainMax[dim];
        }
        DomainExtents domainExtents(cubeDomain,physDomain);
        SubDomain subDomain(domainExtents);

//    subDomain.addObject(geoObj);

        // Retain Function
        auto functionToRetain = [&](const double *physCoords,double physSize) {
            return (subDomain.functionToRetain(physCoords, physSize));
        };

        /// Create DA
        DistTREE distTree;
        DA * octDA = createSubDA(distTree, functionToRetain, uniLevel, eleOrder, 0.3);
//        subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
        performRefinement(octDA, distTree, domainExtents, bdLevel, subDomain);
//    performNoChangeRefinement(octDA, distTree, domainExtents, 5, subDomain);

        SubDomainBoundary boundary(&subDomain,octDA,domainExtents);

        GenerateMeshInfo meshInfo(octDA,distTree.getTreePartFiltered(),domainExtents,&boundary);
        meshInfo.generate();
//    meshInfo.print("GMSH/nodes.txt","GMSH/elements.txt","GMSH/boundary.txt");
        meshInfo.printGMSH("mesh.msh");
        std::cout << " done!" << std::endl;
        IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), "test", "test", domainExtents);

        NeighborSearch neighborSearch(octDA,distTree.getTreePartFiltered(),domainExtents,&boundary);
        neighborSearch.generate();  // generate neighbor search info



//        delete mshGeo;
    }
#endif
#if(DIM==3)
//    {
//        VOXEL::Sphere sphere(idata.sphereCenter, idata.sphereRadius, idata.sphereRetainSide);
//        CreateMesh(idata,sphere,uniLevel, bdLevel, eleOrder,"sphere");
//    }
//
//    {
//        VOXEL::Box box(idata.boxLowerLeft, idata.boxUpperRight, idata.boxRetainSide);
//        CreateMesh(idata,box,uniLevel, bdLevel, eleOrder,"box");
//    }

    {
        GEOMETRY::STL stl(idata.stlFileName.c_str());
        Point<DIM> temp2(0.0);
        GEOMETRY::Geometry *stlGeo = new GEOMETRY::Geometry(&stl, temp2, idata.stlRetainSide);
//        CreateMesh(idata,stlGeo,uniLevel, bdLevel, eleOrder,"stl");
//        delete stlGeo;

//      CreateMesh(idata,mshGeo,uniLevel, bdLevel, eleOrder,"msh");



        DomainInfo cubeDomain, physDomain;
        for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
            cubeDomain.min[dim] = idata.cubeDomainMin[dim];
            cubeDomain.max[dim] = idata.cubeDomainMax[dim];
            physDomain.min[dim] = idata.physDomainMin[dim];
            physDomain.max[dim] = idata.physDomainMax[dim];
        }
        DomainExtents domainExtents(cubeDomain,physDomain);
        SubDomain subDomain(domainExtents);

//    subDomain.addObject(geoObj);

        // Retain Function
        auto functionToRetain = [&](const double *physCoords,double physSize) {
            return (subDomain.functionToRetain(physCoords, physSize));
        };

        /// Create DA
        DistTREE distTree;
        DA * octDA = createSubDA(distTree, functionToRetain, uniLevel, eleOrder, 0.3);
//        subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
        performRefinement(octDA, distTree, domainExtents, bdLevel, subDomain);
//    performNoChangeRefinement(octDA, distTree, domainExtents, 5, subDomain);

        SubDomainBoundary boundary(&subDomain,octDA,domainExtents);

//        GenerateMeshInfo meshInfo(octDA,distTree.getTreePartFiltered(),domainExtents,&boundary);
//        meshInfo.generate();
//    meshInfo.print("GMSH/nodes.txt","GMSH/elements.txt","GMSH/boundary.txt");
//        meshInfo.printGMSH("mesh.msh");
//        std::cout << " done!" << std::endl;
        IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), "test", "test", domainExtents);

        NodePrint neighborSearch(octDA, distTree.getTreePartFiltered(), domainExtents, &boundary);
        neighborSearch.generate();  // generate neighbor search info















    }
#endif
//    dendrite_finalize();
// std::cout << " done!" << std::endl;
//    return 0;
}