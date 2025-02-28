#include <iostream>
#include <distTree.h>
#include <oda.h>
#include <point.h>
#include <sfcTreeLoop_matvec_io.h>

#include <octUtils.h>
#include <PETSc/IO/petscVTU.h>
#include <DendriteUtils.h>
#include <SDARefine.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include "SSHTEquation.h"
#include <SubDA/SubDomain.h>
#include <Boundary/SubDomainBoundary.h>

using namespace PETSc;



int main(int argc, char * argv[]){
  dendrite_init(argc,argv);
  DENDRITE_UINT eleOrder = 1;
  DENDRITE_UINT ndof = 1;
  m_uiMaxDepth = 25;

  MPI_Comm comm = MPI_COMM_WORLD;
#if(DIM == 3)
  DomainInfo cubeDomain;
  cubeDomain.min.fill(0.0);
  cubeDomain.max.fill(4.0);

  DomainInfo physDomain;
  physDomain.min.fill(0.0);
  physDomain.max.fill(2.0);

#endif
#if(DIM == 2)
  DomainInfo cubeDomain;
  cubeDomain.min.fill(0.0);
  cubeDomain.max.fill(4.0);

  DomainInfo physDomain;
  physDomain.min.fill(0.0);
  physDomain.max[0] = 4.0;
  physDomain.max[1] = 2.0;
#endif


  DomainExtents domainExtents(cubeDomain, physDomain);
  SubDomain subDA(domainExtents);

#if(DIM == 3)
  double centerCoords[DIM] {physDomain.max[0]/(2),physDomain.max[1]/(2),physDomain.max[2]/(2)};
//  double _centerCoords[DIM];
//  static constexpr double RADIUS = 0.45;
//  static double centerSphere[3][3]{{0.5, 0.5, 8.0},
//                                   {0.5, 0.0, 8.0 + RADIUS * 2 + 0.01},
//                                   {0.5, 1.0, 8.0 + RADIUS * 2 + 0.01}};
//  _centerCoords[0] = centerSphere[0][0] ;
//  _centerCoords[1] = centerSphere[0][1] ;
//  _centerCoords[2] = centerSphere[0][2] ;
  subDA.addObject(VOXEL::Sphere(centerCoords,0.45));

//  _centerCoords[0] = centerSphere[1][0] ;
//  _centerCoords[1] = centerSphere[1][1] ;
//  _centerCoords[2] = centerSphere[1][2] ;
//  subDA.addObject(VOXEL::Sphere(_centerCoords,0.45));
//
//  _centerCoords[0] = centerSphere[2][0] ;
//  _centerCoords[1] = centerSphere[2][1] ;
//  _centerCoords[2] = centerSphere[2][2] ;
//  subDA.addObject(VOXEL::Sphere(_centerCoords,0.45));

//  for(int l = 0; l < 2; l++){
//  for(int i = 0; i < 2; i++) {
//
//    _centerCoords[0] = centerCoords[0];
//    _centerCoords[1] = centerCoords[1] - 1.2 * i ;//+ (l%2)*0.6;
//    _centerCoords[2] = centerCoords[2] + l*1.2;
//    if (_centerCoords[1] - 0.5 > 0) {
//      subDA.addObject(VOXEL::Sphere(_centerCoords, 0.45));
//    }
//
//
//    _centerCoords[0] = centerCoords[0];
//    _centerCoords[1] = centerCoords[1] + 1.2 * i ;//+ (l%2)*0.6;
//    _centerCoords[2] = centerCoords[2] + l*1.2;
//    if (_centerCoords[1] + 0.5 < physDomain.max[1]) {
//      subDA.addObject(VOXEL::Sphere(_centerCoords, 0.45));
//    }
//  }
//
//
//  }

  std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
    return (subDA.functionToRetain(physCoords, physSize));
  };

  DistTREE distTree;

  DA * octDA = createSubDA(distTree,functionToRetain,baseLevel,eleOrder);
  const auto & treePart = distTree.getTreePartFiltered();
  subDA.finalize(octDA,treePart,domainExtents);
  std::cout << "Num Local Element = " << octDA->getLocalElementSz() << "\n";
  IO::writeBoundaryElements(octDA,treePart,"subDA","subDA",domainExtents);
  MPI_Barrier(MPI_COMM_WORLD);


  DENDRITE_UINT count = 0;

  for(count = 0; count < 1; count++){
    SDARefine refine(octDA, treePart,domainExtents,subDA, blockLevel, objectLevel);
    DA * newDA = refine.getRefineSubDA(distTree);
    if(newDA == nullptr){
      break;
    }
    std::swap(octDA, newDA);
    delete newDA;

    TALYFEMLIB::PrintInfo(count, " Iteration complete"," ",treePart.size()," ", octDA->getLocalElementSz());
    IO::writeBoundaryElements(octDA,treePart,("RefineSubDA"+std::to_string(count)).c_str(),"subDA",domainExtents);

    count++;
      subDA.finalize(octDA,treePart,domainExtents);

  }

  IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), "BoundaryRefined", "Boundary", domainExtents);
  checkHangingNodes(octDA,distTree.getTreePartFiltered(),domainExtents);
//  checkHangingNodes(octDA,treePart,domainExtents);
  auto sshtEq = new TalyEquation<SSHTEquation, SSHTNodeData>(octDA,treePart, domainExtents);
  LinearSolver *sshtSolver = setLinearSolver(sshtEq, octDA, ndof, false);

  /// Boundary condition
  sshtSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
    Boundary b;
    b.addDirichlet(0, 10);
    return b;
  });

  /// Solve
  sshtSolver->solve();

  const char*varname[]{"T"};
  /// Print files
  petscVectopvtu(octDA, treePart,sshtSolver->getCurrentSolution(), "Solution","ssht", varname, domainExtents, false, false, ndof);
  dendrite_finalize(octDA);
#endif
#if(DIM == 2)
  std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
    return (subDA.functionToRetain(physCoords, physSize));
  };

  DistTREE dtree;

  DA * octDA = createSubDA(dtree,functionToRetain,2,eleOrder);
  std::vector<std::size_t > bdy_index;

  octDA->getBoundaryNodeIndices(bdy_index);
  double coords[DIM];
  for (int i = 0; i < bdy_index.size(); i++) {
    ot::treeNode2Physical(octDA->getTNCoords()[bdy_index[i] + octDA->getLocalNodeBegin()], octDA->getElementOrder(), coords);
    std::cout << coords[0]*4 << " " << coords[1]*4 << "\n";
  }
  subDA.finalize(octDA,dtree.getTreePartFiltered(),domainExtents);
  IO::writeBoundaryElements(octDA,dtree.getTreePartFiltered(),"boundary","boundary",domainExtents);
  dendrite_finalize(octDA);
#endif

}
