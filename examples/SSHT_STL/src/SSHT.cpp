//
// Created by maksbh on 9/18/19.
//

#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <SSHTNodeData.h>
#include <SSHTEquation.h>
#include <SSHTInputData.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <PETSc/IO/petscVTU.h>
#include <Traversal/Analytic.h>
#include <UniMeshInputData.h>
#include <SDARefine.h>

using namespace PETSc;





void performNoChangeRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                               DENDRITE_UINT maxLevel, SubDomain & subDomain){
    // Use this function to perform no refinement but to remove the physical boundary octants in case of
    // retain inside case
    SDARefine refine(octDA, distTree.getTreePartFiltered(),domainInfo);
    DA * newDA = refine.getForceRefineSubDA(distTree,0.01);
    std::swap(octDA, newDA);
    delete newDA;
    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);
}

//template<typename T>
//void CreateMesh(UniMeshInputData &idata, DomainExtents &domainExtents, SubDomain &subDomain, DistTREE &distTree, DA *octDA, DENDRITE_UINT level, DENDRITE_UINT eleOrder, std::string geoName){
//
//    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
//    performNoChangeRefinement(octDA, distTree, domainExtents, 5, subDomain);
//    IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), geoName.c_str(), geoName.c_str(), domainExtents);
////    delete octDA;
//}


template<typename T>
void CreateMesh(SSHTInputData &idata, T &geoObj, DENDRITE_UINT level, DENDRITE_UINT eleOrder, std::string geoName){
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
//  dendrite_init(argc, argv);
//  /// read parameters from config.txt
//  SSHTInputData idata;
//
//  /// use config.txt first
//  std::ifstream configFile("config.txt");
//  DENDRITE_UINT eleOrder = 0;
//  DENDRITE_UINT level = 0;
//  bool mfree = false;
//  if (configFile.good()) {
//    if (!idata.ReadFromFile()) {  /// read from file named "config.txt"
//      throw std::runtime_error("[ERR] Error reading input data, check the config file!");
//    }
//    if (!idata.CheckInputData()) {
//      throw std::runtime_error("[ERR] Problem with input data, check the config file!");
//    }
//    eleOrder = idata.basisFunction;
//    level = idata.mesh_def.refine_lvl;
//    mfree = idata.mfree;
//  } else if (argc < 3) {
//    TALYFEMLIB::PrintStatus("Usage: ", argv[0], " eleOrder level mfree");
//    return -1;
//  } else {
//    eleOrder = std::atoi(argv[1]);
//    level = std::atoi(argv[2]);
//    mfree = std::atoi(argv[3]);
//  }
//
//  TALYFEMLIB::PrintInfo("Total number of processor = ", TALYFEMLIB::GetMPISize());
//  TALYFEMLIB::PrintInfo("size of DendroInt ", sizeof(DendroIntL));
//  TALYFEMLIB::PrintInfo("size of PetscInt ", sizeof(PetscInt));
//
//  TALYFEMLIB::PrintStatus("eleOrder ", eleOrder);
//  TALYFEMLIB::PrintStatus("Level ", level);
//  TALYFEMLIB::PrintStatus("Mfree ", mfree);
//  TALYFEMLIB::PrintStatus("DIM =  ", DIM);
//
//  /// Create DA
//  const char *varname[]{"T"};
//
//  static const DENDRITE_UINT ndof = SSHTNodeData::valueno();
//
//  /// Physical dimensions.
//    DomainInfo cubeDomain, physDomain;
//    cubeDomain.min.fill(0.0);
//    cubeDomain.max.fill(1.0);
////    cubeDomain.max[0] = 0.5;
//    physDomain.min.fill(0.0);
//    physDomain.max.fill(1.0);
////    physDomain.max[0] = 0.5;
//    DomainExtents domainExtents(cubeDomain,physDomain);
//    /// SubDomain Objects
//    DistTREE distTree;
//    SubDomain subDomain(domainExtents);//
//    double center[DIM]{0.5,0.5};
//    VOXEL::Circle circle(center,0.15);
//    subDomain.addObject(circle);
//    std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
//        return (subDomain.functionToRetain(physCoords, physSize));
//    };
//
//
//    /// Create DA
//    DA * octDA = createSubDA(distTree,functionToRetain,level,eleOrder,0.3);



    dendrite_init(argc, argv);

    const char *varname[]{"T"};
    SSHTInputData idata;

    static const DENDRITE_UINT ndof = SSHTNodeData::valueno();


    /// use config.txt first
    std::ifstream configFile("config.txt");
    DENDRITE_UINT eleOrder = 0;
    DENDRITE_UINT level = 0;
    bool mfree = false;
    if (configFile.good()) {
        if (!idata.ReadFromFile()) {  /// read from file named "config.txt"
            throw std::runtime_error("[ERR] Error reading input data, check the config file!");
        }
        if (!idata.CheckInputData()) {
            throw std::runtime_error("[ERR] Problem with input data, check the config file!");
        }
        eleOrder = idata.basisFunction;
        level = idata.refine_lvl;
        mfree = idata.mfree;
    } else {
        TALYFEMLIB::PrintStatus("Fix the config file!");
        return -1;
    }

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




    TALYFEMLIB::PrintInfo("Total number of processor = ", TALYFEMLIB::GetMPISize());
    TALYFEMLIB::PrintInfo("size of DendroInt ", sizeof(DendroIntL));
    TALYFEMLIB::PrintInfo("size of PetscInt ", sizeof(PetscInt));

    TALYFEMLIB::PrintStatus("eleOrder ", eleOrder);
    TALYFEMLIB::PrintStatus("Level ", level);
    TALYFEMLIB::PrintStatus("DIM =  ", DIM);

    // ********************************************************************************************************************
    // Create STL file carved out domain using the same domainExtents

#if(DIM==2)
    {
        GEOMETRY::MSH msh(idata.mshFileName.c_str());
        Point<DIM> temp(0.0);
        GEOMETRY::Geometry *mshGeo = new GEOMETRY::Geometry(&msh, temp, idata.mshRetainSide);
        subDomain.addObject(mshGeo);

    }
#endif


    GEOMETRY::STL stl(idata.stlFileName.c_str());



#if(DIM==3)
    {
        // Use CreateMesh to create a mesh with GEOMETRY::STL as carve out object
        Point<DIM> temp2(0.0);
        GEOMETRY::Geometry *stlGeo = new GEOMETRY::Geometry(&stl, temp2, idata.stlRetainSide);
//        CreateMesh(idata, stlGeo, level, eleOrder, "stl");
        subDomain.addObject(stlGeo);


    }
#endif

    // Retain Function
    auto functionToRetain = [&](const double *physCoords,double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };
    /// Create DA
//        std::string stll;
    DistTREE distTree;
    DA *octDA = createSubDA(distTree, functionToRetain, level, eleOrder, 0.3);

    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
    performNoChangeRefinement(octDA, distTree, domainExtents, 5, subDomain);
//    IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), stl.c_str(), stl.c_str(), domainExtents);









//
//
//#if(DIM==2)
//    {
//        VOXEL::Circle circle(idata.circleCenter, idata.circleRadius, idata.circleRetainSide);
//
//    }
//
//    {
//        VOXEL::Box box(idata.rectLowerLeft, idata.rectUpperRight, idata.rectRetainSide);
//        CreateMesh(idata,box,level,eleOrder,"rectangle");
//    }
//    {
//        GEOMETRY::MSH msh(idata.mshFileName.c_str());
//        Point<DIM> temp(0.0);
//        GEOMETRY::Geometry *mshGeo = new GEOMETRY::Geometry(&msh, temp, idata.mshRetainSide);
//        CreateMesh(idata,mshGeo,level,eleOrder,"msh");
//        delete mshGeo;
//    }
//#endif
//#if(DIM==3)
//    {
//        // Use CreateMesh to create a mesh with VOXEL::Sphere as carve out object
//        VOXEL::Sphere sphere(idata.sphereCenter, idata.sphereRadius, idata.sphereRetainSide);
//        CreateMesh(idata, sphere, level, eleOrder, "sphere");
//    }
//
//    {
//        // Use CreateMesh to create a mesh with VOXEL::Box as carve out object
//        VOXEL::Box box(idata.boxLowerLeft, idata.boxUpperRight, idata.boxRetainSide);
//        CreateMesh(idata, box, level, eleOrder, "box");
//    }
//
//    {
//        // Use CreateMesh to create a mesh with GEOMETRY::STL as carve out object
//        GEOMETRY::STL stl(idata.stlFileName.c_str());
//        Point<DIM> temp2(0.0);
//        GEOMETRY::Geometry *stlGeo = new GEOMETRY::Geometry(&stl, temp2, idata.stlRetainSide);
//        CreateMesh(idata, stlGeo, level, eleOrder, "stl");
//        delete stlGeo;
//    }
//#endif
//    dendrite_finalize();
////    dendrite_finalize(octDA);
////    dendrite_finalize(octDA2);
////    dendrite_finalize(octDA3);

















  /// Equation setup
  auto sshtEq = new TalyEquation<ADEquation, SSHTNodeData>(octDA, distTree.getTreePartFiltered(), domainExtents);
  LinearSolver *sshtSolver = setLinearSolver(sshtEq, octDA,distTree, ndof, mfree);
  if (configFile.good()) {
    /// apply solver parameter from config.txt
    idata.solverOptionsSSHT.apply_to_petsc_options("-ssht_");
    {
      KSP m_ksp = sshtSolver->ksp();
      KSPSetOptionsPrefix(m_ksp, "ssht_");
      KSPSetFromOptions(m_ksp);
    }
  }

  /// Boundary condition
  sshtSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
    Boundary b;
    static constexpr double eps = 1e-14;
    double x = pos.x();
    double y = pos.y();

    bool on_wall = (fabs(x - physDomain.min[0]) < eps) ||
        (fabs(x - physDomain.max[0]) < eps) ||
        (fabs(y - physDomain.min[1]) < eps) ||
        (fabs(y - physDomain.max[1]) < eps);
#if (DIM >= 3)
    double z = pos.z();
    on_wall = on_wall || (fabs(z - physDomain.min[2]) < eps) || (fabs(z - physDomain.max[2]) < eps);
#endif
#if (DIM == 4)
    double t = pos.t();
    on_wall = on_wall ||  (fabs(t) < eps) || (fabs(t - physDomain.min[3]) < eps);
#endif
    if (on_wall) {
      b.addDirichlet(0, 0);
    }
    return b;
  });

  /// Solve
  sshtSolver->solve();



  /// Print files
  petscVectopvtu(octDA, distTree.getTreePartFiltered(), sshtSolver->getCurrentSolution(), "ssht", varname, domainExtents, false, false, ndof);

  /// L2 error
  const auto analytic_sol = [](const TALYFEMLIB::ZEROPTV &pos, const DENDRITE_UINT dof, const DENDRITE_REAL time) {
#if (DIM == 2)
    return sin(M_PI * pos.x()) * sin(M_PI * pos.y());
#elif (DIM == 3)
    return sin(M_PI * pos.x()) * sin(M_PI * pos.y()) * sin(M_PI * pos.z());
#endif
  };
  VecInfo v(sshtSolver->getCurrentSolution(), 1, 0);
  Analytic sshtAnalytic(octDA, distTree.getTreePartFiltered(), v, analytic_sol, physDomain);
  sshtAnalytic.getL2error();

  /// Save vector file (only for regression test)
  if (idata.dump_vec) {
    petscDumpFilesforRegressionTest(octDA, sshtSolver->getCurrentSolution(), "solution_vec.vec");
  }

  delete sshtEq;
  delete sshtSolver;

  dendrite_finalize(octDA);
//    dendrite_finalize();

}