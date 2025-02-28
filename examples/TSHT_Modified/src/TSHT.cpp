//
// Created by maksbh on 5/20/20.
//

#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <HTEquation.h>
#include <HTNodeData.h>
#include <HTInputData.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <PETSc/IO/petscVTU.h>
#include <Traversal/Analytic.h>
#include <unistd.h>
#include "SDARefine.h"
#include "SDACoaren.h"

using namespace PETSc;




void performNoChangeRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                               DENDRITE_UINT maxLevel, SubDomain & subDomain){
    // Use this function to perform no refinement but to remove the physical boundary octants in case of
    // retain inside case
    SDARefine refine(octDA, distTree.getTreePartFiltered(),domainInfo,0);
    DA * newDA = refine.getForceRefineSubDA(distTree,0.01);
    std::swap(octDA, newDA);
    delete newDA;
    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);
}

void performRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                       SubDomain & subDomain, const int blevel){
    while(true) {
        SDARefine refine(octDA, distTree.getTreePartFiltered(), domainInfo,blevel);

        DA *newDA = refine.getRefineSubDA(distTree, 0.01,RefinementStrategy::FULL_DOMAIN,ot::RemeshPartition::SurrogateInByOut);
        if(newDA == nullptr){
            break;
        }
        std::swap(octDA, newDA);
        delete newDA;
        subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);
    }
}

void performCoarsening(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                       SubDomain & subDomain, const int blevel){
    while(true) {
        SDARefine refine(octDA, distTree.getTreePartFiltered(), domainInfo,blevel);

        DA *newDA = refine.getRefineSubDA(distTree, 0.01,RefinementStrategy::FULL_DOMAIN,ot::RemeshPartition::SurrogateOutByIn);
        if(newDA == nullptr){
            break;
        }
        std::swap(octDA, newDA);
        delete newDA;
        subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);
    }
}



int main(int argc, char *argv[]) {

  dendrite_init(argc, argv);
  int rank = TALYFEMLIB::GetMPIRank();
  /// read parameters from config.txt
  HTInputData idata;

  /// use config.txt first
  std::ifstream configFile("config.txt");
  DENDRITE_UINT eleOrder = 0;
  DENDRITE_UINT level = 0;
  bool mfree = false;

  std::vector<DENDRITE_REAL> dt;

  std::vector<DENDRITE_REAL> totalT;


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
    dt = idata.dt;

    totalT = idata.totalT;
  } else if (argc < 3) {
    TALYFEMLIB::PrintStatus("Usage: ", argv[0], " eleOrder level mfree");
    return -1;
  } else {
    eleOrder = std::atoi(argv[1]);
    level = std::atoi(argv[2]);
    mfree = std::atoi(argv[3]);
  }

  TALYFEMLIB::PrintInfo("Total number of processor = ", TALYFEMLIB::GetMPISize());
  TALYFEMLIB::PrintInfo("size of DendroInt ", sizeof(DendroIntL));
  TALYFEMLIB::PrintInfo("size of PetscInt ", sizeof(PetscInt));

  TALYFEMLIB::PrintStatus("eleOrder ", eleOrder);
  TALYFEMLIB::PrintStatus("Level ", level);
  TALYFEMLIB::PrintStatus("Mfree ", mfree);
  TALYFEMLIB::PrintStatus("DIM =  ", DIM);




    /// Create DA
    const char *varname[]{"T"};

    static const DENDRITE_UINT ndof = 1;

    /// Physical dimensions.
    DomainInfo cubeDomain, physDomain;
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
        cubeDomain.min[dim] = idata.cubeDomainMin[dim];
        cubeDomain.max[dim] = idata.cubeDomainMax[dim];
        physDomain.min[dim] = idata.physDomainMin[dim];
        physDomain.max[dim] = idata.physDomainMax[dim];
    }


    DomainExtents domainExtents(cubeDomain,physDomain);
    /// SubDomain Objects
    DistTREE distTree;
    SubDomain subDomain(domainExtents);

#if(DIM==3)
    GEOMETRY::STL stl("spacex.stl");
    Point<DIM> temp2(0.0);
    GEOMETRY::Geometry *stlGeo = new GEOMETRY::Geometry(&stl, temp2, RetainSide::IN);
    subDomain.addObject(stlGeo);
#elif(DIM==2)
    double center[DIM]{0.5,0.5};
    VOXEL::Circle circle(center,0.15);
    subDomain.addObject(circle);
#endif

    std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };

    /// Create DA
    DA * octDA = createSubDA(distTree,functionToRetain,level,eleOrder,0.3);



////
//    performNoChangeRefinement(octDA, distTree, domainExtents, 8, subDomain);
//    performRefinement(octDA, distTree, domainExtents, subDomain, 8);
//    performRefinement(octDA, distTree, domainExtents, subDomain, 8+1);
////    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
////    performCoarsening(octDA, distTree, domainExtents, subDomain, 8);
//    performCoarsening(octDA, distTree, domainExtents, subDomain, 8);









  TimeInfo ti(0.0, dt, totalT);


  auto htEq = new TalyEquation<HTEquation, HTNodeData>(octDA, distTree.getTreePartFiltered(), domainExtents, ndof,&ti, false,nullptr, &idata);

  LinearSolver *htSolver = setLinearSolver(htEq, octDA, distTree, ndof, mfree);
  /// octDA, treePart, domainExtents, ndof, &ti, false, nullptr, &idata
  if (configFile.good()) {
    /// apply solver parameter from config.txt
    idata.solverOptionsHT.apply_to_petsc_options("-ht_");
    {
      KSP m_ksp = htSolver->ksp();
      KSPSetOptionsPrefix(m_ksp, "ht_");
      KSPSetFromOptions(m_ksp);
    }
  }

  htSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
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
    if (on_wall) {
      b.addDirichlet(0, 0);
    }
    return b;
  });



  Vec prev_solution, prev_prev_solution;
  octDA->petscCreateVector(prev_solution, false, false, ndof);
  octDA->petscCreateVector(prev_prev_solution, false, false, ndof);
  std::function<void(const double *, double *)> initial_condition = [](const double *x, double *var) {
#if (DIM == 2)
    var[0] = sin(M_PI * x[0]) * sin(M_PI * x[1]);
#elif(DIM == 3)
    var[0] = sin(M_PI * x[0]) * sin(M_PI * x[1]) * sin(M_PI * x[2]);
#endif
  };
  octDA->petscSetVectorByFunction(prev_prev_solution, initial_condition, false, false, ndof);
  VecCopy(prev_prev_solution, prev_solution);
  htEq->setVectors({VecInfo(prev_solution, ndof, U_PRE),
                    VecInfo(prev_prev_solution, ndof, U_PRE_PRE)},
                   SYNC_TYPE::ALL);


  if (configFile.good()) {
    dt = idata.dt;
    totalT = idata.totalT;
  }
  TALYFEMLIB::PrintStatus("Pid = ", getpid());
  TALYFEMLIB::PrintStatus("Starting solve");

  while (ti.getCurrentTime() < ti.getEndTime() - 1e-15) {

    TALYFEMLIB::PrintStatus("Time = ", ti.getCurrentTime());
    htSolver->solve();
    VecCopy(prev_solution, prev_prev_solution);
    VecCopy(htSolver->getCurrentSolution(), prev_solution);
    ti.increment();
    /// Save solution when required
    if ((ti.getCurrentTime() >= (idata.OutputStartTime - 1e-16))
        && (ti.getTimeStepNumber() % idata.OutputInterval == 0)) {
      char fname[PATH_MAX];
      snprintf(fname, sizeof(fname), "%s_%03d", "ht", ti.getTimeStepNumber());
      petscVectopvtu(octDA, distTree.getTreePartFiltered(),htSolver->getCurrentSolution(), fname, varname, domainExtents, false, false, ndof);
      IO::writeBoundaryElements(octDA,distTree.getTreePartFiltered(),"boundary","boundary",domainExtents);
    }
  }

  const auto analytic_sol = [&](TALYFEMLIB::ZEROPTV pos, int dof, double currentTime) {
    if (configFile.good() && idata.forcing) {
#if(DIM == 2)
      return sin(M_PI * pos.x()) * sin(M_PI * pos.y()) * cos(M_PI * 2 * ti.getCurrentTime());
#elif (DIM == 3)
      return sin(M_PI * pos.x()) * sin(M_PI * pos.y()) * sin(M_PI * pos.z()) * cos(M_PI * 2 * ti.getCurrentTime());
#endif
    } else {
#if(DIM == 2)
      return sin(M_PI * pos.x()) * sin(M_PI * pos.y()) * exp(-currentTime);
#elif (DIM == 3)
      return sin(M_PI * pos.x()) * sin(M_PI * pos.y()) * sin(M_PI * pos.z()) * exp(-currentTime);
#endif
    }
  };
  /// Save vector file (only for regression test)
  if (idata.dump_vec) {
    petscDumpFilesforRegressionTest(octDA, htSolver->getCurrentSolution(), "solution_vec.vec");
  }

  VecInfo v(prev_solution, 1, 0);
  Analytic htAnalytic(octDA, distTree.getTreePartFiltered(), v, analytic_sol, domainExtents, ti.getCurrentTime());
  htAnalytic.getL2error();

  /** Cleanup**/
  VecDestroy(&prev_solution);
  VecDestroy(&prev_prev_solution);
  delete htEq;
  delete htSolver;
  dendrite_finalize(octDA);

}