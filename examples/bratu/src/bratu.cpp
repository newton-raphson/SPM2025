//
// Created by maksbh on 9/18/19.
//

#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <BTNodeData.h>
#include <BTEquation.h>
#include <BTInputData.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <PETSc/IO/petscVTU.h>
#include <Traversal/Analytic.h>

using namespace PETSc;
int main(int argc, char *argv[]) {

  dendrite_init(argc, argv);
  int rank = TALYFEMLIB::GetMPIRank();
  DENDRITE_UINT eleOrder = 1;
  DENDRITE_UINT level = 4;

  bool mfree = false;

  DomainInfo physDomain;
  physDomain.min.fill(0.0);
  physDomain.max.fill(0.5);
  DomainExtents domainExtents(physDomain);


  TALYFEMLIB::PrintInfo("Total number of processor = ", TALYFEMLIB::GetMPISize());
  TALYFEMLIB::PrintInfo("size of DendroInt ", sizeof(DendroIntL));
  TALYFEMLIB::PrintInfo("size of PetscInt ", sizeof(PetscInt));

  TALYFEMLIB::PrintStatus("eleOrder ", eleOrder);
  TALYFEMLIB::PrintStatus("Level ", level);
  TALYFEMLIB::PrintStatus("Mfree ", mfree);
  TALYFEMLIB::PrintStatus("DIM =  ", DIM);

  /// create DA
  SubDomain subDomain (domainExtents);
  std::function<ibm::Partition(const double *, double)> functionToRetain = [&](const double *octCoords, double scale) {
    return (subDomain.functionToRetain (octCoords, scale));
  };
  DistTREE  dTree;
  DA * octDA = createSubDA (dTree, functionToRetain, level, eleOrder);


  /// Problem size
  DENDRITE_UINT ndof = 1;


  DomainExtents domain(physDomain);
  /// Equation and solver setup
  auto btEq = new TalyEquation<BTEquation, BTNodeData>(octDA,dTree.getTreePartFiltered(), domain, ndof);
  NonlinearSolver *btSolver =setNonLinearSolver(btEq, octDA,dTree, ndof, mfree);


  const auto analytic_sol = [&](TALYFEMLIB::ZEROPTV pos, int dof, double currentTime) {
    return -log(1 + cos(M_PI * (pos.x())));
  };

  /// For non linear solve : The vector you need and the corresponding dof
  btEq->setVectors({VecInfo(PLACEHOLDER_GUESS, ndof, 0)});

  /// BC setup
  btSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
    Boundary b;
    static constexpr double eps = 1e-14;
    if ((fabs(pos.x() - physDomain.min[0]) < eps) ||
        (fabs(pos.x() - physDomain.max[0]) < eps)) {

      b.addDirichlet(0, analytic_sol(pos, 0, 0));
    }
    return b;
  });
  const char *varname[]{"bt"};

  /// Solve
  btSolver->solve();

  /// Print Files
  petscVectopvtu(octDA, dTree.getTreePartFiltered(),btSolver->getCurrentSolution(), "bratu", varname, domain, false, false, ndof);

  VecInfo v(btSolver->getCurrentSolution(), 1, 0);
  Analytic btAnalytic(octDA, dTree.getTreePartFiltered(), v,  analytic_sol, domain, 0);
  btAnalytic.getL2error();

  dendrite_finalize(octDA);

}