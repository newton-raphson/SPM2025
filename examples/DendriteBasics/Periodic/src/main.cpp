//
// Created by maksbh on 8/12/21.
//

#include <DendriteUtils.h>
#include <SubDA/SubDomain.h>
#include <PETSc/IO/petscVTU.h>
#include <PRefine.h>
#include <ADEquation.h>
#include <TimeInfo.h>
#include <ADNodeData.h>
#include <TalyEquation.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>

using namespace PETSc;
std::array<DENDRITE_UINT,DIM> PERIODS;
int main(int argc, char *argv[]) {
  dendrite_init(argc,argv);
////
  DomainInfo cubeDomain,physDomain;
  cubeDomain.min.fill(0.0);
  cubeDomain.max.fill(1.0);

  physDomain.min.fill(0.0);
  physDomain.max.fill(1.0);
  physDomain.max[1] = 0.5;

  DomainExtents domainExtents(cubeDomain,physDomain);


  PERIODS[0] =  (1u << m_uiMaxDepth);
  PERIODS[1] =  (1u << m_uiMaxDepth)*0.5;
  periodic::PCoord<DENDRITE_UINT , DIM>::periods(PERIODS);

  SubDomain subDomain(domainExtents);
  std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
    return (subDomain.functionToRetain(physCoords, physSize));
  };
  DistTREE dtree;
  unsigned int level = 3;
  DA * octDA = createSubDA(dtree,functionToRetain,level,1,0.1);
  std::cout << octDA->getTotalNodalSz() << " " <<  octDA->getLocalElementSz() << "\n";
  IO::writeBoundaryElements(octDA,dtree.getTreePartFiltered(),"bnd","bnd",domainExtents);
 /* PRefine refine(octDA,dtree.getTreePartFiltered(),domainExtents);
  DA * newDA = refine.getRefineSubDA(dtree);
  std::swap(newDA,octDA);
  delete newDA;
  std::cout << octDA->getTotalNodalSz() << " " <<  octDA->getLocalElementSz() << "\n";
  IO::writeBoundaryElements(octDA,dtree.getTreePartFiltered(),"newbnd","bnd",domainExtents);*/
  std::vector<double> dt = {0.1};
  std::vector<double> totalT = {1};

  std::function<void(const double *, double *)> initial_condition = [](const double *x, double *var) {
    var[0] = exp(-((x[0] - 0.5)*(x[0] - 0.5)/0.01));
  };
  Vec prevSolution;
  octDA->petscCreateVector(prevSolution,false,false,ADNODEDATA_MAX);
  octDA->petscSetVectorByFunction(prevSolution,initial_condition,false,false,ADNODEDATA_MAX);

  TimeInfo ti(0.0, dt, totalT);

  auto adEq =
    new TalyEquation<ADEquation, ADNodeData>(octDA, dtree.getTreePartFiltered(),
                                             subDomain.domainExtents(), ADNODEDATA_MAX, &ti,false,
                                             nullptr/*, &integrandsGenForm*/);

  LinearSolver *adSolver = setLinearSolver(adEq, octDA, ADNODEDATA_MAX, false,true);
  VecCopy(prevSolution,adSolver->getCurrentSolution());
  static const char *varname[]{"C"};
  adEq->setVectors({VecInfo(prevSolution, 1, 0)},
                   SYNC_TYPE::ALL);
  petscVectopvtu(octDA,dtree.getTreePartFiltered(),prevSolution,"init","init",varname,domainExtents,false,false,ADNODEDATA_MAX);
  while(ti.getCurrentTime() < ti.getEndTime()){
    char fname[PATH_MAX];
    char folder[PATH_MAX];
    snprintf(folder, sizeof(folder), "results_%06d", ti.getTimeStepNumber());
    const std::string prefix = "ad";
    snprintf(fname, sizeof(fname), "%s_%06d", prefix.c_str(), ti.getTimeStepNumber());
    static const char * varName[]{"C"};
    PETSc::petscVectopvtu(octDA, dtree.getTreePartFiltered(), prevSolution, folder, fname, varName,
                          subDomain.domainExtents(), false, false, 1);
    ti.print();
    adSolver->solve();
    VecCopy(adSolver->getCurrentSolution(),prevSolution);
    ti.increment();

  }
  petscVectopvtu(octDA,dtree.getTreePartFiltered(),prevSolution,"final","init",varname,domainExtents,false,false,ADNODEDATA_MAX);

  dendrite_finalize(octDA);
}