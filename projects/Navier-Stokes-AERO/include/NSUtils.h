//
// Created by maksbh on 6/13/20.
//

#ifndef DENDRITEKT_NSUTILS_H
#define DENDRITEKT_NSUTILS_H

#include <DataTypes.h>
#include <NSRefine.h>
#include <SDARefine.h>
#include <PETSc/Solver/NonLinearSolver.h>
#include <NSInputData.h>
#include <NSNodeData.h>
#include <Traversal/Analytic.h>
using namespace PETSc;

//void performRefinement(DA *&octDA, const DomainExtents &domainInfo, std::vector<TREENODE> &treeNode, const NSInputData &inputData) {
//  while (true) {
//    NSRefine refine(octDA, treeNode,domainInfo, treeNode.size(), inputData.meshDef.refineLevelBoundary);
//    DA *newDA = refine.getRefineDA(treeNode);
//    if (newDA == NULL) {
//      break;
//    }
//    std::swap(newDA, octDA);
//    delete newDA;
//  }
//}
//
//
//void performRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
//                       DENDRITE_UINT maxLevel, SubDomain & subDomain, NSInputData *inputData){
//    while(true) {
//        SDARefine refine(octDA, distTree.getTreePartFiltered(), domainInfo, maxLevel, inputData);
//
//        DA *newDA = refine.getRefineSubDA(distTree, 0.01,RefinementStrategy::FULL_DOMAIN,ot::RemeshPartition::SurrogateInByOut);
//        if(newDA == nullptr){
//            break;
//        }
//        std::swap(octDA, newDA);
//        delete newDA;
////        subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);
//    }
//
//}
//
//
//void performNoChangeRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
//                               DENDRITE_UINT maxLevel, SubDomain & subDomain, NSInputData *inputData){
//    // Use this function to perform no refinement but to remove the physical boundary octants in case of
//    // retain inside case
////    SDARefine refine(octDA, distTree.getTreePartFiltered(),domainInfo,maxLevel);
//    SDARefine refine(octDA, distTree.getTreePartFiltered(), domainInfo, maxLevel, inputData);
//    DA * newDA = refine.getForceRefineSubDA(distTree,0.01);
//    std::swap(octDA, newDA);
//    delete newDA;
//    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);
//}
//





//void performRefinementSubDA(DA *&octDA, DomainExtents &domainExtents, DistTREE &dTree,
//                            NSHTInputData &inputData, SubDomain *subdomain,const IMGA *imga, const my_kd_tree_t *kd_tree)
//{
//    int no_refine = 0;
//    while (true)
//    {
//
//        SubDomainBoundary subDomainBoundary(subdomain, octDA, domainExtents);
//
//        NSHTRefine refine(octDA, dTree.getTreePartFiltered(), domainExtents, &inputData, &subDomainBoundary,imga,kd_tree);
//        DA *newDA = refine.getRefineSubDA(dTree);
//        //DA *newDA = refine.getRefineSubDA(dTree, 0.03, RefinementStrategy::SUBDOMAIN, ot::RemeshPartition::SurrogateOutByIn);
//        //DA *newDA = refine.getRefineSubDA(dTree, 0.03, RefinementStrategy::FULL_DOMAIN, ot::RemeshPartition::SurrogateInByOut);
//        if (newDA == NULL)
//        {
//            newDA = refine.getForceRefineSubDA(dTree);
//            std::swap(newDA, octDA);
//            break;
//        }
//
//        std::swap(newDA, octDA);
//
//        delete newDA;
//
//        subdomain->finalize(octDA, dTree.getTreePartFiltered(), domainExtents);
//        TALYFEMLIB::PrintStatus("Refinement ", no_refine, ", mesh count (node) = ", octDA->getGlobalNodeSz());
//        no_refine++;
//    }
//}





static double AnalyticalSolution(const TALYFEMLIB::ZEROPTV &pt, int dof, const DENDRITE_REAL time) {
#if(DIM == 3)
  std::cout << "Analytical solution not supported \n";
    return 0;
#else

  if (dof == NSNodeData::VEL_X) {
    return (sin(M_PI * pt.x()) * cos(M_PI * pt.y()) * sin(2 * M_PI * (time)) + 2);
  }
  if (dof == NSNodeData::VEL_Y) {
    return (-cos(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(2 * M_PI * (time)) + 2);
  }
  if (dof == NSNodeData::PRESSURE) {
    return (sin(M_PI * pt.x()) * sin(M_PI * pt.y()) * cos(2 * M_PI * (time)) + 2);
  }
#endif
};
void setInitialCondition(DA *octDA, NSInputData &inputData, Vec &Solution) {
  std::function<void(const double *, double *)> initial_condition = [&](const double *x, double *var) {
    if (not(inputData.ifMMS)) {
      var[NSNodeData::VEL_X] = 0;
      var[NSNodeData::VEL_Y] = 0;
#if(DIM == 3)
      var[NSNodeData::VEL_Z] = 0;
#endif
      var[NSNodeData::PRESSURE] = 0;
    } else {
      var[NSNodeData::VEL_X] = AnalyticalSolution(TALYFEMLIB::ZEROPTV{x[0], x[1], 0.0}, NSNodeData::VEL_X, 0.0);
      var[NSNodeData::VEL_Y] = AnalyticalSolution(TALYFEMLIB::ZEROPTV{x[0], x[1], 0.0}, NSNodeData::VEL_Y, 0.0);
      var[NSNodeData::PRESSURE] = AnalyticalSolution(TALYFEMLIB::ZEROPTV{x[0], x[1], 0.0}, NSNodeData::PRESSURE, 0.0);
    }
  };
  octDA->petscSetVectorByFunction(Solution, initial_condition, false, false, NSNodeData::NS_DOF);
}

#endif //DENDRITEKT_NSUTILS_H
