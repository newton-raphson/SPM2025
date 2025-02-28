//
// Created by maksbh on 11/6/22.
//

#ifndef DENDRITEKT_UTILS_H
#define DENDRITEKT_UTILS_H
#include <oda.h>
#include "Traversal/Refinement.h"
#include "ElementMarker.h"
#include "IO/VTU.h"
#include "BFS.h"
#include <PETSc/IO/petscVTU.h>
void performRefinement(DA *& octDA, DistTREE & distTree, DomainExtents & domainExtents, SubDomain & subDomain){
  std::vector<ot::OCT_FLAGS::Refine> refineFlags(octDA->getLocalElementSz());
  std::fill(refineFlags.begin(), refineFlags.end(), ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
  DistTREE newDistTree, surrDistTree;
  DistTREE::distRemeshSubdomain(distTree, refineFlags, newDistTree, surrDistTree,ot::RemeshPartition::SurrogateOutByIn,0.3);

  DA *newDA = new DA(newDistTree, MPI_COMM_WORLD, octDA->getElementOrder(), 100, 0.3);
  std::swap(newDistTree, distTree);
  std::swap(newDA, octDA);
  delete newDA;
  subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
}


void generateNeighborsOfFalseIntercepted(DA *& octDA, DistTREE & distTree, DomainExtents & domainExtents, SubDomain & subDomain,
                                         std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & elementMarker,
                                         Vec & nodalFalseElement,  bool isAllocated = false){
    using CoordT = typename ot::DA<DIM>::C;
    using ot::RankI;
    OctToPhysical octToPhysical(domainExtents);
    const DENDRITE_UINT nPe = octDA->getNumNodesPerElement();
    const auto &treeNodes = distTree.getTreePartFiltered();
    std::vector<PetscInt> nodeIDs(treeNodes.size() * nPe, -1);

    // Get Global element ID
    {
        const std::vector<RankI> &ghostedGlobalNodeId = octDA->getNodeLocalToGlobalMap();
        const size_t ghostedNodalSz = octDA->getTotalNodalSz();
        const TREENODE *odaCoords = octDA->getTNCoords();
        const bool visitEmpty = false;
        const unsigned int padLevel = 0;
        ot::MatvecBaseIn<DIM, RankI, false> treeLoopIn(ghostedNodalSz,
                                                       1,                // node id is scalar
                                                       octDA->getElementOrder(),
                                                       visitEmpty,
                                                       padLevel,
                                                       odaCoords,
                                                       &(*ghostedGlobalNodeId.cbegin()),
                                                       &(*treeNodes.cbegin()),
                                                       treeNodes.size(),
                                                       *octDA->getTreePartFront(),
                                                       *octDA->getTreePartBack());

        int eleCounter = 0;
        while (!treeLoopIn.isFinished()) {
            const ot::TreeNode<CoordT, DIM> subtree = treeLoopIn.getCurrentSubtree();
            const auto subtreeInfo = treeLoopIn.subtreeInfo();
            if (treeLoopIn.isPre() && subtreeInfo.isLeaf()) {
                const RankI *nodeIdsFlat = subtreeInfo.readNodeValsIn();
                const auto &octCoords = subtreeInfo.getNodeCoords();
                const std::vector<bool> &nodeNonhangingIn = subtreeInfo.readNodeNonhangingIn();
                for (int i = 0; i < nPe; i++) {
                    //if (nodeNonhangingIn[i]) {
                    nodeIDs[eleCounter * nPe + i] = nodeIdsFlat[i];
                    //}
                }
                eleCounter++;
            }
            treeLoopIn.step();
        }
    }
    // Get Global node ID ends

    // Transfer cell information to nodes
    if(!isAllocated) {
        octDA->petscCreateVector(nodalFalseElement, false, false, 1);
    }
    PetscScalar  sum;
    VecSet(nodalFalseElement,0);
    for(int i = 0; i < treeNodes.size();i++){
        if(elementMarker[i].test(ElementMarker::SBM_FALSE_INTERCEPTED)){
            for(int n = 0; n < nPe; n++){
                VecSetValue(nodalFalseElement,nodeIDs[nPe*i + n],NodeMarker::SBM_FALSE_INTERCEPTED_NODES,INSERT_VALUES);
            }
        }
    }
    VecAssemblyBegin(nodalFalseElement);
    VecAssemblyEnd(nodalFalseElement);

    VecSum(nodalFalseElement,&sum);

    // Do one iteration of BFS to find the neighbors of false intercepted
    {
        BFS bfs(octDA, distTree.getTreePartFiltered(), elementMarker,domainExtents,VecInfo(nodalFalseElement, 1, 0), &subDomain);

    }
}
void generateNewMarkers(DA * octDA, DistTREE & distTree, DomainExtents & domainExtents, SubDomain & subDomain, std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & elementMarker,Vec &nodalFalseElement){
    SBMMarker marker(octDA,distTree.getTreePartFiltered(),domainExtents,elementMarker,&subDomain);

    const char *varname[] = {"marker"};
    std::vector<double> printMarker(elementMarker.size());
    for(int i = 0; i < elementMarker.size();i++){
        printMarker[i] = (double )(elementMarker[i].to_ulong());
    }

    IO::writeVecTopVtu(octDA,distTree.getTreePartFiltered(),printMarker.data(),"Marker","marker",varname,domainExtents,true);
    octDA->petscCreateVector(nodalFalseElement, false, false, 1);

    generateNeighborsOfFalseIntercepted(octDA, distTree, domainExtents,  subDomain,elementMarker,nodalFalseElement, true);
    for(int i = 0; i < elementMarker.size();i++){
        printMarker[i] = (double )(elementMarker[i].to_ulong());
    }
    IO::writeVecTopVtu(octDA,distTree.getTreePartFiltered(),printMarker.data(),"MarkerFinal","marker",varname,domainExtents,true);
    PETSc::petscVectopvtu(octDA,distTree.getTreePartFiltered(),nodalFalseElement,"MarkerNeighbors","Nodes",varname,domainExtents,
                   false);

}

#endif //DENDRITEKT_UTILS_H
