//
// Created by maksbh on 3/26/23.
//

#ifndef DENDRITEKT_REFINEMENTUTILS_H
#define DENDRITEKT_REFINEMENTUTILS_H

#include "Refinement/BoundaryRefine.h"

static void getnewDABoundaryRefine(DA * & octDA,DistTREE & distTree,DomainExtents & domainExtents, int numRefineTimes){
    for(int i = 0 ; i < numRefineTimes; i++){
        BoundaryRefine boundaryRefine(octDA,distTree.getTreePartFiltered(),domainExtents);
        DA * newDA = boundaryRefine.getRefineSubDA(distTree,0.1,RefinementStrategy::FULL_DOMAIN,ot::RemeshPartition::SurrogateOutByIn);
        std::swap(octDA,newDA);
        delete newDA;
    }

}

static void detectHangingNodes(DA * octDA, DistTREE & distTree){
    const size_t sz = octDA->getTotalNodalSz();
    auto partFront = octDA->getTreePartFront();
    auto partBack = octDA->getTreePartBack();
    const auto tnCoords = octDA->getTNCoords();
    using ot::RankI;
    using CoordT = typename ot::DA<DIM>::C;
    const DENDRITE_UINT nPe = octDA->getNumNodesPerElement();
    const auto &treeNodes = distTree.getTreePartFiltered();
    std::vector<PetscInt> nodeIDs(treeNodes.size() * nPe, -1);
    struct NodeCoords{
        double nodeCoords[DIM];
        NodeCoords(const double * coords){
            std::memcpy(nodeCoords,coords,DIM*sizeof(double ));
        }
        bool operator==(const NodeCoords & a) const {
            double sum = 0.0;
            for(int d = 0; d < DIM; d++){
                sum += std::abs(a.nodeCoords[d] - this->nodeCoords[d]);
            }
            if(FEQUALS(sum,0.0)) {
                return true;
            }
            else {
                return false;
            }
        }
        bool operator < (const NodeCoords & a) const {
            for(int d = 0; d < DIM; d++){
                if( this->nodeCoords[d]  < a.nodeCoords[d] ){
                    return true;
                }
                else if(this->nodeCoords[d]  > a.nodeCoords[d]){
                    return false;
                }
            }
        }

        bool operator > (const NodeCoords & a) const {
            return !(*this < a);
        }
    };
    std::vector<NodeCoords> set_coords;
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
                    if (nodeNonhangingIn[i]) {

                    }
                    else{
                        NodeCoords coords(&octCoords[i*DIM]);

                        set_coords.push_back(coords);
                    }
                }
                eleCounter++;
            }
            treeLoopIn.step();
        }
    }
    std::sort(set_coords.begin(),set_coords.end());
    auto it = std::unique(set_coords.begin(),set_coords.end());
    set_coords.erase(it,set_coords.end());
    for(auto coord:set_coords ){
        for(int d = 0; d < DIM-1; d++){
            std::cout << coord.nodeCoords[d] << ",";
        }
        std::cout << coord.nodeCoords[DIM -1 ] ;
        std::cout << "\n";

    }
}
#endif //DENDRITEKT_REFINEMENTUTILS_H
