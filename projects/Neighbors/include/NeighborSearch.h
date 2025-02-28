//
// Created by mehdi on 8/7/23.
//


#include "OctToPhysical.h"
#include <talyfem/grid/femelm.h>
#include "DataTypes.h"
#include "mathUtils.h"
#include "TalyDendroSync.h"
#include "Boundary/SubDomainBoundary.h"
#include "DataTypes.h"
#include "talyfem/grid/elem_types/elem2dbox.h"
#include "talyfem/grid/elem_types/elem2dtriangle.h"
#include "sfcTreeLoop_matvec.h"
#include "talyfem/utils/timers.h"
#include <mpi.h>

#ifndef DENDRITEKT_NEIGHBORSEARCH_H
#define DENDRITEKT_NEIGHBORSEARCH_H

#endif //DENDRITEKT_NEIGHBORSEARCH_H


class NodePrint {


    static constexpr unsigned int numBits = intPow(3,DIM);
    void initiateOrder1();



    std::bitset<numBits> bits_;
    static constexpr int MAX_COORDS_PER_ELEMENT = 4;
    static constexpr int QUAD_COORDS_PER_ELEMENT = 4;
    static constexpr int TRI_COORDS_PER_ELEMENT = 3;
    static constexpr int MAX_COORDS_PER_FACE = 3;
    static constexpr int MIDDLE_NODE = 4;
    static constexpr int NUM_FACE = 4;


protected:
    enum GMSH_GRID:int{
        GMSH_TRIANGLE = 2,
        GMSH_QUAD = 3,
        GMSH_LINE = 1
    };

    std::vector<std::array<double,DIM>> m_vertexPosition;
    std::vector<std::tuple<GMSH_GRID,std::vector<PetscInt>>> gridElements;
    std::vector<std::tuple<GMSH_GRID,int,std::vector<PetscInt>>> boundary_elements;


    DA * m_octDA;
    OctToPhysical m_octToPhysical;
    const std::vector<TREENODE> & m_treePart;
    DENDRITE_REAL * target_coords;
    DENDRITE_REAL * neigh_coords;
    SubDomainBoundary * m_subDomainBoundary;


    TALYFEMLIB::ELEM *m_elem;
    TRAVERSAL_TYPE m_traversalType;
    DENDRITE_UINT m_level_target;
    DENDRITE_UINT m_level_neigh;
    bool m_BoundaryOctant = false;
    TALYFEMLIB::GRID m_grid;
    int m_childNum;

    DENDRITE_UINT nPe_;


public:
    NodePrint(DA * octDA, const std::vector<TREENODE> & treePart, const DomainExtents & domainExtents, SubDomainBoundary * subDomainBoundary);
    void generate();
    void checkNeighbors(const ot::RankI *target_nodeIdsFlat, const ot::RankI *neigh_nodeIdsFlat, DENDRITE_UINT m_level_target, DENDRITE_UINT m_level_neigh,
                          bool &EastRefine, bool &WestRefine, bool &NorthRefine, bool &SouthRefine, bool &FrontRefine, bool &BackRefine);
    void checkNeighborsCoord(const ot::RankI *target_nodeIdsFlat, const ot::RankI *neigh_nodeIdsFlat, DENDRITE_UINT m_level_target, DENDRITE_UINT m_level_neigh);
    bool CoordCheck(int target_idx, int neigh_idx);

};


NodePrint::NodePrint(DA * octDA, const std::vector<TREENODE> & treePart, const DomainExtents & domainExtents, SubDomainBoundary * subDomainBoundary)
        :m_octDA(octDA), m_octToPhysical(domainExtents), m_treePart(treePart),m_subDomainBoundary(subDomainBoundary){
    target_coords = new DENDRITE_REAL[intPow(3, DIM) * DIM];
    neigh_coords = new DENDRITE_REAL[intPow(3, DIM) * DIM];
    nPe_ = m_octDA->getNumNodesPerElement();
}


void NodePrint::generate() {
    using CoordT = typename ot::DA<DIM>::C;
    using ot::RankI;

#ifdef BUILD_WITH_PETSC
    using ScalarT = PetscScalar;
    using IndexT = PetscInt;
#else
    using ScalarT = DendroScalar;
    using IndexT = long long unsigned;
#endif
    const size_t ghostedNodalSz = m_octDA->getTotalNodalSz();
    const ot::TreeNode<CoordT, DIM> *odaCoords = m_octDA->getTNCoords();
    const std::vector<RankI> &ghostedGlobalNodeId = m_octDA->getNodeLocalToGlobalMap();
    std::vector<std::tuple<std::array<int,2>,unsigned int>> boundaryFaceID;



    ot::MatvecBase<DIM, RankI> treeloop(ghostedNodalSz, 1, m_octDA->getElementOrder(), odaCoords, &(*ghostedGlobalNodeId.cbegin()),
                                        &(*this->m_treePart.cbegin()),
                                        this->m_treePart.size(),
                                        *m_octDA->getTreePartFront(),
                                        *m_octDA->getTreePartBack());










//    MPITimer timer("loop object");
//    timer.Start();
//
//
//
//
//    timer.Stop();
//    timer.PrintTotalTimeSeconds();


    int counter = 0;
    while (!treeloop.isFinished()) {
        if (treeloop.isPre() && treeloop.subtreeInfo().isLeaf()) {

            ot::MatvecBase<DIM, RankI> treeloop2(ghostedNodalSz, 1, m_octDA->getElementOrder(), odaCoords,
                                                &(*ghostedGlobalNodeId.cbegin()),
                                                 &(*this->m_treePart.cbegin()),
                                                 this->m_treePart.size(),
                                                 *m_octDA->getTreePartFront(),
                                                 *m_octDA->getTreePartBack());



            const double *target_nodeCoordsFlat = treeloop.subtreeInfo().getNodeCoords();
            const RankI *target_nodeIdsFlat = treeloop.subtreeInfo().readNodeValsIn();
            const bool isBoundaryOctant = treeloop.subtreeInfo().isElementBoundary();
            bits_ = treeloop.subtreeInfo().getLeafBitsetInfo();
            m_level_target = treeloop.getCurrentSubtree().getLevel();
            std::memcpy(target_coords, target_nodeCoordsFlat, sizeof(DENDRITE_REAL) * nPe_ * DIM);





            std::cout << "ID:" << counter <<  "  is Boundary: " << isBoundaryOctant << std::endl;


            MPITimer timer("loop object");
            timer.Start();

            bool EastRefine = true;
            bool WestRefine = true;
            bool NorthRefine = true;
            bool SouthRefine = true;
            bool FrontRefine = true;
            bool BackRefine = true;
            while (!treeloop2.isFinished()) {
                if (treeloop2.isPre() && treeloop2.subtreeInfo().isLeaf()) {
                    const double *neigh_nodeCoordsFlat = treeloop2.subtreeInfo().getNodeCoords();
                    const RankI *neigh_nodeIdsFlat = treeloop2.subtreeInfo().readNodeValsIn();
                    m_level_neigh = treeloop2.getCurrentSubtree().getLevel();

                    std::memcpy(neigh_coords, neigh_nodeCoordsFlat, sizeof(DENDRITE_REAL) * nPe_ * DIM);


                    checkNeighbors(target_nodeIdsFlat, neigh_nodeIdsFlat, m_level_target, m_level_neigh, EastRefine, WestRefine, NorthRefine,  SouthRefine, FrontRefine, BackRefine);






                    treeloop2.next();
                } else
                    treeloop2.step();
            }

            timer.Stop();
//            timer.PrintTotalTimeSeconds();

            counter ++;



            treeloop.next();
        } else
            treeloop.step();
    }
}


bool NodePrint::CoordCheck(int target_idx, int neigh_idx){
    if (target_coords[target_idx * DIM] == neigh_coords[neigh_idx * DIM] &&
    target_coords[target_idx * DIM + 1] == neigh_coords[neigh_idx * DIM + 1]
#if(DIM==3)
        && target_coords[target_idx * DIM + 2] == neigh_coords[neigh_idx * DIM + 2]
#endif
            ){
        return true;
    }
    else{
        return false;
    }



}




void NodePrint::checkNeighbors(const ot::RankI *target_nodeIdsFlat, const ot::RankI *neigh_nodeIdsFlat, DENDRITE_UINT m_level_target, DENDRITE_UINT m_level_neigh,
                               bool &EastRefine, bool &WestRefine, bool &NorthRefine, bool &SouthRefine, bool &FrontRefine, bool &BackRefine){


    /// East Check
    if (CoordCheck(1, 0) && EastRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "east" << std::endl;
            EastRefine = false;

        }
        else {
            std::cout << "east 1" << std::endl;
        }
    }
    if (CoordCheck(3, 2) && EastRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "east" << std::endl;
            EastRefine = false;
        }
        else {

            std::cout << "east 2" << std::endl;

        }
    }


    /// West Check
    if (CoordCheck(0, 1) && WestRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "west" << std::endl;
            WestRefine = false;
        }
        else {
            std::cout << "west 1" << std::endl;
        }
    }
    if (CoordCheck(2, 3) && WestRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "west" << std::endl;
            WestRefine = false;
        }
        else {
            std::cout << "west 2" << std::endl;
        }
    }




    /// North Check
    if (CoordCheck(2, 0) && NorthRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "north" << std::endl;
            NorthRefine = false;
        }
        else {
            std::cout << "north 1" << std::endl;
        }
    }
    if (CoordCheck(3, 1) && NorthRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "north" << std::endl;
            NorthRefine = false;
        }
        else {
            std::cout << "north 2" << std::endl;
        }
    }



    /// South Check
    if (CoordCheck(0, 2) && SouthRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "south" << std::endl;
            SouthRefine = false;
        }
        else {
            std::cout << "south 1" << std::endl;
        }
    }
    if (CoordCheck(1, 3) && SouthRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "south" << std::endl;
            SouthRefine = false;
        }
        else {
            std::cout << "south 2" << std::endl;
        }
    }



#if(DIM==3)
    /// East Check
    if (CoordCheck(7, 6) && EastRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "east" << std::endl;
            EastRefine = false;

        }
        else {
            std::cout << "east 3" << std::endl;
        }
    }
    if (CoordCheck(5, 4) && EastRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "east" << std::endl;
            EastRefine = false;
        }
        else {

            std::cout << "east 4" << std::endl;

        }
    }


    /// West Check
    if (CoordCheck(4, 5) && WestRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "west" << std::endl;
            WestRefine = false;
        }
        else {
            std::cout << "west 3" << std::endl;
        }
    }
    if (CoordCheck(6, 7) && WestRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "west" << std::endl;
            WestRefine = false;
        }
        else {
            std::cout << "west 4" << std::endl;
        }
    }




    /// North Check
    if (CoordCheck(6, 4) && NorthRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "north" << std::endl;
            NorthRefine = false;
        }
        else {
            std::cout << "north 3" << std::endl;
        }
    }
    if (CoordCheck(7, 5) && NorthRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "north" << std::endl;
            NorthRefine = false;
        }
        else {
            std::cout << "north 4" << std::endl;
        }
    }



    /// South Check
    if (CoordCheck(4, 6) && SouthRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "south" << std::endl;
            SouthRefine = false;
        }
        else {
            std::cout << "south 3" << std::endl;
        }
    }
    if (CoordCheck(5, 7) && SouthRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "south" << std::endl;
            SouthRefine = false;
        }
        else {
            std::cout << "south 4" << std::endl;
        }
    }


    /// Front Check
    if (CoordCheck(2, 6) && FrontRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "front" << std::endl;
            FrontRefine = false;
        }
        else {
            std::cout << "front 1" << std::endl;
        }
    }
    if (CoordCheck(3, 7) && FrontRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "front" << std::endl;
            FrontRefine = false;
        }
        else {
            std::cout << "front 2" << std::endl;
        }
    }

    if (CoordCheck(0, 4) && FrontRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "front" << std::endl;
            FrontRefine = false;
        }
        else {
            std::cout << "front 3" << std::endl;
        }
    }
    if (CoordCheck(1, 5) && FrontRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "front" << std::endl;
            FrontRefine = false;
        }
        else {
            std::cout << "front 4" << std::endl;
        }
    }



    /// Back Check
    if (CoordCheck(5, 1) && BackRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "back" << std::endl;
            BackRefine = false;
        }
        else {
            std::cout << "back 1" << std::endl;
        }
    }
    if (CoordCheck(4, 0) && BackRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "back" << std::endl;
            BackRefine = false;
        }
        else {
            std::cout << "back 2" << std::endl;
        }
    }

    if (CoordCheck(7, 3) && BackRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "back" << std::endl;
            BackRefine = false;
        }
        else {
            std::cout << "back 3" << std::endl;
        }
    }
    if (CoordCheck(6, 2) && BackRefine){

        if (m_level_target >= m_level_neigh){
            std::cout << "back" << std::endl;
            BackRefine = false;
        }
        else {
            std::cout << "back 4" << std::endl;
        }
    }













#endif


}







void NodePrint::checkNeighborsCoord(const ot::RankI *target_nodeIdsFlat, const ot::RankI *neigh_nodeIdsFlat, DENDRITE_UINT m_level_target, DENDRITE_UINT m_level_neigh){

//                std::cout << "0: " << target_nodeIdsFlat[0] << "    1: " << target_nodeIdsFlat[1] << "      2: " << target_nodeIdsFlat[2] <<
//             "    3: " << target_nodeIdsFlat[3] << std::endl;

    // East Check
    bool EastRefine = true;
    if (target_nodeIdsFlat[1] == neigh_nodeIdsFlat[0] && EastRefine){


        std::cout << "target::   0: " << target_nodeIdsFlat[0] << "    1: " << target_nodeIdsFlat[1] << "      2: "
                  << target_nodeIdsFlat[2] <<
                  "    3: " << target_nodeIdsFlat[3] << "  lvl:  " << m_level_target << std::endl;


        std::cout << "neigh::  0: " << neigh_nodeIdsFlat[0] << "    1: " << neigh_nodeIdsFlat[1] << "      2: "
                  << neigh_nodeIdsFlat[2] <<
                  "    3: " << neigh_nodeIdsFlat[3] << "  lvl:  " << m_level_neigh << std::endl;

        if (m_level_target >= m_level_neigh){
            std::cout << "east" << std::endl;
            EastRefine = false;
        }
        else {
            std::cout << "east 1" << std::endl;
        }
    }
    if (target_nodeIdsFlat[3] == neigh_nodeIdsFlat[2] && EastRefine){
        if (m_level_target >= m_level_neigh){

            std::cout << "target::   0: " << target_nodeIdsFlat[0] << "    1: " << target_nodeIdsFlat[1] << "      2: "
                      << target_nodeIdsFlat[2] <<
                      "    3: " << target_nodeIdsFlat[3] << "  lvl:  " << m_level_target << std::endl;


            std::cout << "neigh::  0: " << neigh_nodeIdsFlat[0] << "    1: " << neigh_nodeIdsFlat[1] << "      2: "
                      << neigh_nodeIdsFlat[2] <<
                      "    3: " << neigh_nodeIdsFlat[3] << "  lvl:  " << m_level_neigh << std::endl;


            std::cout << "east" << std::endl;
            EastRefine = false;
        }
        else {
            std::cout << "target::   0: " << target_nodeIdsFlat[0] << "    1: " << target_nodeIdsFlat[1] << "      2: "
                      << target_nodeIdsFlat[2] <<
                      "    3: " << target_nodeIdsFlat[3] << "  lvl:  " << m_level_target << std::endl;


            std::cout << "neigh::  0: " << neigh_nodeIdsFlat[0] << "    1: " << neigh_nodeIdsFlat[1] << "      2: "
                      << neigh_nodeIdsFlat[2] <<
                      "    3: " << neigh_nodeIdsFlat[3] << "  lvl:  " << m_level_neigh << std::endl;

            std::cout << "east 2" << std::endl;
        }
    }


}